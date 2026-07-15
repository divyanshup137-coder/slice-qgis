from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsProcessingContext,
    QgsProcessingFeedback,
    QgsProcessing,
    QgsFeature,
    QgsFields,
    QgsField,
    QgsGeometry,
    QgsWkbTypes,
    QgsFeatureSink,
    QgsPointXY,
    QgsRaster,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsSpatialIndex
)
from qgis import processing
import math
from collections import defaultdict, deque


class SLIndexAlgorithm(QgsProcessingAlgorithm):
    """SL index from a DEM by consecutive contour pairs along extracted channels.

    SL = (DH / DL) * L, evaluated between successive contour crossings.

    Channels are extracted headwater-to-outlet and processed longest-first. The
    longest channel keeps all of its SL points and claims its segments. Each
    later channel generates points only on segments not already claimed by a
    longer channel, and points whose contour pair straddles the confluence are
    dropped so that no point spans a junction. The result is a full network in
    which every channel contributes only its unique tributary reach, with no
    duplicated points on shared downstream segments.

    TRUNK_STREAMS:
        0  process all channels (full network, longest to shortest)
        N  process only the N longest independent channels

    A channel sharing more than 50 percent of its segments with an
    already-selected trunk is skipped so that near-duplicate channels do not
    consume a trunk slot.
    """

    DEM            = "DEM"
    THRESH         = "THRESH"
    INTERVAL       = "INTERVAL"
    TRUNK_STREAMS  = "TRUNK_STREAMS"

    OUTPUT_STREAMS  = "OUTPUT_STREAMS"
    OUTPUT_CONTOURS = "OUTPUT_CONTOURS"
    OUTPUT_SL_FINAL = "OUTPUT_SL_FINAL"

    # coordinate snapping precision (decimal places)
    PRECISION = 6

    # QGIS processing framework interface

    def tr(self, s):
        return QCoreApplication.translate("SLIndex", s)

    def createInstance(self):
        return SLIndexAlgorithm()

    def name(self):
        return "sl_index"

    def displayName(self):
        return self.tr("SL index (SLICE)")

    def group(self):
        return self.tr("Tectonic geomorphology")

    def groupId(self):
        return "tectonicgeomorphology"

    def shortDescription(self):
        return self.tr("Stream Length-gradient (SL) index from a DEM (SLICE).")

    def shortHelpString(self):
        return self.tr(
            "SLICE computes the Stream Length-gradient (SL) index along "
            "channels extracted from a DEM.\n\n"
            "SL = (DH / DL) * L, measured between successive contour "
            "crossings on each channel.\n\n"
            "Channels are processed longest-first: the longest keeps all of "
            "its SL points, and each shorter channel contributes points only "
            "on its own tributary reach, so shared downstream segments are "
            "not counted twice and no point spans a confluence.\n\n"
            "Parameters\n"
            "  Input DEM: elevation raster in a projected CRS (metres).\n"
            "  Flow accumulation threshold: minimum contributing cells for a "
            "cell to count as channel. Higher values give a sparser network.\n"
            "  Contour interval: vertical spacing (m) of the contours that "
            "define the DH steps.\n"
            "  Number of trunk streams: 0 builds the full network; N keeps "
            "only the N longest independent channels.\n\n"
            "Outputs\n"
            "  Stream network: channel segments split at junctions.\n"
            "  Contours: generated contours, closed loops filtered per stream.\n"
            "  SL index points: one point per contour pair, with DH, DL, L "
            "and SL.\n\n"
            "Requires the GRASS provider (r.watershed, r.stream.extract) and "
            "GDAL (gdal:contour). GeoPackage output is recommended."
        )

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.DEM, self.tr("Input DEM (projected CRS)")))

        self.addParameter(QgsProcessingParameterNumber(
            self.THRESH,
            self.tr("Flow accumulation threshold (cells)"),
            QgsProcessingParameterNumber.Integer, defaultValue=500))

        self.addParameter(QgsProcessingParameterNumber(
            self.INTERVAL,
            self.tr("Contour interval (m)"),
            QgsProcessingParameterNumber.Double, defaultValue=50.0))

        self.addParameter(QgsProcessingParameterNumber(
            self.TRUNK_STREAMS,
            self.tr("Number of trunk streams (0 = full network, "
                     "N = N longest independent streams)"),
            QgsProcessingParameterNumber.Integer,
            defaultValue=0, minValue=0))

        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_STREAMS,
            self.tr("Stream network (split segments)"),
            QgsProcessing.TypeVectorLine))

        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_CONTOURS,
            self.tr("Contours (closed-loop filtered)"),
            QgsProcessing.TypeVectorLine))

        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_SL_FINAL,
            self.tr("SL index points"),
            QgsProcessing.TypeVectorPoint))

    # Geometry helpers

    def snap_coordinate(self, coord, precision=None):
        if precision is None:
            precision = self.PRECISION
        return round(coord, precision)

    def snap_point(self, point, precision=None):
        if precision is None:
            precision = self.PRECISION
        return (self.snap_coordinate(point.x(), precision),
                self.snap_coordinate(point.y(), precision))

    def is_line_geometry(self, geom):
        if geom is None or geom.isEmpty():
            return False
        return geom.wkbType() in [
            QgsWkbTypes.LineString, QgsWkbTypes.LineString25D,
            QgsWkbTypes.MultiLineString, QgsWkbTypes.MultiLineString25D
        ]

    def extract_line_features(self, layer, feedback):
        return [f for f in layer.getFeatures()
                if f.geometry() and not f.geometry().isEmpty()
                and self.is_line_geometry(f.geometry())]

    def is_closed_contour(self, geom):
        if geom is None or geom.isEmpty():
            return False
        if geom.isMultipart():
            parts = geom.asMultiPolyline()
            for part in parts:
                if len(part) < 2:
                    return False
                if (abs(part[0].x() - part[-1].x()) > 1e-6 or
                        abs(part[0].y() - part[-1].y()) > 1e-6):
                    return False
            return len(parts) > 0
        line = geom.asPolyline()
        if len(line) < 2:
            return False
        return (abs(line[0].x() - line[-1].x()) < 1e-6 and
                abs(line[0].y() - line[-1].y()) < 1e-6)

    # Network construction

    def manual_split_at_junctions(self, stream_features, feedback):
        feedback.pushInfo("  Splitting at junctions...")
        precision = self.PRECISION
        endpoint_count = defaultdict(int)
        line_data = []

        for idx, feat in enumerate(stream_features):
            geom = feat.geometry()
            if not self.is_line_geometry(geom):
                continue
            if geom.isMultipart():
                for pidx, part in enumerate(geom.asMultiPolyline()):
                    if len(part) >= 2:
                        s = self.snap_point(part[0], precision)
                        e = self.snap_point(part[-1], precision)
                        endpoint_count[s] += 1
                        endpoint_count[e] += 1
                        line_data.append((idx, pidx, part, s, e, feat))
            else:
                line = geom.asPolyline()
                if len(line) >= 2:
                    s = self.snap_point(line[0], precision)
                    e = self.snap_point(line[-1], precision)
                    endpoint_count[s] += 1
                    endpoint_count[e] += 1
                    line_data.append((idx, None, line, s, e, feat))

        junctions = {pt for pt, cnt in endpoint_count.items() if cnt >= 3}
        feedback.pushInfo(f"  Junctions: {len(junctions)}")

        split_features = []
        for idx, pidx, line, start, end, feat in line_data:
            jpos = []
            for i, v in enumerate(line):
                if i in (0, len(line) - 1):
                    continue
                if self.snap_point(v, precision) in junctions:
                    jpos.append(
                        (QgsGeometry.fromPolylineXY(line[:i+1]).length(), i, v))

            if not jpos:
                split_features.append(feat)
                continue

            jpos.sort(key=lambda x: x[0])
            prev = 0
            for _, ipos, jv in jpos:
                part = line[prev:ipos+1]
                if len(part) >= 2:
                    part[-1] = jv
                    nf = QgsFeature(feat.fields())
                    nf.setGeometry(QgsGeometry.fromPolylineXY(part))
                    nf.setAttributes(feat.attributes())
                    split_features.append(nf)
                prev = ipos

            tail = line[prev:]
            if len(tail) >= 2:
                nf = QgsFeature(feat.fields())
                nf.setGeometry(QgsGeometry.fromPolylineXY(tail))
                nf.setAttributes(feat.attributes())
                split_features.append(nf)

        feedback.pushInfo(f"  Split segments: {len(split_features)}")
        return split_features

    def build_graph_from_streams(self, stream_features, feedback):
        graph = defaultdict(list)
        reverse_graph = defaultdict(list)
        node_coords = {}
        segment_data = {}
        precision = self.PRECISION

        for idx, feat in enumerate(stream_features):
            geom = feat.geometry()
            if not self.is_line_geometry(geom):
                continue
            if geom.isMultipart():
                for pidx, part in enumerate(geom.asMultiPolyline()):
                    if len(part) < 2:
                        continue
                    pg = QgsGeometry.fromPolylineXY(part)
                    s = self.snap_point(part[0], precision)
                    e = self.snap_point(part[-1], precision)
                    node_coords[s] = part[0]
                    node_coords[e] = part[-1]
                    sid = f"seg_{idx}_p{pidx}"
                    segment_data[sid] = {'geometry': pg, 'start': s,
                                         'end': e, 'length': pg.length()}
                    graph[s].append((e, sid))
                    reverse_graph[e].append(s)
                continue
            line = geom.asPolyline()
            if len(line) < 2:
                continue
            s = self.snap_point(line[0], precision)
            e = self.snap_point(line[-1], precision)
            node_coords[s] = line[0]
            node_coords[e] = line[-1]
            sid = f"seg_{idx}"
            segment_data[sid] = {'geometry': geom, 'start': s,
                                 'end': e, 'length': geom.length()}
            graph[s].append((e, sid))
            reverse_graph[e].append(s)

        return graph, reverse_graph, node_coords, segment_data

    def find_headwaters(self, graph, reverse_graph):
        all_nodes = set(graph.keys()) | set(reverse_graph.keys())
        return [n for n in all_nodes
                if (n not in reverse_graph or not reverse_graph[n])
                and n in graph]

    def compute_strahler_order(self, graph, reverse_graph, headwaters,
                               segment_data, feedback):
        segment_order = {}
        node_order = {hw: 1 for hw in headwaters}
        upstream_count = {n: len(reverse_graph[n]) for n in reverse_graph}
        queue = deque(headwaters)
        processed = set()

        while queue:
            cur = queue.popleft()
            order = node_order.get(cur, 1)
            if cur not in graph:
                continue
            for nxt, sid in graph[cur]:
                if sid not in processed:
                    segment_order[sid] = order
                    processed.add(sid)
                    upstream_count[nxt] = upstream_count.get(nxt, 0) - 1
                    if upstream_count[nxt] == 0:
                        up_orders = [
                            segment_order[s]
                            for u in reverse_graph.get(nxt, [])
                            for (d, s) in graph.get(u, [])
                            if d == nxt and s in segment_order
                        ]
                        if not up_orders:
                            no = 1
                        elif len(up_orders) == 1:
                            no = up_orders[0]
                        else:
                            mx = max(up_orders)
                            no = mx + 1 if up_orders.count(mx) >= 2 else mx
                        node_order[nxt] = no
                        queue.append(nxt)

        return segment_order, max(segment_order.values()) if segment_order else 0

    def compute_river_ids(self, graph, reverse_graph, segment_data, feedback):
        all_nodes = set(graph.keys()) | set(reverse_graph.keys())
        outlet_nodes = [n for n in all_nodes
                        if n not in graph or not graph[n]]
        feedback.pushInfo(f"  Outlet nodes found: {len(outlet_nodes)}")

        node_river = {}
        river_id = 0
        for outlet in outlet_nodes:
            if outlet in node_river:
                continue
            queue = deque([outlet])
            visited = set()
            while queue:
                node = queue.popleft()
                if node in visited:
                    continue
                visited.add(node)
                node_river[node] = river_id
                for up in reverse_graph.get(node, []):
                    if up not in visited:
                        queue.append(up)
            river_id += 1

        unassigned = 0
        for sid, sinfo in segment_data.items():
            rid = node_river.get(sinfo['end'])
            if rid is None:
                rid = node_river.get(sinfo['start'], 0)
                unassigned += 1
            segment_data[sid]['river_id'] = rid

        feedback.pushInfo(f"  River basins identified: {river_id}")
        if unassigned:
            feedback.pushInfo(
                f"  Warning: {unassigned} segments used start-node fallback")
        return river_id

    # Channel extraction

    def extract_channels_from_headwaters(self, headwaters, graph,
                                         segment_data, feedback):
        feedback.pushInfo("\n  Extracting channels from headwaters...")
        channels = []
        channel_id = 0

        for headwater in headwaters:
            current_node = headwater
            channel_segments = []
            visited_nodes = set()

            while current_node in graph:
                if current_node in visited_nodes:
                    break
                visited_nodes.add(current_node)
                edges = graph[current_node]
                if not edges:
                    break
                next_node, sid = edges[0]
                channel_segments.append(sid)
                current_node = next_node

            if not channel_segments:
                continue

            merged_vertices = []
            for sid in channel_segments:
                seg_line = segment_data[sid]['geometry'].asPolyline()
                merged_vertices.extend(
                    seg_line if not merged_vertices else seg_line[1:])

            if len(merged_vertices) < 2:
                continue

            merged_geom = QgsGeometry.fromPolylineXY(merged_vertices)
            channels.append({
                'channel_id': channel_id,
                'path_id': f"HW{channel_id:03d}",
                'geometry': merged_geom,
                'length': merged_geom.length(),
                'segment_ids': channel_segments,
                'headwater': headwater
            })
            channel_id += 1

        feedback.pushInfo(f"  Channels extracted: {len(channels)}")
        return channels

    def find_segment_for_midpoint(self, mid_distance, channel_segments,
                                  segment_data):
        accumulated = 0.0
        for sid in channel_segments:
            accumulated += segment_data[sid]['length']
            if mid_distance <= accumulated:
                return sid
        return channel_segments[-1]

    # SL point generation along a channel

    def generate_sl_for_channel(self, channel, contour_index, contour_dict,
                                contour_is_closed, dem_provider,
                                segment_data, feedback):
        merged_geom = channel['geometry']
        channel_id  = channel['channel_id']
        path_id     = channel['path_id']

        cand_ids      = contour_index.intersects(merged_geom.boundingBox())
        intersections = []

        for cid in cand_ids:
            cf   = contour_dict[cid]
            cg   = cf.geometry()
            elev = cf['ELEV']
            is_closed = cid in contour_is_closed
            if merged_geom.intersects(cg):
                ig = merged_geom.intersection(cg)
                if ig.isEmpty():
                    continue
                if ig.type() == QgsWkbTypes.PointGeometry:
                    pts = (ig.asMultiPoint() if ig.isMultipart()
                           else [ig.asPoint()])

                    # If closed contour produces 2+ crossings with this
                    # stream (entry + exit through the loop), skip ALL
                    # of them -- this is always a DEM artifact pattern
                    if is_closed and len(pts) >= 2:
                        continue

                    for pt in pts:
                        intersections.append({
                            'point':     pt,
                            'elev':      elev,
                            'geom':      QgsGeometry.fromPointXY(pt),
                            'is_closed': is_closed,
                        })

        if len(intersections) < 2:
            return []

        for ix in intersections:
            ix['dist'] = merged_geom.lineLocatePoint(ix['geom'])
        intersections.sort(key=lambda x: x['dist'])

        # --- Per-stream closed contour filter ---
        # Rule: if a closed contour crossing sits between two open
        # contour crossings (has open upstream AND open downstream),
        # it is redundant / a DEM artifact --> skip it.
        # If it is at the edge (no open contour upstream, i.e. near
        # the headwater on a mountain peak), it is valid --> keep it.
        has_open_before = False
        # First pass: mark whether each intersection has an open one
        # somewhere before it (upstream)
        open_before = []
        for ix in intersections:
            open_before.append(has_open_before)
            if not ix['is_closed']:
                has_open_before = True

        has_open_after = False
        # Second pass (reverse): mark whether each intersection has
        # an open one somewhere after it (downstream)
        open_after = [False] * len(intersections)
        for i in range(len(intersections) - 1, -1, -1):
            open_after[i] = has_open_after
            if not intersections[i]['is_closed']:
                has_open_after = True

        filtered = []
        for i, ix in enumerate(intersections):
            if ix['is_closed'] and open_before[i] and open_after[i]:
                # Closed contour between open contours --> skip
                continue
            filtered.append(ix)

        intersections = filtered

        if len(intersections) < 2:
            return []

        sl_points = []
        for j in range(len(intersections) - 1):
            i1 = intersections[j]
            i2 = intersections[j + 1]

            DH = abs(i2['elev'] - i1['elev'])
            if DH < 1.0:
                continue
            DL = i2['dist'] - i1['dist']
            if DL < 1.0:
                continue

            mid_d    = (i1['dist'] + i2['dist']) / 2.0
            mid_geom = merged_geom.interpolate(mid_d)
            if mid_geom.isEmpty():
                continue

            mid_pt = mid_geom.asPoint()
            sid    = self.find_segment_for_midpoint(
                mid_d, channel['segment_ids'], segment_data)
            river_id = segment_data[sid].get('river_id', 0)

            ident = dem_provider.identify(mid_pt, QgsRaster.IdentifyFormatValue)
            z_mid = (ident.results()[1]
                     if ident.isValid() and 1 in ident.results()
                     else (i1['elev'] + i2['elev']) / 2.0)

            L  = mid_d
            SL = (DH / DL) * L

            if math.isfinite(SL) and SL >= 0:
                sl_points.append({
                    'point':      mid_pt,
                    'path_id':    path_id,
                    'channel_id': channel_id,
                    'segment_id': sid,
                    'river_id':   river_id,
                    'dist':       mid_d,
                    'x':          mid_pt.x(),
                    'y':          mid_pt.y(),
                    'z':          z_mid,
                    'dh':         DH,
                    'dl':         DL,
                    'l':          L,
                    'sl':         SL,
                    'elev_lower': min(i1['elev'], i2['elev']),
                    'elev_upper': max(i1['elev'], i2['elev']),
                })

        return sl_points

    # Progressive network builder

    def build_network_progressively(self, channels, n_trunks,
                                     contour_index, contour_dict,
                                     contour_is_closed,
                                     dem_provider, segment_data,
                                     feedback):
        """
        Core algorithm: process channels from longest to shortest.

        Trunk 1 (longest): all SL points kept, segments claimed.
        Trunk N (N>1): SL points generated, then:
          - Points on already-claimed segments are removed
            (those belong to a longer channel)
          - Points whose contour pair straddles the confluence
            (downstream crossing exceeds the confluence distance)
            are removed to avoid duplicate points at junctions
          - Points sitting entirely within the tributary are kept
          - Remaining segments are claimed

        For TRUNK_STREAMS=0: processes all channels.
        For TRUNK_STREAMS=N: processes N channels, skipping
        near-twins (>50% segment overlap with claimed set).
        """
        sorted_channels = sorted(channels, key=lambda c: c['length'],
                                 reverse=True)

        # Determine how many to process
        process_all = (n_trunks == 0)
        max_trunks = len(sorted_channels) if process_all else n_trunks

        claimed_segments = set()      # segments owned by earlier trunks
        all_accepted_points = []      # final output
        trunk_number = 0              # counter for accepted trunks
        skipped_twins = 0

        for ch_idx, channel in enumerate(sorted_channels):
            if trunk_number >= max_trunks:
                break

            if feedback.isCanceled():
                break

            ch_segs = set(channel['segment_ids'])

            # --- Independence check (skip near-twins) ---
            # For trunk mode (n_trunks > 0), skip channels with >50% overlap.
            # For ALL mode, we still skip channels with 100% overlap
            # (pure duplicates that add nothing).
            if claimed_segments:
                overlap = ch_segs & claimed_segments
                overlap_ratio = len(overlap) / len(ch_segs) if ch_segs else 1.0

                if not process_all and overlap_ratio >= 0.50:
                    skipped_twins += 1
                    if trunk_number < 10 or skipped_twins <= 5:
                        feedback.pushInfo(
                            f"  Skipped: {channel['path_id']} "
                            f"(length={channel['length']:.1f} m, "
                            f"overlap={overlap_ratio*100:.0f}% -- near-twin)")
                    continue

                if process_all and overlap_ratio >= 1.0:
                    # Pure duplicate, every segment already claimed
                    skipped_twins += 1
                    continue

            trunk_number += 1
            trunk_label = f"Trunk{trunk_number}"

            # Reassign path_id to trunk label
            channel['path_id'] = trunk_label

            # --- Generate SL points for this channel ---
            raw_pts = self.generate_sl_for_channel(
                channel, contour_index, contour_dict,
                contour_is_closed, dem_provider, segment_data, feedback)

            if trunk_number == 1:
                # First trunk: keep ALL points
                accepted = raw_pts
                blocked_shared = 0
                dropped_confluence = 0
            else:
                # Subsequent trunks: filter out points on claimed segments
                unclaimed_pts = []
                blocked_shared = 0
                for pt in raw_pts:
                    if pt['segment_id'] in claimed_segments:
                        blocked_shared += 1
                    else:
                        unclaimed_pts.append(pt)

                # --- Confluence-aware filtering ---
                # Find the confluence distance: walk through the channel's
                # segments from headwater to outlet, accumulating distance
                # until we hit the first claimed segment.  That boundary
                # is where this tributary joins the already-built network.
                confluence_dist = 0.0
                for sid in channel['segment_ids']:
                    if sid in claimed_segments:
                        break
                    confluence_dist += segment_data[sid]['length']

                # For each unclaimed point, check whether its contour pair
                # straddles the confluence.  The downstream contour crossing
                # is approximately at (dist + dl/2).  If that crossing is
                # beyond the confluence distance, the point spans the
                # junction and should be dropped.  If both crossings are
                # upstream of the confluence, the point sits entirely
                # within the tributary and is valid.
                dropped_confluence = 0
                accepted_unclaimed = []
                for pt in unclaimed_pts:
                    downstream_crossing = pt['dist'] + pt['dl'] / 2.0
                    if downstream_crossing > confluence_dist:
                        dropped_confluence += 1
                    else:
                        accepted_unclaimed.append(pt)

                accepted = accepted_unclaimed

            # Update path_id for all accepted points
            for pt in accepted:
                pt['path_id'] = trunk_label
                pt['trunk_number'] = trunk_number

            all_accepted_points.extend(accepted)

            # Claim this channel's segments
            claimed_segments.update(ch_segs)

            # Log progress
            if trunk_number <= 10 or trunk_number % 50 == 0 or \
               trunk_number == max_trunks:
                feedback.pushInfo(
                    f"  {trunk_label} ({channel['path_id']}): "
                    f"length={channel['length']:.1f} m, "
                    f"raw={len(raw_pts)}, "
                    f"accepted={len(accepted)}, "
                    f"blocked_shared={blocked_shared}, "
                    f"dropped_confluence={dropped_confluence}")

        feedback.pushInfo(f"\n  Trunks processed  : {trunk_number}")
        if skipped_twins:
            feedback.pushInfo(f"  Skipped (twins)   : {skipped_twins}")
        feedback.pushInfo(
            f"  Total SL points   : {len(all_accepted_points)}")

        return all_accepted_points, trunk_number

    # Main algorithm

    def processAlgorithm(self, parameters, context: QgsProcessingContext,
                         feedback: QgsProcessingFeedback):

        dem = self.parameterAsRasterLayer(
            parameters, self.DEM, context)
        thresh = self.parameterAsInt(
            parameters, self.THRESH, context)
        elev_interval = self.parameterAsDouble(
            parameters, self.INTERVAL, context)
        trunk_streams = self.parameterAsInt(
            parameters, self.TRUNK_STREAMS, context)

        crs = dem.crs()
        cell_size = dem.rasterUnitsPerPixelX()

        feedback.pushInfo("SLICE: SL index by progressive longest-first network builder")
        feedback.pushInfo(f"DEM cell size: {cell_size:.2f} m")
        if trunk_streams > 0:
            feedback.pushInfo(
                f"Trunk streams: {trunk_streams} longest independent")
        else:
            feedback.pushInfo("Trunk streams: ALL (full network)")

        # -- STEP 1: Stream network -------------------------------------------
        feedback.pushInfo("\n[STEP 1] Network Generation and Topology")

        acc = processing.run(
            "grass:r.watershed",
            {"elevation": dem, "threshold": thresh,
             "accumulation": "TEMPORARY_OUTPUT", "-m": True},
            context=context, feedback=feedback
        )["accumulation"]

        stream_path = processing.run(
            "grass:r.stream.extract",
            {"elevation": dem, "accumulation": acc,
             "threshold": thresh, "stream_vector": "TEMPORARY_OUTPUT"},
            context=context, feedback=feedback
        )["stream_vector"]

        temp_streams = QgsVectorLayer(stream_path, "temp_streams", "ogr")
        if not temp_streams.isValid():
            raise QgsProcessingException("Failed to load stream network")

        stream_features = self.extract_line_features(temp_streams, feedback)
        split_features = self.manual_split_at_junctions(
            stream_features, feedback)

        graph, reverse_graph, node_coords, segment_data = \
            self.build_graph_from_streams(split_features, feedback)
        headwaters = self.find_headwaters(graph, reverse_graph)
        feedback.pushInfo(f"  Headwaters: {len(headwaters)}")

        segment_order, max_order = self.compute_strahler_order(
            graph, reverse_graph, headwaters, segment_data, feedback)

        order_counts = defaultdict(int)
        for o in segment_order.values():
            order_counts[o] += 1
        feedback.pushInfo("\n  Strahler order distribution:")
        for o in sorted(order_counts):
            feedback.pushInfo(f"    Order {o}: {order_counts[o]} segments")
        feedback.pushInfo(f"    Maximum order: {max_order}")

        # -- STEP 1b: River basin IDs ------------------------------------------
        feedback.pushInfo("\n[STEP 1b] River Basin ID Assignment")
        n_basins = self.compute_river_ids(
            graph, reverse_graph, segment_data, feedback)

        # -- STEP 2: Contour generation ----------------------------------------
        feedback.pushInfo("\n[STEP 2] Contour Generation (closed-loop filtered)")

        contours_path = processing.run(
            "gdal:contour",
            {"INPUT": dem, "INTERVAL": elev_interval,
             "FIELD_NAME": "ELEV", "OUTPUT": "TEMPORARY_OUTPUT"},
            context=context, feedback=feedback
        )["OUTPUT"]

        contours_layer = QgsVectorLayer(contours_path, "contours", "ogr")
        contour_index = QgsSpatialIndex()
        contour_dict = {}       # ALL contours -- used for SL
        contour_is_closed = {}  # track which contour IDs are closed
        all_contours = []       # ALL contours -- used for output layer
        n_total = 0
        n_closed = 0

        for contour in contours_layer.getFeatures():
            n_total += 1
            all_contours.append(contour)
            contour_index.addFeature(contour)
            contour_dict[contour.id()] = contour
            if self.is_closed_contour(contour.geometry()):
                contour_is_closed[contour.id()] = True
                n_closed += 1

        feedback.pushInfo(f"  Total contours     : {n_total}")
        feedback.pushInfo(f"  Closed contours    : {n_closed}")
        feedback.pushInfo(f"  Open contours      : {n_total - n_closed}")
        feedback.pushInfo(f"  Closed filtered per-stream: between open contours = skip, edge = keep")

        # -- STEP 3: Channel extraction ----------------------------------------
        feedback.pushInfo("\n[STEP 3] Channel Extraction (Headwater-Based)")

        channels = self.extract_channels_from_headwaters(
            headwaters, graph, segment_data, feedback)
        channels.sort(key=lambda c: c['length'], reverse=True)

        if channels:
            feedback.pushInfo(f"  Longest : {channels[0]['length']:.1f} m")
            feedback.pushInfo(f"  Shortest: {channels[-1]['length']:.1f} m")
            feedback.pushInfo(f"  Total   : {len(channels)}")

        # -- STEP 4: Progressive network build ---------------------------------
        feedback.pushInfo(
            "\n[STEP 4] Progressive Network Build (longest-first)")

        dem_provider = dem.dataProvider()

        all_sl_data, n_trunks_processed = self.build_network_progressively(
            channels, trunk_streams,
            contour_index, contour_dict,
            contour_is_closed,
            dem_provider, segment_data,
            feedback)

        # -- OUTPUTS -----------------------------------------------------------
        feedback.pushInfo("\n[OUTPUTS] Writing Results")

        # Stream segment layer
        stream_fields = QgsFields()
        stream_fields.append(QgsField("seg_id", QVariant.String))
        stream_fields.append(QgsField("RiverID", QVariant.Int))
        stream_fields.append(QgsField("order", QVariant.Int))
        stream_fields.append(QgsField("length", QVariant.Double))

        stream_sink, stream_dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_STREAMS, context,
            stream_fields, QgsWkbTypes.LineString, crs)

        for sid, sinfo in segment_data.items():
            feat = QgsFeature(stream_fields)
            feat.setGeometry(sinfo['geometry'])
            feat.setAttributes([
                sid,
                sinfo.get('river_id', 0),
                segment_order.get(sid, 1),
                sinfo['length'],
            ])
            stream_sink.addFeature(feat, QgsFeatureSink.FastInsert)

        # Contour layer
        contours_sink, contours_dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_CONTOURS, context,
            contours_layer.fields(), QgsWkbTypes.LineString, crs)

        for cf in all_contours:
            contours_sink.addFeature(cf, QgsFeatureSink.FastInsert)

        # SL point layer
        sl_fields = QgsFields()
        sl_fields.append(QgsField("PathID", QVariant.String))
        sl_fields.append(QgsField("BranchID", QVariant.Int))
        sl_fields.append(QgsField("StreamID", QVariant.String))
        sl_fields.append(QgsField("RiverID", QVariant.Int))
        sl_fields.append(QgsField("X", QVariant.Double))
        sl_fields.append(QgsField("Y", QVariant.Double))
        sl_fields.append(QgsField("Z", QVariant.Double))
        sl_fields.append(QgsField("DH", QVariant.Double))
        sl_fields.append(QgsField("DL", QVariant.Double))
        sl_fields.append(QgsField("L", QVariant.Double))
        sl_fields.append(QgsField("SL", QVariant.Double))

        sl_sink, sl_dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_SL_FINAL, context,
            sl_fields, QgsWkbTypes.Point, crs)

        for d in all_sl_data:
            feat = QgsFeature(sl_fields)
            feat.setGeometry(QgsGeometry.fromPointXY(d['point']))
            feat.setAttributes([
                d['path_id'],
                d['channel_id'],
                d['segment_id'],
                d['river_id'],
                d['x'], d['y'], d['z'],
                d['dh'], d['dl'], d['l'], d['sl'],
            ])
            sl_sink.addFeature(feat, QgsFeatureSink.FastInsert)

        feedback.pushInfo("\nComplete")
        feedback.pushInfo(f"  Headwaters        : {len(headwaters)}")
        feedback.pushInfo(f"  River basins      : {n_basins}")
        feedback.pushInfo(f"  Channels extracted: {len(channels)}")
        feedback.pushInfo(f"  Trunks processed  : {n_trunks_processed}")
        feedback.pushInfo(f"  Closed contours   : {n_closed} (filtered per-stream)")
        feedback.pushInfo(f"  SL points (final) : {len(all_sl_data)}")
        if trunk_streams > 0:
            feedback.pushInfo(
                f"  Mode: {trunk_streams} trunk streams (independent)")
        else:
            feedback.pushInfo("  Mode: full network (all channels)")
        feedback.pushInfo(
            "  Method: progressive longest-first network builder")

        return {
            self.OUTPUT_STREAMS:  stream_dest_id,
            self.OUTPUT_CONTOURS: contours_dest_id,
            self.OUTPUT_SL_FINAL: sl_dest_id
        }