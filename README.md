# SLICE: Stream Length-Gradient Index by Constant Elevation

An open-source QGIS Processing plugin for automated network-wide computation of the constant-elevation (ΔH) Stream Length-Gradient index.

## Description
SLICE automates the computation of the Stream Length-Gradient (SL) index across a complete drainage network directly from a single Digital Elevation Model (DEM). It utilizes the constant-elevation interval sampling method (constant ΔH), calculating the index between successive contour crossings. The index is evaluated as SL = (ΔH / ΔL) * L.

To prevent spurious data generation at downstream confluences, the plugin employs a longest-first, confluence-aware deduplication strategy. The longest channel retains all its calculated points, while shorter tributaries contribute only their unique upstream sections. 

This tool is intended for active-tectonics, tectonic-geomorphology, and knickpoint analysis.

## Requirements
*   **QGIS**: Version 3.40 LTR or later.
*   **GRASS GIS**: Must be enabled in the QGIS Processing Toolbox (requires `r.watershed` and `r.stream.extract`)
*   **GDAL**: Must be enabled for contour generation (`gdal:contour`).
*   **Input Data**: A DEM in a projected Coordinate Reference System (measured in meters).

## Installation Directory
To install the plugin manually:
1. Download this repository as a `.zip` file and extract it, or clone it via git.
2. Rename the extracted folder to `slice`.
3. Move the `slice` folder into your QGIS Python plugins directory:
   * **Windows:** `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\`
   * **Linux:** `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
   * **macOS:** `~/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins/`
4. Restart QGIS.
5. Go to **Plugins > Manage and Install Plugins...**, click on **Installed**, and check the box next to **SLICE** to enable it.

## How to Run It
1. In QGIS, open the **Processing Toolbox** (gear icon).
2. Expand the **Tectonic geomorphology** group and double-click the **SL index (SLICE)** algorithm.
3. Fill in the required parameters:
   * **Input DEM**: Select your projected DEM raster.
   * **Flow accumulation threshold (cells)**: Minimum contributing cells to define a channel (default is 500).
   * **Contour interval (m)**: Vertical spacing for the ΔH steps (default is 50.0).
   * **Number of trunk streams**: Enter `0` to build the full network, or `N` to process only the *N* longest independent channels.
4. Click **Run**.

## Outputs
The plugin outputs three layers (GeoPackage output is recommended):
*   **Stream network**: Channel segments split at junctions.
*   **Contours**: Generated contours, with closed-loop artifacts filtered out per stream.
*   **SL index points**: One point per contour pair containing calculated attributes including ΔH, ΔL, L, and SL.

## License
This project is licensed under the GNU General Public License v3.0 (GPL-3.0)[cite: 2].
