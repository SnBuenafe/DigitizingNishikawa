# DigitizingNishikawa
Repository for Digitizing Nishikawa et al (1985) Data

This repository digitizes average spawning ground data on the following species from 1956-1981:
1. Yellowfin tuna
2. Swordfish
3. Albacore
4. Skipjack tuna
5. Blujefin tuna
6. Southern bluefin tuna
7. Bigeye tuna
8. Blue marlin
9. Black marlin
10. Striped marlin
11. Sailfish
12. Shortbill spearfish
13. Frigate tuna
14. Little tuna
15. Bonitos
16. Slender tuna
17. Longfin escolar
18. Sauries

Tows were done seasonally (Jan-Mar, Apr-Jun, Jul-Sept, Oct-Dec).

Number of tows and volumes of water strained are reported and digitized.

The digitized data of all 18 species in text delimited (CSV), raster, and vector files are found in `/Output`.

To replicate the process of generating these files go through the R scripts sequentially:

* `01_Digitizing_Data_CSV.R`: Converts all .csv files from `/RawData` to one long csv file for larval distributions and one for towing effort
* `02_Digitizing_Data_VectorAndRaster.R`: Creates the vector (sf objects saved as .RDS) and raster (GEOTIFF files saved as .tif) files
* `03_Making_Maps.R`: Replicates the species plots shown in the data paper.

Save FAO data in `Data/` folder to use this code. FAO data can be downloaded from:

* CWP 1x1 grid code: ["Areal Grid System"](https://data.apps.fao.org/map/catalog/srv/eng/catalog.search#/metadata/cwp-grid-map-1deg_x_1deg?tab=search)
* FAO Major Fishing Area shapefile: ["Major Fishing Areas"](https://www.fao.org/fishery/geonetwork/srv/eng/catalog.search#/metadata/ac02a460-da52-11dc-9d70-0017f293bd28)