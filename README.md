# GIS2RHESSys 1mUrban

We assume that 1-m LULC and elevation data are available for the catchment.

updated feature:
1) workflows_xxx.sh is the main script users need to edit and review. The script is helping users on setup GRASS GIS and RHESSys model using standardized procedures. However, users are responsible to understand the script and procedure and make customization for specific model setup.
2) Catchment outlet can be set by Lat/Long, UTM coordinates, and shapefile.
3) linearElevationEdit.R is a small R script to edit elevation data, e.g., bridge removal. Users need to find the editing area/location/line. Details on this will be posted later.
4) zone_cluster.R is another R script to set up zone object in RHESSys. It analyzes the slope, aspect, elevation, and hillslopes to form clusters, and these clusters are the suggested zone areas.


