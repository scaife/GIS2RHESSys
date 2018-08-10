#!/bin/bash
#
# assume new fully bundled GRASS 7.4 is installed from http://grassmac.wikidot.com/downloads
# assume R and rgrass7 package are installed.
# assume ssurgo_extraction.R is downloaded and saved together with this script.
#
# path setup below is for Mac OS
# adding R paths
export PATH=$PATH:/usr/local/bin
export PATH=$PATH:/Library/Frameworks/R.framework/Resources
# adding GRASS bin paths
export PATH=$PATH:/Applications/GRASS-7.4.0.app/Contents/Resources/bin:/Applications/GRASS-7.4.0.app/Contents/Resources/scripts:/etc:/usr/lib
# adding GRASS environment variables: GISBASE_USER and GISBASE_SYSTEM in Grass.sh
export PATH=$PATH:~/Library/GRASS/7.4/Modules/bin:~/Library/GRASS/7.4/Modules/scripts
# adding GDAL and EPSG code lookup paths (Running gdal-config --datadir shows where GDAL searches for gcs.csv)
export GDAL_DATA=/Applications/GRASS-7.4.0.app/Contents/Resources/share/gdal
############################################################################################################
############################################################################################################
###
### procedure I - setup GRASS dataset
PROJDIR='/Users/laurencelin/Downloads/SLB' # full path to the project location;
EPSGCODE='EPSG:26918' # need to (manually) lookup the EPSG code for NAD83 UTM ##N for the catchment
RESOLUTION=3 #spatial resolution (meters) of the grids
RHESSysNAME='rhessys_SLB3m' # e.g., rhessys_baisman10m
GISDBASE=$PROJDIR/grass_dataset
LOCATION_NAME=$RHESSysNAME
LOCATION=$GISDBASE/$LOCATION_NAME
MAPSET=PERMANENT
grass74 -c $EPSGCODE -e $LOCATION
###
### procedure II - setup catchment elevation
### ... rescale source elevation data to the defined $RESOLUTION (optional;no need if source data is at the same $RESOLUTION)
LOCATIONDEM=$GISDBASE/elevationRAW
grass74 $LOCATIONDEM/$MAPSET --exec g.region raster=demRAW
grass74 $LOCATIONDEM/$MAPSET --exec g.region res=$RESOLUTION -a -p
grass74 $LOCATIONDEM/$MAPSET --exec r.resamp.stats -w input=demRAW output=dem$RESOLUTION'm'
grass74 $LOCATIONDEM/$MAPSET --exec r.out.gdal --overwrite input=dem$RESOLUTION'm' output=$PROJDIR/raw_data/dem$RESOLUTION'm.tif' format=GTiff
### ... import the (rescaled) elevation data into "$LOCATION/$MAPSET"
grass74 $LOCATION/$MAPSET --exec r.in.gdal -o -e --overwrite input=$PROJDIR/raw_data/dem$RESOLUTION'm.tif' output=dem_
grass74 $LOCATION/$MAPSET --exec g.region raster=dem_
### ... edit the DEM for bridge removal (optional)
grass74 $LOCATION/$MAPSET --exec r.mapcalc --overwrite expression="xmap = x()"
grass74 $LOCATION/$MAPSET --exec r.mapcalc --overwrite expression="ymap = y()"
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/LinearElevationEdit.R xmap ymap dem_
###
### procedure III - setup catchment soil
### ... reproject & import (reprojected) soil data to "$LOCATION/$MAPSET" and convert shape-polygon to raster soil (USDA) classes
downloadedSSURGOdirectory='MD005'
LOCATIONSOIL=$GISDBASE/soilRAW
grass74 $LOCATION/$MAPSET --exec v.proj --overwrite location=soilRAW mapset=PERMANENT input=ssurgo output=ssurgo
grass74 $LOCATION/$MAPSET --exec v.to.rast --overwrite input=ssurgo use=cat output=soil_ssurgo
grass74 $LOCATION/$MAPSET --exec v.db.select map=ssurgo separator=comma file=$PROJDIR/$RHESSysNAME/soil_cat_mukey.csv
#Rscript $PROJDIR/$RHESSysNAME/ssurgo_extraction-master/ssurgo_extraction.R $PROJDIR/raw_data/$downloadedSSURGOdirectory
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/ssurgo_extraction-master/ssurgo_soiltexture2gis.R $PROJDIR/$RHESSysNAME/soil_cat_mukey.csv $PROJDIR/raw_data/$downloadedSSURGOdirectory/soil_mukey_texture.csv
### ... use soil to define riparian are in this case (optional)
grass74 $LOCATION/$MAPSET --exec r.mapcalc --overwrite expression="riparian = if(soil_texture==8,1,null())"
###
### procedure IV - delineate watershad & basic rhessys spatial structures
### ... if input is Lat/Long
#gageLat='39.47947' # catchment outlet WSG84 Lat (decimal degree)
#gageLong='-76.67803' # catchment outlet WSG84 Long (decimal degree; includes the negative sign if applied)
#declare $(grass74 $LOCATION/$MAPSET --exec m.proj -i coordinates=$gageLong,$gageLat separator=space | awk '{print "xyCoord=" $1 "," $2}')
#echo $xyCoord | grass74 $LOCATION/$MAPSET --exec v.in.ascii in=- out=outlet x=1 y=2 separator=, --overwrite
## ... if input is UTM X Y
xyCoord='345614.98,4359826.84'
echo $xyCoord | grass74 $LOCATION/$MAPSET --exec v.in.ascii in=- out=outlet x=1 y=2 separator=, --overwrite
### ... if input is a shapefile point
#gageShapefile='usgs.shp'
#grass74 $LOCATION/$MAPSET --exec v.in.ogr --overwrite input=$PROJDIR/raw_data/$gageShapefile output=outlet location=outletRAW
#LOCATIONOUTLET=$GISDBASE/outletRAW
#grass74 $LOCATION/$MAPSET --exec v.proj --overwrite location=outletRAW mapset=PERMANENT input=outlet output=outlet
# declare $(v.out.ascii input=outlet type=point format=point separator=space | awk '{print "xyCoord=" $1 "," $2}')
thres='620000' # meter squre 62 ha
expectedDrainageArea=826191 # meter squre (allow 6% error)
grass74 $LOCATION/$MAPSET --exec sh $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/grass_delineation.sh $(bc <<< "$thres/$RESOLUTION/$RESOLUTION") $(bc <<< "0.98*$expectedDrainageArea/$RESOLUTION/$RESOLUTION") $(bc <<< "1.02*$expectedDrainageArea/$RESOLUTION/$RESOLUTION")
### ... set zone (optional)
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/zone_cluster.R dem slope aspect hill
###
### procedure V - (1-m) LULC, vectorize "patch" from "$LOCATION/$MAPSET" and calculate
grass74 $LOCATION/$MAPSET --exec r.to.vect input=patch output=patch type=area
LOCATIONLULC=$GISDBASE/lulcRAW
grass74 $LOCATIONLULC/$MAPSET --exec v.proj location=$LOCATION_NAME mapset=PERMANENT input=patch output=patch$RESOLUTION'm'
grass74 $LOCATIONLULC/$MAPSET --exec v.to.rast input=patch$RESOLUTION'm' output=patch$RESOLUTION'm' use=attr attribute_column=value
grass74 $LOCATIONLULC/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/aggregate_lulcFrac.R patch$RESOLUTION'm' lulcRAW $PROJDIR/$RHESSysNAME/lulcFrac$RESOLUTION'm.csv'
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/aggregate_lulcFrac_write2GIS.R patch $PROJDIR/$RHESSysNAME/lulcFrac$RESOLUTION'm.csv'
###
### procedure VI - import road network
LOCATIONROAD=$GISDBASE/roadRAW
grass74 $LOCATION/$MAPSET --exec v.proj --overwrite location=roadRAW mapset=PERMANENT input=roads output=roadline
### ... reduce road vector database and set buffer for roads.
grass74 $LOCATION/$MAPSET --exec r.to.vect input=basin output=basin type=area
grass74 $LOCATION/$MAPSET --exec v.clip -r input=roadline clip=basin output=cliproadline
grass74 $LOCATION/$MAPSET --exec v.buffer input=cliproadline type=line output=roadbuffer distance=3
grass74 $LOCATION/$MAPSET --exec v.to.rast --overwrite input=roadbuffer output=roads use=cat
###
### procedure VII - extract
### ... customized to use riparian information
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/g2w_subgridLULC.R $PROJDIR 101 clim/xx.base $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/vegCollection.csv $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/soilCollection.csv $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/lulcCollection.csv $PROJDIR/$RHESSysNAME/worldfiles/worldfile.csv $PROJDIR/$RHESSysNAME/worldfiles/worldfile.hdr $PROJDIR/$RHESSysNAME/defs
### ... flow table
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/createSubsurfaceFlowRouting.R $PROJDIR/$RHESSysNAME/flows/subsurfaceflow.txt
Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/LIB_RHESSys_writeTable2World.R na $PROJDIR/$RHESSysNAME/worldfiles/worldfile.csv $PROJDIR/$RHESSysNAME/worldfiles/worldfile






