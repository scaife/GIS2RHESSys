#!/bin/bash

export PATH=$PATH:/usr/local/bin
export PATH=$PATH:/Library/Frameworks/R.framework/Resources
export PATH=$PATH:/Applications/GRASS-7.4.0.app/Contents/Resources/bin:/Applications/GRASS-7.4.0.app/Contents/Resources/scripts:/etc:/usr/lib
export PATH=$PATH:~/Library/GRASS/7.4/Modules/bin:~/Library/GRASS/7.4/Modules/scripts
export GDAL_DATA=/Applications/GRASS-7.4.0.app/Contents/Resources/share/gdal

### procedure I - setup GRASS dataset
PROJDIR='Users/scaife/Dropbox/Projects/CrossSite/drought/ws14'
EPSGCODE='EPSG:26917'
RESOLUTION=10
RHESSysNAME='rhessys'
GISDBASE=$PROJDIR/grass_dataset
LOCATION_NAME=$RHESSysNAME
LOCATION=$GISDBASE/$LOCATION_NAME
MAPSET=PERMANENT
grass74 -c $EPSGCODE -e $LOCATION
###
### procedure II - setup catchment elevation
### ... rescale source elevation data to the defined $RESOLUTION (optional;no need if source data is at the same $RESOLUTION)
downloadedDEMfile='DEM.tif' # filename
grass74 $LOCATION/$MAPSET --exec r.in.gdal -e --overwrite input=$PROJDIR/raw_data/$downloadedDEMfile output=demRAW location=elevationRAW
LOCATIONDEM=$GISDBASE/elevationRAW
grass74 $LOCATIONDEM/$MAPSET --exec g.region raster=demRAW
grass74 $LOCATIONDEM/$MAPSET --exec g.region res=$RESOLUTION -a -p
grass74 $LOCATIONDEM/$MAPSET --exec r.resamp.stats -w input=demRAW output=dem$RESOLUTION'm'
grass74 $LOCATIONDEM/$MAPSET --exec r.out.gdal --overwrite input=dem$RESOLUTION'm' output=$PROJDIR/raw_data/dem$RESOLUTION'm.tif' format=GTiff
### ... import the (rescaled) elevation data into "$LOCATION/$MAPSET"; dem -> dem_ for the next step
grass74 $LOCATION/$MAPSET --exec r.in.gdal -o -e --overwrite input=$PROJDIR/raw_data/dem$RESOLUTION'm.tif' output=dem
grass74 $LOCATION/$MAPSET --exec g.region raster=dem
### ... edit the DEM for bridge removal (optional)
# grass74 $LOCATION/$MAPSET --exec r.mapcalc --overwrite expression="xmap = x()"
# grass74 $LOCATION/$MAPSET --exec r.mapcalc --overwrite expression="ymap = y()"
# grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/LinearElevationEdit.R xmap ymap dem_
###
### procedure III - setup catchment soil
## ... reproject soil data
downloadedSSURGOdirectory='NC113'
grass74 $LOCATION/$MAPSET --exec v.in.ogr --overwrite input=$PROJDIR/raw_data/$downloadedSSURGOdirectory/spatial/soilmu_a_"$(echo $downloadedSSURGOdirectory | tr '[A-Z]' '[a-z]')".shp output=ssurgo location=soilRAW
LOCATIONSOIL=$GISDBASE/soilRAW
### ... import (reprojected) soil data to "$LOCATION/$MAPSET" and convert shape-polygon to raster soil (USDA) classes
grass74 $LOCATION/$MAPSET --exec v.proj --overwrite location=soilRAW mapset=PERMANENT input=ssurgo output=ssurgo
grass74 $LOCATION/$MAPSET --exec v.to.rast --overwrite input=ssurgo use=cat output=soil_ssurgo

## need to write a RHESSysNAME directory
mkdir $PROJDIR/$RHESSysNAME
grass74 $LOCATION/$MAPSET --exec v.db.select map=ssurgo separator=comma file=$PROJDIR/$RHESSysNAME/soil_cat_mukey.csv
git clone git://github.com/scaife/ssurgo_extraction $PROJDIR/$RHESSysNAME/ssurgo_extraction-master
Rscript $PROJDIR/$RHESSysNAME/ssurgo_extraction-master/ssurgo_extraction.R $PROJDIR/raw_data/$downloadedSSURGOdirectory
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/ssurgo_extraction-master/ssurgo_soiltexture2gis.R $PROJDIR/$RHESSysNAME/soil_cat_mukey.csv $PROJDIR/raw_data/$downloadedSSURGOdirectory/soil_mukey_texture.csv
### ... use soil to define riparian are in this case (optional)
# -- CIS --- Need to define riparian. You can bypass by setting to stream map. 
# grass74 $LOCATION/$MAPSET --exec r.mapcalc --overwrite expression="riparian = if(soil_texture==8,1,null())"
###
### procedure IV - delineate watershad & basic rhessys spatial structures
### ... if input is Lat/Long
#gageLat='39.47947' # catchment outlet WSG84 Lat (decimal degree)
#gageLong='-76.67803' # catchment outlet WSG84 Long (decimal degree; includes the negative sign if applied)
#declare $(grass74 $LOCATION/$MAPSET --exec m.proj -i coordinates=$gageLong,$gageLat separator=space | awk '{print "xyCoord=" $1 "," $2}')
#echo $xyCoord | grass74 $LOCATION/$MAPSET --exec v.in.ascii in=- out=test x=1 y=2 separator=, --overwrite
### ... if input is a shapefile point
gageShapefile='coweeta_weirs_shp.shp'
grass74 $LOCATION/$MAPSET --exec v.in.ogr --overwrite input=$PROJDIR/raw_data/$gageShapefile output=outlet location=outletRAW
LOCATIONOUTLET=$GISDBASE/outletRAW
grass74 $LOCATIONOUTLET/$MAPSET --exec v.extract --overwrite input=outlet type=point where="COMMENT = 'Weir 14'" output=gage
grass74 $LOCATION/$MAPSET --exec v.proj --overwrite location=outletRAW mapset=PERMANENT input=gage output=outlet
# declare $(v.out.ascii input=outlet type=point format=point separator=space | awk '{print "xyCoord=" $1 "," $2}')
thres='10000' # area in squared meters
expectedDrainageArea=640000 # meter square (allow 2% error)
git clone git://github.com/scaife/GIS2RHESSys $PROJDIR/$RHESSysNAME/GIS2RHESSys-master
grass74 $LOCATION/$MAPSET --exec sh $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/grass_delineation.sh $(bc <<< "$thres/$RESOLUTION/$RESOLUTION") $(bc <<< "0.98*$expectedDrainageArea/$RESOLUTION/$RESOLUTION") $(bc <<< "1.02*$expectedDrainageArea/$RESOLUTION/$RESOLUTION")
### ... set zone (optional)
# grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/zone_cluster.R dem slope aspect hill # define this later
###

# procedure - define landcover
downloadedLULCfile='NLCD.tif' # full path to the downloaded <file>
grass74 $LOCATION/$MAPSET --exec r.in.gdal -e --overwrite input=$PROJDIR/raw_data/$downloadedLULCfile output=lulcRAW location=lulcRAW
LOCATIONLULC=$GISDBASE/lulcRAW
grass74 $LOCATIONLULC/$MAPSET --exec g.region raster=lulcRAW
grass74 $LOCATIONLULC/$MAPSET --exec g.region res=$RESOLUTION -a -p
grass74 $LOCATIONLULC/$MAPSET --exec r.resamp.interp --overwrite input=lulcRAW output=lulc$RESOLUTION'm' method=nearest
grass74 $LOCATIONLULC/$MAPSET --exec r.out.gdal --overwrite input=lulc$RESOLUTION'm' output=$PROJDIR/raw_data/lulc$RESOLUTION'm.tif' format=GTiff
grass74 $LOCATION/$MAPSET --exec r.in.gdal -o -e --overwrite input=$PROJDIR/raw_data/lulc$RESOLUTION'm.tif' output=LULCcode
grass74 $LOCATION/$MAPSET --exec g.region raster=LULCcode
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/NLCD2RHESSys.R

### procedure V - LULC, vectorize "patch" from "$LOCATION/$MAPSET" and calculate
# grass74 $LOCATION/$MAPSET --exec r.to.vect input=patch output=patch type=area
# grass74 $LOCATION/$MAPSET --exec r.in.gdal -e --overwrite input=$PROJDIR/raw_data/$downloadedLULCfile output=lulcRAW location=lulcRAW
# LOCATIONLULC=$GISDBASE/lulcRAW
# grass74 $LOCATIONLULC/$MAPSET --exec r.colors map=lulcRAW color=random
# grass74 $LOCATIONLULC/$MAPSET --exec v.proj location=$LOCATION_NAME mapset=PERMANENT input=patch output=patch$RESOLUTION'm'
# grass74 $LOCATIONLULC/$MAPSET --exec v.to.rast input=patch$RESOLUTION'm' output=patch$RESOLUTION'm' use=attr attribute_column=value
# grass74 $LOCATIONLULC/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/aggregate_lulcFrac.R patch$RESOLUTION'm' lulcRAW $PROJDIR/$RHESSysNAME/lulcFrac$RESOLUTION'm.csv'
# grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/aggregate_lulcFrac_write2GIS.R patch $PROJDIR/$RHESSysNAME/lulcFrac$RESOLUTION'm.csv'
###
### procedure VI - import road network shape file
# downloadedROADfile='roads.shp'
# grass74 $LOCATION/$MAPSET --exec v.in.ogr --overwrite input=$PROJDIR/raw_data/$downloadedROADfile output=roads location=roadRAW
# LOCATIONROAD=$GISDBASE/roadRAW
# grass74 $LOCATION/$MAPSET --exec v.proj --overwrite location=roadRAW mapset=PERMANENT input=roads output=roads
# grass74 $LOCATION/$MAPSET --exec v.to.rast --overwrite input=roads@PERMANENT output=roads use=cat

### procedure VI - import road network TIF
downloadedROADfile='roads.tif'
grass74 $LOCATION/$MAPSET --exec r.in.gdal --overwrite input=$PROJDIR/raw_data/$downloadedROADfile output=roads location=roadRAW
LOCATIONROAD=$GISDBASE/roadRAW
grass74 $LOCATIONROAD/$MAPSET --exec g.region raster=roads
grass74 $LOCATIONROAD/$MAPSET --exec g.region res=$RESOLUTION -a -p
grass74 $LOCATIONROAD/$MAPSET --exec r.resamp.stats -w input=roads output=roads$RESOLUTION'm'
grass74 $LOCATIONROAD/$MAPSET --exec r.out.gdal --overwrite input=roads$RESOLUTION'm' output=$PROJDIR/raw_data/roads$RESOLUTION'm.tif' format=GTiff
grass74 $LOCATION/$MAPSET --exec r.in.gdal -o -e --overwrite input=$PROJDIR/raw_data/roads$RESOLUTION'm.tif' output=roads
grass74 $LOCATION/$MAPSET --exec g.region raster=roads

downloadedISOHYETfile='isohyet.tif'
grass74 $LOCATION/$MAPSET --exec r.in.gdal --overwrite input=$PROJDIR/raw_data/$downloadedISOHYETfile output=isohyet location=isohyetRAW
LOCATIONISOHYET=$GISDBASE/isohyetRAW
grass74 $LOCATIONISOHYET/$MAPSET --exec g.region raster=isohyet
grass74 $LOCATIONISOHYET/$MAPSET --exec g.region res=$RESOLUTION -a -p
grass74 $LOCATIONISOHYET/$MAPSET --exec r.resamp.stats -w input=isohyet output=isohyet$RESOLUTION'm'
grass74 $LOCATIONISOHYET/$MAPSET --exec r.out.gdal --overwrite input=isohyet$RESOLUTION'm' output=$PROJDIR/raw_data/isohyet$RESOLUTION'm.tif' format=GTiff
grass74 $LOCATION/$MAPSET --exec r.in.gdal -o -e --overwrite input=$PROJDIR/raw_data/isohyet$RESOLUTION'm.tif' output=isohyet
grass74 $LOCATION/$MAPSET --exec g.region raster=isohyet

### procedure VII - define riparian
grass74 $LOCATION/$MAPSET --exec r.mapcalc --overwrite expression="riparian = if(str>0,1,null())"

### procedure XX climate data
climateBaseFile='clim/cwt.base'
climateBaseID=101

### procedure VII - extract
### ... customized to use riparian information
mkdir $PROJDIR/$RHESSysNAME/defs
mkdir $PROJDIR/$RHESSysNAME/worldfiles
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/g2w.R $PROJDIR $climateBaseID $climateBaseFile $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/vegCollection.csv $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/soilCollection.csv $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/lulcCollection.csv $PROJDIR/$RHESSysNAME/worldfiles/worldfile.csv $PROJDIR/$RHESSysNAME/worldfiles/worldfile.hdr $PROJDIR/$RHESSysNAME/defs

Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/LIB_RHESSys_writeTable2World.R na $PROJDIR/$RHESSysNAME/worldfiles/worldfile.csv $PROJDIR/$RHESSysNAME/worldfiles/worldfile

mkdir $PROJDIR/$RHESSysNAME/flows
grass74 $LOCATION/$MAPSET --exec Rscript $PROJDIR/$RHESSysNAME/GIS2RHESSys-master/createFlowRouting.R $PROJDIR/$RHESSysNAME/flows/flowtable.txt

# running rhessys
start_date='2000 1 1 1'
end_date='2001 10 1 1'
RHESSysOUTPUTDIR='output'
mkdir $PROJDIR/$RHESSysNAME/$RHESSysOUTPUTDIR

# -b only basin output; -gwtoriparian receiving groundwater and put in stream
RHESSys_flags = "-b -newcaprise -capr 0.001 -capMax 0.01 -slowDrain"
RHESSys_worldfile = "-w worldfiles/worldfile -whdr worldfiles/worldfile.hdr"
RHESSys_flowtable = "-r flows/flowtable.txt"
RHESSys_tecfile = "-t tecfiles/tec_daily.txt"
RHESSys_output = '-pre '+RHESSysOUTPUTDIR+'/rhessys'

# write tec file
mkdir $PROJDIR/$RHESSysNAME/tecfiles
tecfile=$PROJDIR/$RHESSysNAME/tecfiles/tec_daily.txt
echo '2001 1 1 1 print_daily_on' >$tecfile

# download rhessys
git clone --single-branch --branch Feb27_2019 git://github.com/laurencelin/RHESSys5.20.EastCoastUS $PROJDIR/$RHESSysNAME/src
unzip $PROJDIR/$RHESSysNAME/src/rhessys5.20.zip -d $PROJDIR/$RHESSysNAME/src

# cd to makefile and make

cp $PROJDIR/$RHESSysNAME/src/rhessys5.20/rhessys/rhessys5.20.0.develop $PROJDIR/$RHESSysNAME

# change directory into modle folder and run from there. 
./rhessys5.20.0.develop -st 2000 1 1 1 -ed 2001 1 1 1 -b -newcaprise -capr 0.001 -capMax 0.01 -slowDrain -w worldfiles/worldfile -whdr worldfiles/worldfile.hdr -r flows/flowtable.txt -t tecfiles/tec_daily.txt -pre output/test
