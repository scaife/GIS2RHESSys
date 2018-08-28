arg=commandArgs(T)

library(rgrass7)

rast = readRAST(arg[1])
mask = !is.na(rast@data[[1]])

lulcFrac = read.csv(arg[2])
	
	# may need to customize rules below 
	# 1 = water
	# 2 = wetland
	# 3 = tree canopy
	# 4 = shrub
	# 5 = lawn
	# 6 = barren
	# 7 = impervious structure
	# 8 = impervious surface
	# 9 = impervious road
	# 10 = tree canopy over structure
	# 11 = tree canopy over impervious surface
	# 12 = tree canopy over roads
	
	deciduousCode = paste('lulc',c(2,3,10,11,12),sep='')
	shurbCode = paste('lulc',c(4),sep='')
	lawCode = paste('lulc',c(5),sep='')
	impCode = paste('lulc',c(7,8,9,10,11,12),sep='')
	roofCode = paste('lulc',c(7,10),sep='')
	parkCode = paste('lulc',c(8,11),sep='')
	roadCode = paste('lulc',c(9,12),sep='')
	
	gisOrder = match(rast@data[[1]][mask], lulcFrac$patchID)
	deciduousCode = deciduousCode[deciduousCode%in%colnames(lulcFrac)]
	shurbCode = shurbCode[shurbCode%in%colnames(lulcFrac)]
	lawCode = lawCode[lawCode%in%colnames(lulcFrac)]
	impCode = impCode[impCode%in%colnames(lulcFrac)]
	
	
	rast$forestFrac = rep(0,length(rast@data[[1]]))
	if(length(deciduousCode)>1) rast$forestFrac[mask] = (rowSums(lulcFrac[, deciduousCode])/lulcFrac$total)[gisOrder]
	if(length(deciduousCode)==1) rast$forestFrac[mask] = (lulcFrac[, deciduousCode]/lulcFrac$total)[gisOrder]
	writeRAST(rast,'forestFrac',zcol='forestFrac',overwrite=T)
	
	rast$shrubFrac = rep(0,length(rast@data[[1]]))
	if(length(shurbCode)>1) rast$shrubFrac[mask] = (rowSums(lulcFrac[, shurbCode])/lulcFrac$total)[gisOrder]
	if(length(shurbCode)==1) rast$shrubFrac[mask] = (lulcFrac[, shurbCode]/lulcFrac$total)[gisOrder]
	writeRAST(rast,'shrubFrac',zcol='shrubFrac',overwrite=T)

	rast$lawnFrac = rep(0,length(rast@data[[1]]))
	if(length(lawCode)>1) rast$lawnFrac[mask] = (rowSums(lulcFrac[, lawCode])/lulcFrac$total)[gisOrder]
	if(length(lawCode)==1) rast$lawnFrac[mask] = (lulcFrac[, lawCode]/lulcFrac$total)[gisOrder]
	writeRAST(rast,'lawnFrac',zcol='lawnFrac',overwrite=T)

	rast$impFrac = rep(0,length(rast@data[[1]]))
	if(length(impCode)>1) rast$impFrac[mask] = (rowSums(lulcFrac[, impCode])/lulcFrac$total)[gisOrder]
	if(length(impCode)==1) rast$impFrac[mask] = (lulcFrac[, impCode]/lulcFrac$total)[gisOrder]
	writeRAST(rast,'impFrac',zcol='impFrac',overwrite=T)



	rast$roofFrac = rep(0,length(rast@data[[1]]))
	if(length(roofCode)>1) rast$roofFrac[mask] = (rowSums(lulcFrac[, roofCode])/lulcFrac$total)[gisOrder]
	if(length(roofCode)==1) rast$roofFrac[mask] = (lulcFrac[, roofCode]/lulcFrac$total)[gisOrder]
	writeRAST(rast,'roofFrac',zcol='roofFrac',overwrite=T)

	rast$parkFrac = rep(0,length(rast@data[[1]]))
	if(length(parkCode)>1) rast$parkFrac[mask] = (rowSums(lulcFrac[, parkCode])/lulcFrac$total)[gisOrder]
	if(length(parkCode)==1) rast$parkFrac[mask] = (lulcFrac[, parkCode]/lulcFrac$total)[gisOrder]
	writeRAST(rast,'parkFrac',zcol='parkFrac',overwrite=T)

	rast$roadFrac = rep(0,length(rast@data[[1]]))
	if(length(roadCode)>1) rast$roadFrac[mask] = (rowSums(lulcFrac[, roadCode])/lulcFrac$total)[gisOrder]
	if(length(roadCode)==1) rast$roadFrac[mask] = (lulcFrac[, roadCode]/lulcFrac$total)[gisOrder]
	writeRAST(rast,'roadFrac',zcol='roadFrac',overwrite=T)
	
	









