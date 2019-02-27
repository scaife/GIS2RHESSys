## flow table
options(scipen=999)
arg=commandArgs(T)
library(rgrass7)
library(rgdal)
gis = gmeta()

DtoR = pi/180
RtoD = 1/DtoR
BASEMENT_DEPTH = 3 #meter
PAVEDROAD_DEPTH = 0.3 #meter
# bounded by GIS mask

	## ... basic maps
	basinMap = 'basin'
	hillslopeMap = 'hill'
	zoneMAP = 'zone'
	xMap = 'xmap'
	yMap = 'ymap'
	patchMAP = 'patch'
	rowMap = 'rowmap' ##<<--- raster calculator row()
	colMap = 'colmap' ##<<--- raster calculator col()
	demMap = 'dem'
	slopeMap = 'slope' ## unit in degree
	
	rast0 = readRAST(c(basinMap, hillslopeMap, zoneMAP, patchMAP, rowMap, colMap),NODATA=0)
		mask = !is.na(rast0@data[[1]])
		basin = rast0@data[[1]][mask]
		hill = rast0@data[[2]][mask]
		zone = rast0@data[[3]][mask]
		patch = rast0@data[[4]][mask]
		rows = rast0@data[[5]][mask]
		cols = rast0@data[[6]][mask]
	rast1 = readRAST(c(demMap, xMap, yMap, slopeMap))
		dem = rast1@data[[1]][mask]
		xx = rast1@data[[2]][mask]
		yy = rast1@data[[3]][mask]
		slope = rast1@data[[4]][mask]
		
	## ... stream, drainage maps
	roadMap = 'roads'
	streamMap = 'strShort' # 'str62ha2', 'reducedChannel'
	streamFullextension = 'str' #
	streamSurfaceDrain = 'str' #
	drainMap = 'drain'
	road_impFracMAP = 'roadFrac' ##<<------------ for road width estimate
	
	rast2 = readRAST(c(streamMap, streamFullextension, streamSurfaceDrain))	
		stream = rast2@data[[1]][mask]
		fullstreamExt = rast2@data[[2]][mask]
		additionalStreamDrain = rep(NA,sum(mask)) #rast2@data[[2]][mask]
		
	rast3 = readRAST(roadMap)		
		road = rast3@data[[1]][mask]
	
	rast4 = readRAST(drainMap)		
		drain = abs(rast4@data[[1]][mask])
		
	rast5 = readRAST(road_impFracMAP)
		road_impFrac = rast5@data[[1]][mask]
			
	## ... ACTION maps
	riparianMAP = 'riparian' #<---- GW to riparian action
	sewercoverMAP = 'sewercover' #<---- water table drain [0-3 m]
	stormdrainMAP = 'roadExit' # "logically sinigicant" storm drain entries along the road network / certain locations
	roofMAP = 'roofFrac'
	parkingMAP = 'parkFrac'
	roadDrainMAP = 'roadDrain'
	lawnMAP = 'lawnFrac'
	naturalMAP = 'naturalLand'
	
	roofNextStreamFile = 'roofFracPtLink.txt'
	parkingNextStreamFile = 'parkFracPtLink.txt'
	stormdrainNextStreamFile = 'roadExitPtLink.txt'
	
	rast6 = readRAST(c(riparianMAP, naturalMAP))
		riparian = rast6@data[[1]][mask]
		natural = rast6@data[[2]][mask] #<<-------------------
		
	rast7 = readRAST(sewercoverMAP)
		sewercover = rast7@data[[1]][mask]
		
	rast8 = readRAST(stormdrainMAP)
		stormdrain = rast8@data[[1]][mask]		
		
	rast9 = readRAST(c(roofMAP, parkingMAP,lawnMAP))
		roofFrac = rast9@data[[1]][mask]	
		parkingFrac = rast9@data[[2]][mask]	
		lawnFrac = rast9@data[[3]][mask]	
		
	rast10 = readRAST(roadDrainMAP)
		roadDrain = rast10@data[[1]][mask]	

		
	## assume grids are squares 
	cellarea = gis$nsres * gis$ewres
	cellsize = sqrt(cellarea) 
	flatDEMdrop = tan(DtoR*0.143)*cellsize # only 0.25m drop per 100m.
	roadWidth = 5
	roadWidth = ifelse(cellsize>=9.99,
		5, # meters (default)
		ifelse(roadWidth>cellsize, cellsize, roadWidth)
	 )# ifelse
	directEdge = cellsize*0.5
	diagonalEdge = cellsize*sqrt(0.5)	
		
	# 1.  2. 3.  4. 5.  6. 7.  8. (GRASS from current drainTO code order)
	# NE, N, NW, W, SW, S, SE, E
	colneighbor = c(1,	0,	-1,	-1,	-1,	0,	1,	1)	
	rowneighbor = c(-1,	-1,	-1,	0,	1,	1,	1,	0)	
	directEdgeIndex = c(2,4,6,8)
	indirectEdgeIndex = c(1,3,5,7)
	 
	maxCol = max(cols,na.rm=T) 
	maskRC = rows*maxCol+cols #paste(rows, cols,sep=':') ## row*[max col]+col (yes: unique ID)
	maskRC_string2Patch_num <- new.env(hash=T)
	list2env(setNames(as.list(patch),maskRC),envir=maskRC_string2Patch_num) #<<---- native R hash
	gridSurroundRC = sapply(rows, FUN=function(x){x+rowneighbor})*maxCol+sapply(cols, FUN=function(x){x+colneighbor})
	gridSurroundRC[!(gridSurroundRC %in% maskRC)] = -1
	
	# part 1: gathering information to temporary files
	#patch_title = c('patchID','dem','xx','yy','basin','hill','zone','rr','cc','grid','strQ','roadQ','accgrid','mslope','Mslope')

	fullLength = seq(1,length.out=length(patch))
	patch_info_dem = tapply(dem,INDEX=patch,mean)
	orderedPatch = as.numeric(names(patch_info_dem[order(patch_info_dem,decreasing=T)])) ### patch could be longer than 'orderedPatch'
	outputOrder = match(patch, orderedPatch) # has the same length as 'patch'
		# test = tapply(patch, outputOrder,mean); sum(test==orderedPatch)==length(orderedPatch)
	maskRC_string2outputOrder_num <- new.env(hash=T)
	list2env(setNames(as.list(outputOrder),maskRC),envir=maskRC_string2outputOrder_num) #<<---- native R hash
	maskRC_string2maskRC_num <- new.env(hash=T)
	list2env(setNames(as.list(maskRC),maskRC),envir= maskRC_string2maskRC_num) #<<---- native R hash
	
	print('starting step I')
	patchInfo = tapply(fullLength, INDEX=outputOrder, FUN=function(ii){
		
		return <- c(
			mean(patch[ii]), 			#1 patchID
			mean(dem[ii]),				#2 elevation
			mean(xx[ii]),				#3 x coordinate
			mean(yy[ii]),				#4 y coordinate
			mean(basin[ii]),			#5 basinID
			mean(hill[ii]),				#6 hillID
			mean(zone[ii]),				#7 zoneID
			mean(rows[ii]),				#8 row index (from left to right)
			mean(cols[ii]),				#9 col index (from top to bottom)
			length(ii),					#10 num of cells
			sum(!is.na(stream[ii])),	#11 strQ (checking whether patch contains stream grids)
			sum(!is.na(road_impFrac[ii]) & road_impFrac[ii]>0),	#12  (checking whether patch contains road grids); roadQ sum(!is.na(road[ii]))
			tan(mean(slope[ii])*DtoR),	#13 average slope
			tan(max(slope[ii])*DtoR),	#14 max slope
			mean(road_impFrac[ii],na.rm=T), 	#15 for road width calculation & average road Frac
			sum(!is.na(riparian[ii])),	#16 (checking whether patch contains riparian grids)
			sum(!is.na(sewercover[ii])), #17 (checking whether patch contains sewercover grids)
			sum(!is.na(stormdrain[ii])), #18 (checking whether patch contains stormdrain entry grids)
			mean(roofFrac[ii],na.rm=T), #19 average roof Frac 
			mean(parkingFrac[ii],na.rm=T), #20 average parking Frac 
			sum(!is.na(fullstreamExt[ii])),	#21 strQ (checking whether patch contains fullstreamExt grids); use to correct road intercepts
			sum(!is.na(additionalStreamDrain[ii])), # 22 non-stream cell but need surface water drain
			mean(lawnFrac[ii],na.rm=T), # 23 lawnFrac for irrigate
			mean(natural[ii],na.rm=T) # 24 natural
			);
	})#tapply <--- this output is a list of c() in outputOrder
	patch_info_lowest = patchInfo[[ length(patchInfo) ]]
	 
	 
	 
	 
	 
	## part 2: sort by 'elevation' & finding neighbor 
	print('starting step II') 
		
	## .......... Neighbour
		# ii=which(orderedPatch==24110) #500
		patchNeighbourRC_edge = tapply(fullLength, INDEX=outputOrder, FUN=function(jj){
			withinPatchGridRC = rows[jj]*maxCol+cols[jj]; # within
	
			hold = as.vector(gridSurroundRC[directEdgeIndex,jj]);
			hold[hold%in% withinPatchGridRC] = -1
			
			hold2 = as.vector(gridSurroundRC[indirectEdgeIndex,jj]);
			hold2[hold2%in% withinPatchGridRC] = -1
	
			#return <- c(table(hold[hold>0])* directEdge, table(hold2[hold2>0])* diagonalEdge) 
			return <- c( tapply(hold[hold>0],hold[hold>0],length)*directEdge, tapply(hold2[hold2>0],hold2[hold2>0],length)*diagonalEdge) 
			 
		})
			
		patchNeighbourRC_LEN = seq_along(patchNeighbourRC_edge) ## <<--------------------------- ordered patch aggregated neighbours
		patchNeighbourRC = sapply(patchNeighbourRC_LEN, function(ii){
			sapply(names(patchNeighbourRC_edge[[ii]]), function(x){maskRC_string2maskRC_num[[ x ]]})
		})
		patchNeighbourPatch = sapply(patchNeighbourRC_LEN, function(ii){
			sapply(names(patchNeighbourRC_edge[[ii]]), function(x){maskRC_string2Patch_num[[ x ]]})
		})
		patchNeighbourPatchIndex = sapply(patchNeighbourRC_LEN, function(ii){
			sapply(names(patchNeighbourRC_edge[[ii]]), function(x){maskRC_string2outputOrder_num[[ x ]]})
		})


	## .......... prefer Neighbour
		patchPreferNeighbourRC = tapply(fullLength, INDEX=outputOrder, FUN=function(ii){
			withinPatchGridRC = rows[ii]*maxCol+cols[ii]; # within
			drainTO_index = cbind(drain[ii],ii)
			hold3 = as.vector(gridSurroundRC[ drainTO_index ])
			hold3[hold3%in% withinPatchGridRC] = -1
			
			return <- sapply(names(tapply(hold3[hold3>0], hold3[hold3>0],length)),function(x){maskRC_string2maskRC_num[[ x ]]})
		})
	
	

	## .......... surface (road drain along the road)
		#roadDrain = rast10@data[[1]][mask]
		patchRoadDrainRC = tapply(fullLength, INDEX=outputOrder, FUN=function(ii){
			withinPatchGridRC = rows[ii]*maxCol+cols[ii]; # within
			if(is.na(roadDrain[ii])){
				return <- -1
			}else if(roadDrain[ii]<0){
				return <- -1
			}else{
				drainTO_index = cbind(roadDrain[ii],ii)
				hold3 = as.vector(gridSurroundRC[ drainTO_index ])
				hold3[hold3%in% withinPatchGridRC] = -1
				
				return <- sapply(names(tapply(hold3[hold3>0], hold3[hold3>0],length)),function(x){maskRC_string2maskRC_num[[ x ]]})
			}
		})# tapply
		
		patchRoadDrainPatchIndex = sapply(fullLength, function(ii){
			if( patchRoadDrainRC[ii]>0 ){
				return <- sapply(toString(patchRoadDrainRC[ii]), function(x){maskRC_string2outputOrder_num[[ x ]]})
			}else{
				return <- -1
			}
		})# 3054
	
	
	## .......... storm drain (with roads / some locations)
		stormdrainNextStream = read.table(paste(stormdrainNextStreamFile,sep=''),header=T)	
		stormdrainNextStream_patchIndex = sapply(seq_len(dim(stormdrainNextStream)[1]),function(ii){
			return <- which.min( abs(xx-stormdrainNextStream$to_x[ii]) + abs(yy-stormdrainNextStream$to_y[ii]) )
		})# sapply
		
		stormdrainNextStream_patchIndex_LEN = 1:length(stormdrainNextStream_patchIndex)
		patchStormdrainPatchIndex = tapply(fullLength, INDEX=outputOrder, FUN=function(ii){
			# ii = 3074; outputOrder[35]
			# ii = 35
			patchID = patch[ii][1];
			cond = stormdrainNextStream$from_cat == patchID
			if(sum(cond)>0){
				return <- unique(stormdrainNextStream_patchIndex[stormdrainNextStream_patchIndex_LEN[cond]])
			}else{
				return <- -1
			}
		})# tapply
	
	
	## .......... surface (roof drains to road)
		roofNextStream = read.table(paste(roofNextStreamFile,sep=''),header=T)	
		roofNextStream_patchIndex = sapply(seq_len(dim(roofNextStream)[1]),function(ii){
			return <- which.min( abs(xx-roofNextStream$to_x[ii]) + abs(yy-roofNextStream$to_y[ii]) )
		})# sapply
		
		roofNextStream_patchIndex_LEN = 1:length(roofNextStream_patchIndex)
		patchRoofPatchIndex = tapply(fullLength, INDEX=outputOrder, FUN=function(ii){
			# ii = 3074; outputOrder[35]
			# ii = 35
			patchID = patch[ii][1];
			cond = roofNextStream$from_cat == patchID
			if(sum(cond)>0){
				return <- unique(roofNextStream_patchIndex[roofNextStream_patchIndex_LEN[cond]])
			}else{
				return <- -1
			}
		})# tapply
	
	
	## .......... surface (parking drains to road)
		parkingNextStream = read.table(paste(parkingNextStreamFile,sep=''),header=T)	
		parkingNextStream_patchIndex = sapply(seq_len(dim(parkingNextStream)[1]),function(ii){
			return <- which.min( abs(xx-parkingNextStream$to_x[ii]) + abs(yy-parkingNextStream$to_y[ii]) )
		})# sapply
		
		parkingNextStream_patchIndex_LEN = 1:length(parkingNextStream_patchIndex)
		patchParkingPatchIndex = tapply(fullLength, INDEX=outputOrder, FUN=function(ii){
			# ii = 3074; outputOrder[35]
			# ii = 35
			patchID = patch[ii][1];
			cond = parkingNextStream$from_cat == patchID
			if(sum(cond)>0){
				return <- unique(parkingNextStream_patchIndex[parkingNextStream_patchIndex_LEN[cond]])
			}else{
				return <- -1
			}
		})# tapply	
	

	
	## ...........................................
	## .......... writing out flow table
	arg=c('flows/subTest9d.txt','flows/surfTest9d.txt')
	 
	subsurfaceflow_table_buff <- file(arg[1],'w') # open a file connection
	surfaceflow_table_buff <- file(arg[2],'w') # open a file connection
	cat( length(patchInfo), '\n', file=subsurfaceflow_table_buff) #,sep='\n'
	cat( length(patchInfo), '\n', file=surfaceflow_table_buff) #,sep='\n'
	
	#silent = sapply(patchNeighbourRC_LEN, function(ii){
	for(ii in patchNeighbourRC_LEN){

		withinNeighbourRC_edge = patchNeighbourRC_edge[[ii]] 	
		withinNeighbourRC = patchNeighbourRC[[ii]]					
		withinNeighbourRC_prefer = rep(0,length(withinNeighbourRC))			##<<------
			withinNeighbourRC_prefer[withinNeighbourRC%in%patchPreferNeighbourRC[[ii]] ] = 1
		index4neighbour = patchNeighbourPatchIndex[[ii]] 
		 			 	
		current_patch_info = patchInfo[[ii]]
			# 0 = land (default)
			# 1 = class::stream
			# 2 = class::road
			# 3 = actionSTORMDRAIN
			# 5 = actionGWDRAIN
			# 7 = actionRIPARIAN
			# 11 = actionSEWER
			# 13 = actionIRRIGRATION 
			# 17 = actionPIPEDRAIN
			
			# current_patch_info[12] = roadFrac
			# current_patch_info[16] = riparian -> 7 riparian
			# current_patch_info[17] = sewercover -> 11 sewer drain
			# current_patch_info[18] = stormdrain -> stormdrain
			# current_patch_info[12] = road grid -> as LAND but no GW drain
			# current_patch_info[15] = road Frac -> as LAND but no GW drain
			# current_patch_info[19] = roof Frac -> as LAND but no GW drain
			# current_patch_info[20] = parking Frac -> as LAND but no GW drain
			# current_patch_info[21] = stream extension
			# current_patch_info[23] = lawnFrac
			
		## actionCode is mostly for subsurface processes; surface storm drain see below.
		actionCode = 	ifelse(current_patch_info[15]>0|current_patch_info[19]>0|current_patch_info[20]>0|current_patch_info[12]>0,1,5) * # GW drain
						ifelse(current_patch_info[17]>0,11,1) * # sewer drain (top 3-m)
						ifelse(current_patch_info[16]>0,7,1) * # riparian
						ifelse(current_patch_info[18]>0,ifelse(current_patch_info[21]>0,1,3), 1) * #storm drain <<--- map "roadExit" -- no in use
						ifelse(current_patch_info[23]>0,13,1) * # lawn irrigration
						ifelse(current_patch_info[12]>0,17,1) # subsurface storm drain / pipe drain (top 1-m) under the road
							
		drainage_type = ifelse(current_patch_info[11]>0, 1, # class::stream
						ifelse(actionCode>1, actionCode,0)
						)			
		
				
		
		neighbourLength = 1:length(withinNeighbourRC)	
		neighbourOrder = match(withinNeighbourRC,unique(withinNeighbourRC))
		
		allNeighbourInfo = simplify2array(tapply(neighbourLength, INDEX=neighbourOrder, function(jj){
			## exploring information between "current" and neighbour(jj)
				# index4neighbour[jj][1] # index of neighbour(jj) in "patchInfo" list
				# withinNeighbourRC_edge[jj] # all edges between current and neighbour(jj)
				# withinNeighbourRC_prefer[jj] # all prefers between current and neighbour(jj)
			
			neighbor_patch_info = patchInfo[[ index4neighbour[jj][1] ]];
			idiffDEM = current_patch_info[2]-neighbor_patch_info[2]
			idiffDEM = ifelse(idiffDEM<0,0, idiffDEM)
			
			return <- c(
				neighbor_patch_info[c(1,7,6)], #patchID, zone, hill [1,2,3]
				## ... distance
				sqrt((neighbor_patch_info[3]-current_patch_info[3])^2 + 
					(neighbor_patch_info[4]-current_patch_info[4])^2), # distance [4]
				## ... rise
				idiffDEM, #rise (local prefer) [5] # zero correct
				ifelse( mean(withinNeighbourRC_prefer[jj])>0, flatDEMdrop, 0), #rise (regional prefer) [6] # zero correct
				## ... shared edge
				sum(withinNeighbourRC_edge[jj]), #edge [7]
				## info
				neighbor_patch_info[c(19,15,20,23)]# roof 8 /road 9 /parking 10 / lawn 11
			)
		}))#tapply <<--- not in a right order
		
		## local prefer 
		slope_jj_l = allNeighbourInfo[5,]/allNeighbourInfo[4,] # rise / distance
		gamma_jj_l = slope_jj_l*allNeighbourInfo[7,] # edge (width)
		
		## regional prefer 
		slope_jj_r = allNeighbourInfo[6,]/allNeighbourInfo[4,] # rise / distance
		gamma_jj_r = slope_jj_r*allNeighbourInfo[7,]
		
		cc1 = sum(gamma_jj_l)==0 # use gamma_jj_r
		cc2 = sum(gamma_jj_l < gamma_jj_r)==0 	## T: use gamma_jj_l	
		cc3 = sum(gamma_jj_l) > sum(gamma_jj_r)	
		if(cc1){
			gamma_jj = gamma_jj_r
			selectedFlow2neigbour = slope_jj_r>0
		}else if(cc2){
			gamma_jj = gamma_jj_l
			selectedFlow2neigbour = slope_jj_l>0
		}else if(cc3){
			gamma_jj = gamma_jj_l/sum(gamma_jj_l)*0.3 + (gamma_jj_r>0)*0.7
			selectedFlow2neigbour = slope_jj_r>0 | slope_jj_l>0
		}else{
			gamma_jj = gamma_jj_l + gamma_jj_r
			selectedFlow2neigbour = slope_jj_r>0 | slope_jj_l>0
		}	
		
		## ... neighbour gamma fraction
		neighbor_frac_gamma = gamma_jj/ifelse(sum(gamma_jj)>0,sum(gamma_jj),1)
		
		## ... total_gamma
		total_perimeter = sum( allNeighbourInfo[7, selectedFlow2neigbour] )
		total_gamma = sum(gamma_jj)/total_perimeter*current_patch_info[10]*cellarea; # currrent CF calculation
		if(drainage_type==1) total_gamma = current_patch_info[13]*current_patch_info[10]*cellarea; # special for stream
			
	
		cat(
			paste(current_patch_info[c(1,7,6)], collapse=' '),
			paste(sprintf('%.1f',current_patch_info[c(8,9,2)]), collapse=' '),
			sprintf('%.2f',1.0),
			sprintf('%.2f', ifelse(!is.na(current_patch_info[19]),current_patch_info[19],0)*BASEMENT_DEPTH + ifelse(!is.na(current_patch_info[15]),current_patch_info[15],0)*PAVEDROAD_DEPTH + ifelse(!is.na(current_patch_info[20]),current_patch_info[20],0)*PAVEDROAD_DEPTH),
			drainage_type, 
			total_gamma, length(withinNeighbourRC),'\n', file=subsurfaceflow_table_buff,sep=' ')
		
		cat( paste(
			allNeighbourInfo[1,],
			allNeighbourInfo[2,],
			allNeighbourInfo[3,], 
			sprintf('%.5f',neighbor_frac_gamma),
			sprintf('%.2f',allNeighbourInfo[7,]/allNeighbourInfo[4,]),
			sprintf('%.2f',allNeighbourInfo[7,]),sep=' '), file=subsurfaceflow_table_buff,sep='\n')
		
		if(drainage_type==2) cat (
			patch_info_lowest[1], 
			patch_info_lowest[7], 
			patch_info_lowest[6], 
			cellsize*current_patch_info[15],'\n', file=subsurfaceflow_table_buff,sep=' ')  #*current_patch_info[16]
			
			
		
		#-------------- surface -------------------#
		if( (current_patch_info[19]>0 | current_patch_info[20]>0 | current_patch_info[15]>0 | current_patch_info[22]>0) & current_patch_info[11]==0 ){
			# debug: surface water is the "detention" in the model; not calculated by gamma/total_gamma
			# debug: road grids (by frac) has only one drain direction; storm drain will overwrite every direction; 
			# debug: parking/roof has partial direction by their frac of the grid
			# this scheme (roof/parking to road and then storm) does not work in RHESSys because surface water only move by one grid per day!
			# correction: we need to teleport surface water from root/parking/road to stream (storm drain targets)
			
			# surface
			# what needs to be surface: 
			# roof / parking / road 
			
			stormsurfacedrainFrac = c(
				ifelse(current_patch_info[19]>0 & current_patch_info[21]==0,current_patch_info[19],0),  # roof Frac
				ifelse(current_patch_info[20]>0 & current_patch_info[21]==0,current_patch_info[20],0),  # parking Frac
				ifelse(current_patch_info[15]>0 & current_patch_info[21]==0,current_patch_info[15],0),  # road Frac
				ifelse(current_patch_info[17]>0 & current_patch_info[21]==0,0.5,0) # use [17] sewercover for surface drain area (extended from roof/road/parking); 
			); names(stormsurfacedrainFrac) = c('roof','parking','road','extenddrain')
			if(sum(stormsurfacedrainFrac[1:3])>=1){ stormsurfacedrainFrac[1:3] = stormsurfacedrainFrac[1:3]/sum(stormsurfacedrainFrac[1:3]); stormsurfacedrainFrac[4]=0; }else if(stormsurfacedrainFrac[4]>0 & sum(stormsurfacedrainFrac[1:4])>=1){ stormsurfacedrainFrac[4] = 1 - sum(stormsurfacedrainFrac[1:3]); }
			
			normal_neighborNum = length(neighbor_frac_gamma)
			normalFrac = 1 - sum(stormsurfacedrainFrac) - ifelse(current_patch_info[22]>0,0.6,0);    # current_patch_info[21] >= current_patch_info[22] !!
			if(normalFrac<0){ print(paste(ii, normalFrac)); normalFrac = 0;}
			normal_neighbor_frac_gamma = (neighbor_frac_gamma * normalFrac)
			
			## stop routing from roof to roof on surface
			normal_neighbor_frac_gammaSUM = sum(normal_neighbor_frac_gamma)
			normal_neighbor_frac_gamma = normal_neighbor_frac_gamma * (1.0 - allNeighbourInfo[8,]) # roofFrac
			if(sum(normal_neighbor_frac_gamma)>0){normal_neighbor_frac_gamma = normal_neighbor_frac_gamma/sum(normal_neighbor_frac_gamma)*normal_neighbor_frac_gammaSUM;}else{ normal_neighbor_frac_gamma = rep(0,length(normal_neighbor_frac_gamma)) }
			
			if(current_patch_info[15]>0 & current_patch_info[21]==0) normal_neighborNum = normal_neighborNum + 1;
			if(current_patch_info[19]>0 & current_patch_info[21]==0) normal_neighborNum = normal_neighborNum + 1; ## direct apply to patch[]
			if(current_patch_info[20]>0 & current_patch_info[21]==0) normal_neighborNum = normal_neighborNum + 1; ## direct apply to patch[]	 
			if( stormsurfacedrainFrac['extenddrain']>0 ) normal_neighborNum = normal_neighborNum + 1; ## direct apply to patch[]	     
			if(current_patch_info[22]>0) normal_neighborNum = normal_neighborNum + 1;
			## 21 is the full stream extension;  11 is the modeled stream grid; 22 is non-modeled stream grid bounded by full stream extension
			## note that put all the flow to the outlet for now because surface routing in RHESSys is not correct (not similar to my stream model).
			
			cat(
				paste(current_patch_info[c(1,7,6)], collapse=' '),
				paste(sprintf('%.1f',current_patch_info[c(8,9,2)]), collapse=' '),
				sprintf('%.2f',1.0),
				sprintf('%.2f', ifelse(!is.na(current_patch_info[19]),current_patch_info[19],0)*BASEMENT_DEPTH + ifelse(!is.na(current_patch_info[15]),current_patch_info[15],0)*PAVEDROAD_DEPTH + ifelse(!is.na(current_patch_info[20]),current_patch_info[20],0)*PAVEDROAD_DEPTH),
				drainage_type, 
				total_gamma, normal_neighborNum, '\n', file=surfaceflow_table_buff,sep=' ') 
				
			cat( paste(
				allNeighbourInfo[1,],
				allNeighbourInfo[2,],
				allNeighbourInfo[3,], 
				sprintf('%.5f',normal_neighbor_frac_gamma),
				sprintf('%.2f',allNeighbourInfo[7,]/allNeighbourInfo[4,]),
				sprintf('%.2f',allNeighbourInfo[7,]),sep=' '), file=surfaceflow_table_buff,sep='\n')	
			
			# current_patch_info[15]>0	road / storm drain
			if(current_patch_info[15]>0 & current_patch_info[21]==0) cat(
				patch_info_lowest[1], 
				patch_info_lowest[7], 
				patch_info_lowest[6], 
				sprintf('%.5f', stormsurfacedrainFrac['road']),
				sprintf('%.2f',1.0),
				sprintf('%.2f',1.0),'\n', file=surfaceflow_table_buff,sep=' ')	
				
			# current_patch_info[19]>0 roof
			if(current_patch_info[19]>0 & current_patch_info[21]==0) cat(
				patch_info_lowest[1], 
				patch_info_lowest[7], 
				patch_info_lowest[6], 
				sprintf('%.5f', stormsurfacedrainFrac['roof']),
				sprintf('%.2f',1.0),
				sprintf('%.2f',1.0),'\n', file=surfaceflow_table_buff,sep=' ')	
				
			# current_patch_info[20]>0 parking
			if(current_patch_info[20]>0 & current_patch_info[21]==0) cat(
				patch_info_lowest[1], 
				patch_info_lowest[7], 
				patch_info_lowest[6], 
				sprintf('%.5f', stormsurfacedrainFrac['parking']),
				sprintf('%.2f',1.0),
				sprintf('%.2f',1.0),'\n', file=surfaceflow_table_buff,sep=' ')
				
			# current_patch_info[17]>0 & current_patch_info[21]==0
			if(stormsurfacedrainFrac['extenddrain']>0) cat(
				patch_info_lowest[1], 
				patch_info_lowest[7], 
				patch_info_lowest[6], 
				sprintf('%.5f', stormsurfacedrainFrac['extenddrain']),
				sprintf('%.2f',1.0),
				sprintf('%.2f',1.0),'\n', file=surfaceflow_table_buff,sep=' ')
				
			if(current_patch_info[22]>0 & current_patch_info[11]==0) cat(
				patch_info_lowest[1], 
				patch_info_lowest[7], 
				patch_info_lowest[6], 
				sprintf('%.5f',0.6),
				sprintf('%.2f',1.0),
				sprintf('%.2f',1.0),'\n', file=surfaceflow_table_buff,sep=' ')			
				
				
		}else{
			# same as subsurface flow
			
			normal_neighbor_frac_gamma = neighbor_frac_gamma
			
			## stop routing from roof to roof on surface
			normal_neighbor_frac_gammaSUM = sum(normal_neighbor_frac_gamma)
			normal_neighbor_frac_gamma = normal_neighbor_frac_gamma * (1.0 - allNeighbourInfo[8,]) # roofFrac
			if(sum(normal_neighbor_frac_gamma)>0){normal_neighbor_frac_gamma = normal_neighbor_frac_gamma/sum(normal_neighbor_frac_gamma)*normal_neighbor_frac_gammaSUM;}else{ normal_neighbor_frac_gamma = rep(0,length(normal_neighbor_frac_gamma)) }
			
			cat(
				paste(current_patch_info[c(1,7,6)], collapse=' '),
				paste(sprintf('%.1f',current_patch_info[c(8,9,2)]), collapse=' '),
				sprintf('%.2f',1.0),
				sprintf('%.2f', ifelse(!is.na(current_patch_info[19]),current_patch_info[19],0)*BASEMENT_DEPTH + ifelse(!is.na(current_patch_info[15]),current_patch_info[15],0)*PAVEDROAD_DEPTH + ifelse(!is.na(current_patch_info[20]),current_patch_info[20],0)*PAVEDROAD_DEPTH),
				drainage_type, 
				total_gamma,length(withinNeighbourRC),'\n', file=surfaceflow_table_buff,sep=' ') 
				
			cat( paste(
				allNeighbourInfo[1,],
				allNeighbourInfo[2,],
				allNeighbourInfo[3,], 
				sprintf('%.5f', normal_neighbor_frac_gamma),
				sprintf('%.2f',allNeighbourInfo[7,]/allNeighbourInfo[4,]),
				sprintf('%.2f',allNeighbourInfo[7,]),sep=' '), file=surfaceflow_table_buff,sep='\n')
		}# 
		
		
		#return <- 1;
	#})#tapply
	}# for loop ii 
	close(subsurfaceflow_table_buff)
	close(surfaceflow_table_buff)
	
	## part 3: write out
	#write(subsurfaceflow_table_buff, arg[1], ncolumns=1)
		
		
