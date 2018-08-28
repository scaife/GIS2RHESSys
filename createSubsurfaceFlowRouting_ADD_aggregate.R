## flow table
options(scipen=999)
arg=commandArgs(T)
library(rgrass7)
library(rgdal)
gis = gmeta()

DtoR = pi/180
RtoD = 1/DtoR



	flowtable = read.table('flows/subsurfaceflow.txt',fill=T,skip=1)
	cond = !is.na(flowtable[,11])
	patchID = flowtable[cond,1]
	
	
	rast5 = readRAST(c('patch','aggforest','agglawn'))	
		patch = rast5@data[[1]][mask]
		aggforest = rast5@data[[2]][mask]
		agglawn = rast5@data[[3]][mask]
		
		
	matchOrder = match(patchID, patch)
	aggforestCode = aggforest[matchOrder] + 0; aggforestCode[is.na(aggforestCode)]=0
	agglawnCode = agglawn[matchOrder] + 1; agglawnCode[is.na(agglawnCode)]=0
	aggCode = aggforestCode + agglawnCode
			
	flowtable$aggID = rep(NA,dim(flowtable)[1])
	flowtable$aggIDLOC = rep(NA,dim(flowtable)[1])	
		
	flowtable$aggID[cond] = aggCode
	flowtable$aggIDLOC[cond] = aggCode
	
	write.table(flowtable,'flows/subsurfaceflow_agg.txt',row.names=F,col.names=F,quote=F, na='')
	