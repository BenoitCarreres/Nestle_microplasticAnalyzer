###############################################################################
#### Author: Benoit Carreres
###############################################################################
############             Nestle microplastic analyzer             #############
# Tool developed by Société des Produits Nestlé S.A. (Nestlé Research) to process Raman spectra (Horiba) data and identify microplastic particles.
# This tool was developed as part of a larger application to identify microplastic raman spectra. This minimalized version was generated to allow reproducibility of our results published in Nature Scientific Reports. The necessary data and random forest database are available on Zenodo (link bellow). The R script should take an input csv file that will provide the necessary information to analyzed the published dataset.
# Users can also modify this code to perform other analysis of their own HORIBA generated *.spc files.
#
# This software was built and tested on R 4.0.2
#
# Article available in Nature Scientific Reports: -link-
# Data and random forest model at Zenodo: -link-
###############################################################################
# This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
###############################################################################

library("hyperSpec")
library("ranger")
library("splines")
library("igraph")

label <<- list()## Fix a bug in the required version of hyperSpec

read.spc <- function(file, classes=NULL){
  message("reading spc file: ",file)
  spc = suppressWarnings(hyperSpec::read.spc(file))
  message("reading done, preprocessing")
  spc = spc.dataCorrection(spc)
  #
  if(!is.null(classes)){
    spc$predictions = classes
  }else if (is.null(spc$predictions) ){
    spc$predictions = "NMP"
  }
  message("preprocessing done, adding plate positions")
  spc = addPlatePositions(spc)
  message("adding plate positions done, exiting read.spc.list")
  if(hyperSpec::chk.hy(spc)) {
      return(spc)
  } else {
      stop("SPC object is not valid, please check the function for errors")
  }
}
######################################################
addPlatePositions <- function(spc, roundPositions=FALSE) {
  ### this cannot be done on a splited spc object!
  if(is.null(spc@data$z)) spc@data$z=0
  if(is.null(spc@data$z.end)) spc@data$z.end=0
  if(is.null(spc@data$w)) spc@data$w=0
  #
  spc@data$cols = as.numeric(factor(spc@data$z))
  spc@data$rows = as.numeric(factor(spc@data$w))
  #
  if((max(spc$cols) * max(spc$rows)) != nrow(spc)) {
    warning("There is a missmatch in #row and #col.")
    rowToRm = which(table(spc$rows) != median(table(spc$rows)))
    colToRm = which(table(spc$cols) != median(table(spc$cols)))
    if(length(rowToRm)>0 & length(colToRm)>0){
      if(length(rowToRm) < length(colToRm)){### if p with one dim, the other has many issues (=number of more/missing measurements)
        rowToRm = which(table(spc$rows) != median(table(spc$rows)))### Usually the problem is less and the latest row or col.
        warning("Cleaning up row: ", rowToRm)
        spc = spc[which(spc$rows!=rowToRm),]
      } else {
        colToRm = which(table(spc$cols) != median(table(spc$cols)))
        warning("Cleaning up col: ", colToRm)
        spc = spc[which(spc$cols!=colToRm),]
      }
    }
  }
  #
  if(hyperSpec::chk.hy(spc)) {
    return(spc)
  } else {
    stop("SPC object is not valid, please check the function for errors")
  }
}

spc.dataCorrection <- function(spc){
  spc = spc.spline(spc)
  spc@data$spc.raw = spc@data$spc
  spc = spc.bl_fp(spc)
  spc = hyperSpec::spc.smooth.spline(spc, all.knots=TRUE, spar=0.4)
  spc = spc.scale(spc)
  return(spc)
}
spc.bl_fp <- function(spc, lineType="cor") {
  mat = spc@data$spc
  bl = baseline::baseline.fillPeaks(mat, lambda=0.5, hwi=10, it=3, int=100)
  if(lineType=="cor") {
    mat = bl$corrected ###
  } else if(lineType=="bl") {
    mat = bl$baseline
  }
  spc@data$spc = mat
  return(spc)
}
spc.spline = function(spc){
  if(identical(spc@wavelength,newx)) return(spc)
  newspc <- t(apply(spc@data$spc, 1, .spline, x = spc@wavelength, newx = newx))
  if (any(is.na(newspc))){
    msg = paste0("NAs generated during SPLINE correction. Are all the values within the analyzed range? ", min(newx)," - ",max(newx),"cm-1.")
    warning(msg)
  }
  spc@data$spc <- newspc
  spc@wavelength <- newx
  colnames(spc@data$spc) = spc@wavelength
  validObject(spc)
  spc
}
.spline <- function(x, y, newx) {
    s=splinefun(x, y, method="monoH.FC")
    s(newx)
}
spc.scale <- function(spc) {
  spc@data$spc = t(scale(t(spc@data$spc), center = FALSE, scale = TRUE))
  if(!is.null(spc$scaleFactor)){
    message("scaleFactor found, multiplying new factor to existing one")
    newFactor = attr(spc@data$spc, 'scaled:scale')
    spc$scaleFactor = spc$scaleFactor * newFactor
  } else {
    spc$scaleFactor = attr(spc@data$spc, 'scaled:scale')
  }
  message("removing attribute: scaled:scale")
  attr(spc@data$spc, 'scaled:scale') = NULL
  return(spc)
}
##
spc.addParticleStats <- function(spc, maxDist=2, by=NULL, LQthrs=NULL){
  if(is.null(LQthrs) && !is.null(attr(spc@data$predictions, "LQthrs"))) LQthrs = attr(spc@data$predictions, "LQthrs")
  if(is.null(LQthrs)) LQthrs = sapply(sort(unique(spc@data$predictions))[!sort(unique(spc@data$predictions))=="NMP"], function(x) 0.4)
  ### Setting the new grouping by attribute value
  by = "Area_um2"### If nothing, default. this case sould not happen though

  ### Different breaks for different column to make the grouping (by)
  sizeBreaks = read.table("_SizeCountBreaks.txt", sep="=")
  message("Declared break points for: ",paste(sizeBreaks[,1], collapse=", "))
  if(!any(sizeBreaks[,1] == by)) {
    warning("/!\\ Specified breaks not found, using hardcoded default break points.\nSize counts by:", by)
    breaks = c(0,5,10,20,50,100,500,Inf)## This should no occur unless the file has been broken
  } else {
    breaks.str = sizeBreaks[sizeBreaks[,1] == by,2]
    breaks = as.numeric(unlist(strsplit(breaks.str, split=",", fixed=TRUE)))
    breaks = c(0,breaks,Inf)
    # message("Breaks = ",paste(breaks, collapse=", "))
  }
  ###
  ### Consider mapping data
  ###### Different distance measurement: see Chebyshev distance
  if((max(spc$cols) * max(spc$rows)) != nrow(spc)) {
    stop("There is a missmatch in #row and #col. This should have been corrected earlier in the process")
  }
  spc = spc[order(as.numeric(paste0(sprintf("%04d",spc$rows),sprintf("%04d",spc$cols)))),]
  ## require(igraph)
  # g <- graph.lattice(c(length(unique(spc@data$rows)),length(unique(spc@data$cols))), nei=maxDist)### set to 2, allows diagonals and jump by 1 block (no jump in diagonal!).
  g <- graph.lattice(c(length(unique(spc@data$cols)),length(unique(spc@data$rows))), nei=maxDist)
  edgelist <- get.edgelist(g)
  retain <- spc@data$predictions[edgelist[,1]] == spc@data$predictions[edgelist[,2]]
  g <- delete.edges(g, E(g)[!retain])
  #### lyt <- layout.auto(g) #To visualize
  #### plot(g, layout=lyt) #To visualize
  # Now calculate the stats
  spc@data$particleID = clusters(g)$membership
  spc@data$particleID[spc@data$predictions == "NMP"] <- 0
  spc@data$particleID[spc@data$particleID != 0] = as.numeric(factor(spc@data$particleID[spc@data$particleID != 0]))
  spc@data$particleID = paste0("#",sprintf("%06d",spc@data$particleID))### worked with integers, but paulo wanted #0001 format
  df=data.frame(table(spc@data$particleID)); colnames(df)=c("particleID","count")
  df$Area_um2 = df$count * spc.measurementArea(spc)
  df$meanScore = round(sapply(split(spc@data$predictionsScore, spc@data$particleID),mean),digits=2)
  ##########
  particles = merge(unique(data.frame(particleID=spc@data$particleID,type=spc@data$predictions)), df)
  particles = particles[particles$particleID!="#000000",]
  particles$type = factor(particles$type, levels=names(LQthrs))### Else NMP is still there and gives pb
  numParticles = table(particles$type, exclude = "NMP")### I think this is useless at this point
  particles.type.surfaceGroups = sapply(split.data.frame(particles, as.vector(particles$type[particles$particleID!="#000000"])), function(x) table(cut(x[,by], breaks=breaks, dig.lab=10)))
  if(is.null(dim(particles.type.surfaceGroups)) & dim(numParticles)==0){
    particles.type.surfaceGroups = setNames(data.frame(matrix(ncol=ncol(spc@data$predictionTable),nrow = 1)), colnames(spc@data$predictionTable))
    particles.type.surfaceGroups[1,] = 0
    rownames(particles.type.surfaceGroups) = "total"
  } else {
    particles.type.surfaceGroups = rbind(particles.type.surfaceGroups, total=numParticles)
  }
  attr(spc@data$predictions, "particles") <- particles
  attr(spc@data$predictions, "particlesStats") <- particles.type.surfaceGroups
  attr(spc@data$predictions, "LQthrs") <- LQthrs
  attr(spc@data$predictions, "groupBy") <- by
  return(spc)
}
#
spc.rfPrediction <- function(spc, rf) {
  spc = spc.checkAndFixHeaders(spc)### fixing here, because this is the place where it can bother.
  vals = spc2cdata(spc,groups=NULL)
  if(any(!colnames(vals) %in% paste0("X",newx))){### Non existing values are problematic, missing values should be fine.
    warning("Possible spectral value missmatch, most-probably header problem in spc@data$spc")
  }
  pred = predict(rf, vals, num.threads= parallel::detectCores())
  spc = rf.probs.postprocess(pred, spc)### Call function that will add info to the spc object, wether it is probs table or a best fit category
  return(spc)
}
spc.checkAndFixHeaders <- function(spc) {
  if(is.null(colnames(spc@data$spc))) {
    warning("Correcting spc@data$spc header names missing (unkown reason, possibly spc file issue)")
    colnames(spc@data$spc) = spc@wavelength
  }
  if(is.null(colnames(spc@data$spc.raw))) {
    warning("Correcting spc@data$spc.raw header names missing (unkown reason, possibly spc file issue)")
    colnames(spc@data$spc.raw) = spc@wavelength
  }
  return(spc)
}
spc2cdata <- function(spc, groups="predictions", allData=FALSE){
  ##  use groups=NULL if no metadata is wanted (just want the spectral data)
  mat = as.matrix(spc)
  # colnames(mat) = paste0(colnames(mat),"_cm-1")### not cecessary, X will be added when creating a DF
  df = data.frame(mat)
  if(!allData) {
    cdata = cbind(type = spc@data[,groups], df)
  } else if(allData) {
    if(is.list(spc$predictionScore)) spc$predictionScore = unlist(spc$predictionScore)
    if (!is.null(spc$spc.raw)) {###normally alway there
      mat.raw = as.matrix(spc$spc.raw)## If raw data is wanted. I do not think it is necessary.
      df.raw = data.frame(mat.raw)
      colnames(df.raw)=paste0("raw.",colnames(df.raw))
      df=cbind(df.raw,df)
    }
    if (!is.null(spc$predictionTable)) {
      df.pt=spc$predictionTable
      colnames(df.pt)=paste0("pred.",colnames(df.pt))
      df = cbind(df.pt,df)
    }
    cdata = cbind(spc@data[!names(spc@data) %in% c("spc","spc.raw","predictionTable")], df)#, df.raw, )
  }
  return(cdata)
}
rf.probs.postprocess <- function(predictions, spc) {
  if(class(predictions$predictions) == "matrix") {
    message("postprocessing RF Probabilities results")
    spc$predictionTable = predictions$predictions
    preds = colnames(spc$predictionTable)[apply(spc$predictionTable, 1, which.max)]
    # preds[preds == "u1"] = "NMP" old stuff for testing special group
    spc$predictions = preds
    spc$predictionsScore = apply(spc$predictionTable, 1, max)
  } else {
    message("postprocessing RF Classification results")
    spc$predictions = predictions$predictions
  }
  return(spc)
}
spc.measurementArea <- function(spc) {
  ### possible improvement abs(z.end:z)^2
  ### report z and w (interspaces) to the reports?
  Zs = sort(unique(spc$z))
  Ws = sort(unique(spc$w))
  z = round(abs(Zs[2]-Zs[1]))
  w = round(abs(Ws[2]-Ws[1]))
  z*w
}

################################################################################
newx = seq(559,1990,3)###closest to the number of measurements
rfFile = "DB_V2.10.RF.rds"
rf.ref = readRDS(rfFile)
# Spiked milk samples
milkThresFile="MilkSpiked_thresholds.csv"
waterThresFile="WaterSpiked_thresholds.csv"
analyseFileList <- function(thresFile) {
  pref = gsub("_thresholds.csv","",gsub(".*/","",thresFile))
  folder = gsub("[^/]*$","",thresFile)
  thresTable = read.csv(thresFile)
  countsTable = apply(thresTable,1,function(line){
    fileName=line["file"]
    file=paste0(folder, pref, "/", gsub(" .+$","",fileName), "/SPC/",fileName)
    if(!file.exists(file)) {
      message("NOT found, skipping: ", file)
      return()
    }
    spc = read.spc(file)
    spc@data$filename = file
    rownames(spc) = paste0("spl_",rownames(spc))### better to keep it here in case the import functions of spc becomes same for ref and sample
    # Use RF model to classify signals
    spc = spc.rfPrediction(spc, rf.ref)
    spc$rf_model = rfFile
    #
    spc = spc.addParticleStats(spc)
    LQthrs = attr(spc@data$predictions, "LQthrs")### automaticall set by the stat function, here collecting it just in case, need to be given to the addParticleStats when modified
    LQthrs$PA = line["PA"]
    LQthrs$PE = line["PE"]
    LQthrs$PMMA = line["PMMA"]
    LQthrs$PP = line["PP"]
    LQthrs$PS = line["PS"]

    spc.mp = spc[spc$predictions != "NMP",]
    spc.mp$predictions = factor(spc.mp$predictions, levels=names(LQthrs))### so table() returns zero counts for missing values
    parts = attr(spc@data$predictions, "particles")
    HQP = table(parts$type[parts$meanScore >= LQthrs[parts$type]])
    HQP = HQP[names(HQP) %in% c("PA","PE","PMMA","PP","PS")]
    print(HQP)
    return(HQP)
  })
  countsTable = t(countsTable)## apply returns the table relatively transposed
  rownames(countsTable) = thresTable[,1]
  write.csv(countsTable, file=paste0(folder,pref,"_counts.csv"))
}
########## Run analysis #########
# Spiked water samples
analyseFileList(milkThresFile)
analyseFileList(waterThresFile)

#
# # save object to rds file
# fileName.rds = gsub("\\.spc$",".rds",input$UploadSample_spc$name[fileNum])
# # saveRDS(spc, fileName.rds)
