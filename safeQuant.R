#!/usr/bin/Rscript

# Author: ahrnee-adm
# Adaptation: Georgia Angelidou (v.2.4)
###############################################################

############################################################### INIT ###############################################################
#### DEPENDANCIES

suppressWarnings(suppressPackageStartupMessages(library("stringr", quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library("affy", quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library("limma", quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(gplots, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(seqinr, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(corrplot, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(optparse, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(data.table, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(ggplot2, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(ggrepel, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(readxl, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(hash, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(scales, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(pcaMethods, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(colorspace, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(gridExtra, quiet=T)))
suppressWarnings(suppressPackageStartupMessages(library(crayon, quiet=T)))

# Note: Don't add the / or \ after the R folder name. Otherwise it will not access the correct path
#sourceDirOSX <- ""
sourceDirOSX <- "D:/proteomics/Github/SafeQuant_p/R"
sourceDirTPP <-  "C:/Proteomics/SafeQuant-master_Geo/R"

# first check if dev or tpp mode
if(file.exists(sourceDirOSX) | file.exists(sourceDirTPP)){

	sourceDir <- ifelse(file.exists(sourceDirOSX),sourceDirOSX,sourceDirTPP)

	source(paste(sourceDir,"/ExpressionAnalysis.R",sep=""))
	source(paste(sourceDir,"/SafeQuantAnalysis.R",sep=""))
	source(paste(sourceDir,"/Graphics.R",sep=""))
	source(paste(sourceDir,"/GGGraphics.R",sep=""))
	source(paste(sourceDir,"/IdentificationAnalysis.R",sep=""))
	source(paste(sourceDir,"/Parser.R",sep=""))
	source(paste(sourceDir,"/TMT.R",sep=""))
	source(paste(sourceDir,"/UserOptions.R",sep=""))
	source(paste(sourceDir,"/OutputTableFormat.R",sep=""))

}else if("SafeQuant" %in%  installed.packages()[,1]){ # used installed SafeQuant

	cat("Loading SafeQuant Library \n")
	library("SafeQuant")

}else{
	stop("SafeQuant Package not installed\n")
}

rm(sourceDirOSX, sourceDirTPP)

VERSION <- "2.4"

### USER CMD LINE OPTIONS
userOptions <- getUserOptions(version=VERSION)

### USER CMD LINE OPTIONS END

if(userOptions$verbose) print(userOptions$proteinQuant)

### SUPRESS WARNINGS
if(!userOptions$verbose){
	options(warn=-1)
}

#### DEPENDENCIES

if(userOptions$verbose) print(userOptions)

############################################################### PARSING ###############################################################

if(userOptions$verbose) cat("PARSING INPUT FILE \n")

# get file type
fileType <- .getFileType(userOptions$inputFile)

### Progenesis Export
if(fileType %in% c("ProgenesisProtein","ProgenesisFeature","ProgenesisPeptide")){

	# default
	expDesign <- getExpDesignProgenesisCsv(userOptions$inputFile)

	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,expDesign)
	}

	if(fileType == "ProgenesisProtein"){
		cat("INFO: PARSING PROGENESIS PROTEIN EXPORT FILE ",userOptions$inputFile, "\n" )
		eset <- parseProgenesisProteinCsv(file=userOptions$inputFile,expDesign=expDesign)

	}else if(fileType == "ProgenesisPeptide"){

		#"ProgenesisFeature"
		cat("INFO: PARSING PROGENESIS PEPTIDE EXPORT FILE ",userOptions$inputFile, "\n" )
	  # gives some general information about the different files all together
		eset <- parseProgenesisPeptideMeasurementCsv(file=userOptions$inputFile,expDesign=expDesign, exclusivePeptides=userOptions$FExclusivePeptides)

	}else{ 	#"ProgenesisFeature"

		cat("INFO: PARSING PROGENESIS FEATURE EXPORT FILE ",userOptions$inputFile, "\n" )
		eset <- parseProgenesisFeatureCsv(file=userOptions$inputFile,expDesign=expDesign)
	}

# Scaffold Export (TMT data)
}else if(fileType == "ScaffoldTMT"){
	cat("INFO: PARSING SCAFFOLD RAW EXPORT FILE ",userOptions$inputFile, "\n" )

	# get default experimental design
	# six plex or ten plex ?
	# use default experimental design unless specified by the user
	nbPlex <- .getNbPlex(userOptions$inputFile)
	if(nbPlex == 6){
		# 6-plex default: 1,2,3:4,5,6
		expDesign <- data.frame(condition=paste("Condition",sort(rep(c(1,2),3)),sep=""),isControl=sort(rep(c(T,F),3),decreasing=T) )
	}else{
		# 10-plex default is "1,4,7,10:2,5,8:3,6,9"
		expDesign <- data.frame(condition=paste("Condition",c(1,2,3,1,2,3,1,2,3,1),sep=""),isControl=c(T,F,F,T,F,F,T,F,F,T) )
	}

	eset <- parseScaffoldRawFile(file=userOptions$inputFile,expDesign=expDesign)

	if(userOptions$TAdjustRatios){
		#if((nbPlex == 10)){ # only possible for tmt-10 plex
			nbCalMixSpectra <- sum( (fData(eset)$proteinName %in% names(CALIBMIXRATIOS)))
			if(nbCalMixSpectra < 100) stop("Not enough Calibration Mix spectra were found: ",nbCalMixSpectra, "\n ")
			cat("INFO: FOUND  ", nbCalMixSpectra ," Calibration Mix spectra\n")
			esetCalibMix <- .getCalibMixEset(eset)

			# discard calibration mix proteins
			eset <- eset[!(fData(eset)$proteinName %in% names(CALIBMIXRATIOS)),]

			intAdjObj <- .intensityAdjustment(eset, esetCalibMix)

		#}else{
		#	stop("Ratio Correction Not implemented for TMT 6-plex")
		#}
	}

	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,expDesign)
	}

#	eset <- parseScaffoldRawFile(file=userOptions$inputFile,expDesign=expDesign)
	# apply specified experimental design
	eset <- createExpressionDataset(expressionMatrix=exprs(eset)[,rownames(expDesign)],expDesign=expDesign,featureAnnotations=fData(eset))

	if(!is.na(userOptions$scaffoldPTMSpectrumReportFile)){

		cat("INFO: ADDING SCAFFOLD PTM ANNOTATIONS \n")
		eset <- addScaffoldPTMFAnnotations(eset,userOptions$scaffoldPTMSpectrumReportFile)

	}

}else if(fileType == "MaxQuantProteinGroup"){

	# get user specified experimental design
	if(!is.na(userOptions$expDesignTag)){
		# user specified
		expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,data.frame(condition=paste("Condition",1:1000),isControl=c(F,1000)), file=userOptions$inputFile, fileT =fileType)
		}else{
		stop(red("Please Specify Experimental Design"))
		}
  cat("INFO: PARSING MaxQuant PROTEIN EXPORT FILE ",userOptions$inputFile, "\n" )
	eset <- parseMaxQuantProteinGroupTxt(userOptions$inputFile,expDesign=expDesign, method="auc", sr_flag = userOptions$SRawDataAnalysis)


}else if(fileType == "GenericCSV"){

}else if (fileType == "SpectronautProteinGroup"){
  if(!is.na(userOptions$expDesignTag)){
    expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,data.frame(condition=paste("Condition",1:1000),isControl=c(F,1000)), file=userOptions$inputFile, fileT=fileType)
  }else{
    stop(red("Please Specify Experimental Design File"))
  }

  cat("INFO: PARSING Spectronaut PROTEIN EXPORT FILE ",userOptions$inputFile, "\n" )
  eset <- parseSpectronautProteinGroupTxt(userOptions$inputFile,expDesign=expDesign, method="auc")
}else if (fileType == "DiaNNProteinGroup"){
  if(!is.na(userOptions$expDesignTag) && str_detect(userOptions$expDesignTag, ".txt$")){
    expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,data.frame(condition=paste("Condition",1:1000),isControl=c(F,1000)), file=userOptions$inputFile, fileT=fileType)
  }else{
    stop(red("Please Specify Experimental Design File"))
  }

  cat("INFO: PARSING Dia-NN PROTEIN EXPORT FILE ",userOptions$inputFile, "\n" )
  eset <- parseDiaNNProteinGroupTxt(userOptions$inputFile,expDesign=expDesign, method="auc")

}else if (fileType == "DIANN_Peptide"){
  expDesign <- expDesignTagToExpDesign(userOptions$expDesignTag,data.frame(condition=paste("Condition",1:1000),isControl=c(F,1000)), file=userOptions$inputFile, fileT=fileType)
  cat("INFO: PARSING Dia-NN PEPTIDE EXPORT FILE ",userOptions$inputFile, "\n" )
  eset <- parseDIANNPeptideMeasurementCsv(file=userOptions$inputFile,expDesign=expDesign, exclusivePeptides=userOptions$FExclusivePeptides)

}else{
	stop("Unknown File Type", userOptions$inputFile)
}


# test option, limit number of entries
if(userOptions$test){
	eset <- eset[1:nrow(eset) %in% sample(1:min(300,c(nrow(eset))),replace=F),]
}

# parse .fasta file
if(!is.na(userOptions$proteinFastaFile)){
	cat("INFO: PARSING PROTEIN SEQUENCE DB ",userOptions$proteinFastaFile, "\n" )
	### read protein db
	proteinDB <- read.fasta(userOptions$proteinFastaFile,seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

	# dirty fix check if ACs in Progenesis file are stripped
	if(isStrippedACs(sample(fData(eset)$proteinName,100))){
		cat("INFO: RE-FORMATTING ACCESSION NUMBERS\n")
		names(proteinDB) <- stripACs(names(proteinDB))
	}
}

############################################################### CREATE DATA MODEL ###############################################################

if(userOptions$verbose) print(eset)
if(userOptions$verbose) print(pData(eset))
if(userOptions$verbose) print(names(fData(eset)))

#### CREATE FEATURE DATA AND FILTER (pre-rollup)

# generic
if (fileType %in% c("DiaNNProteinGroup", "DIANN_Peptide")){
  filter <- data.frame(
    con=isCon(fData(eset)$ac)	# contaminants
    ,ac = !(grepl(userOptions$selectedProteinName,fData(eset)$proteinName,ignore.case=T)) # protein ac
  )
}else{
  filter <- data.frame(
		con=isCon(fData(eset)$proteinName)	# contaminants
		,ac = !(grepl(userOptions$selectedProteinName,fData(eset)$proteinName,ignore.case=T)) # protein ac
  )
}

# do not filter TMT data
if("pMassError" %in% names(fData(eset))  &&  (fileType != "ScaffoldTMT") ){
	### applicable to Progenesis feature Exports
	if(is.na(userOptions$precursorMassFilter)){ # if not user specified
		# automatically get precursor limits, X * sd of 50% top scoring
		userOptions$precursorMassFilter <- getMeanCenteredRange(fData(eset)$pMassError[fData(eset)$idScore > quantile(fData(eset)$idScore,na.rm=T)[3]],nbSd = 3)
		filter <- cbind(filter, pMassError=
						(fData(eset)$pMassError < userOptions$precursorMassFilter[1])
						| (fData(eset)$pMassError > userOptions$precursorMassFilter[2]) # precursor mass tolerance
		)
	}
}

if("ptm" %in% names(fData(eset))){

	# add motif-X and ptm coordinates
	if(exists("proteinDB")){

		cat("INFO: EXTRACTING PTM COORDINATES AND MOTIFS\n")
		#format 1) progensis  2) scaffold
		eset <- .addPTMCoord(eset,proteinDB,motifLength=6, isProgressBar=T,format= (fileType == "ScaffoldTMT") +1)

	}
	filter <- cbind(filter
			, ptm = !(grepl(userOptions$selectedModifName,as.character(fData(eset)$ptm),ignore.case=T))
			, nbPtmsPerPeptide = (fData(eset)$nbPtmsPerPeptide > userOptions$maxNbPtmsPerPeptide) )

}

if("peptide" %in% names(fData(eset))){
	filter <- cbind(filter
			, peptideLength =nchar(as.character(fData(eset)$peptide)) < userOptions$minPeptideLength
			, charge =  fData(eset)$charge == 1 # discard singly charged
	)
}

if(!("nbPeptides" %in% names(fData(eset))) & fileType %in% c("ProgenesisPeptide", "DIANN_Peptide")){
  ### set nb peptides per protein
  eset <- setNbPeptidesPerProtein(eset)
}
if(("nbPeptides" %in% names(fData(eset)))){
  filter <- cbind(filter,nbPeptides=(fData(eset)$nbPeptides < userOptions$minNbPeptidesPerProt))
}


# do not filter TMT data
#if(("idScore" %in% names(fData(eset))) && (fileType != "ScaffoldTMT")){

if(("idScore" %in% names(fData(eset)))){
	eset <- addIdQvalues(eset)
	filter <- cbind(filter,qvalue=fData(eset)$idQValue > userOptions$fdrCutoff)
}

# set pre-rollup filters
eset <- .setFilter(eset,filter=filter)
rm(filter)


### make sure at least 1 feature pass the filter
if(sum(!fData(eset)$isFiltered,na.rm=T) == 0){
	stop("CHECK FILTER SETTINGS. ALL FEATURES WERE FILTERED OUT")
}

#### CREATE FEATURE DATA AND FILTER END
### SET ANCHOR PROTEINS
fData(eset)$isNormAnchor <- grepl(userOptions$normAC,fData(eset)$proteinName)

if(userOptions$verbose){
	cat("\nNB. ANCHOR PROTEINS: ")
	cat(sum(fData(eset)$isNormAnchor))
	cat("\n")
	print(fData(eset)$proteinName[fData(eset)$isNormAnchor])
	cat("\n")
}

### SET ANCHOR PROTEINS END

############################################################### EXPRESSION ANALYSIS ###############################################################
# create paired experiemntal design()
if(userOptions$ECorrelatedSamples){
	eset <- createPairedExpDesign(eset)
}
# add filters etc to adjusted expressionSet
# update expDesign of intAdjObj$esetAd
if(exists("intAdjObj")){
	fData(intAdjObj$esetAdj) <- fData(eset)
	pData(intAdjObj$esetAdj) <- pData(eset)
	exprs(intAdjObj$esetAdj) <- exprs(intAdjObj$esetAdj)[,colnames(exprs(eset))]
}
### non-pairwise stat test
statMethod <- c("")
if(userOptions$SNonPairWiseStatTest) statMethod <- c("all") #@TODO what about 'naRep'
if(userOptions$SRawDataAnalysis){ # No Normalization
	esetNorm <- eset
	if(exists("intAdjObj")) intAdjObj$esetAdjNorm <- intAdjObj$esetAdj
}else{ # NORMALIZE
	method <- c("global","median")
	# norm based on sum if norm anchor is specified
	if(sum(fData(eset)$isNormAnchor) < nrow(eset)) method <- c("global","sum")

	esetNorm <- sqNormalize(eset, method=method)
	if(exists("intAdjObj")) intAdjObj$esetAdjNorm <- sqNormalize(intAdjObj$esetAdj, method=method )

}

### MISSING VALUES IMPUTATION
# baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(esetNorm)[,1])),promille=5)
# exprs(esetNorm)[  is.na(exprs(esetNorm)) | (exprs(esetNorm) <= 0)  ] <- 0
# exprs(esetNorm) <- exprs(esetNorm) + baselineIntensity
# if(exists("intAdjObj")){
# 	baselineIntensity <- getBaselineIntensity(as.vector(unlist(exprs(intAdjObj$esetAdjNorm)[,1])),promille=5)
# 	exprs(intAdjObj$esetAdjNorm)[  is.na(exprs(intAdjObj$esetAdjNorm)) | (exprs(intAdjObj$esetAdjNorm) <= 0)  ] <- 0
# 	exprs(intAdjObj$esetAdjNorm) <- exprs(intAdjObj$esetAdjNorm) + baselineIntensity
# }
# ToDo: 18/05/2022
# esetNorm_only <- esetNorm
esetNorm = sqImpute(esetNorm,method=userOptions$SMissingValuesImutationMethod)

if(exists("intAdjObj")){
   intAdjObj$esetAdjNorm = sqImpute(intAdjObj$esetAdjNorm,method=userOptions$SMissingValuesImutationMethod )
}
# Note: this overwrites the option when the SNonPairWiseStatTest acivated and sets it as all
if (userOptions$medianInfo){
  statMethod <- c(statMethod, "median")
}else{
  statMethod <- c(statMethod, 'mean')
}

if((fileType == "ProgenesisProtein") |  (fileType == "MaxQuantProteinGroup") |
   (fileType == "SpectronautProteinGroup") | (fileType == "DiaNNProteinGroup")){
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered  | isDecoy(fData(esetNorm)$proteinName)
	if (userOptions$log2_i){
	  sqaProtein <- safeQuantAnalysis(esetNorm, method=statMethod, fcThrs=userOptions$ratioCutOff, log2v = FALSE)
	}else{
	  sqaProtein <- safeQuantAnalysis(esetNorm, method=statMethod, fcThrs=userOptions$ratioCutOff)
	}

	sqaProtein <- safeQuantAnalysis(esetNorm, method=statMethod, fcThrs=userOptions$ratioCutOff)
}else if((fileType == "ScaffoldTMT") && is.na(userOptions$scaffoldPTMSpectrumReportFile)){

	# roll-up protein level
	cat("INFO: ROLL-UP PROTEIN LEVEL\n")

	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered |  isDecoy(fData(esetNorm)$proteinName)
	# correct TMT ratios
	if(userOptions$TAdjustRatios){
		fData(intAdjObj$esetAdjNorm)$isFiltered <- fData(esetNorm)$isFiltered
		intAdjObjProt <- intAdjObj
		intAdjObjProt$esetAdjNorm <- rollUp(intAdjObj$esetAdjNorm,featureDataColumnName= c("proteinName"))
		sqaProtein <- safeQuantAnalysis(rollUp(esetNorm,featureDataColumnName= c("proteinName")), method=statMethod,intensityAdjustmentObj=intAdjObjProt, fcThrs=userOptions$ratioCutOff )
	}else{
		sqaProtein <- safeQuantAnalysis(rollUp(esetNorm,featureDataColumnName= c("proteinName")), method=statMethod , fcThrs=userOptions$ratioCutOff)
	}

	fData(sqaProtein$eset)$isFiltered <- fData(sqaProtein$eset)$isFiltered | isDecoy(fData(sqaProtein$eset)$proteinName) | (fData(sqaProtein$eset)$nbPeptides <  userOptions$minNbPeptidesPerProt)

}else{

	# roll-up peptide level
	cat("INFO: ROLL-UP PEPTIDE LEVEL\n")

	# correct TMT ratios
	if(userOptions$TAdjustRatios){
		cat("WARN: Ratio correction not yet implemented in this anlysis mode \n")
	}
  # Geo: here should add the method options and change from median to mean
	esetPeptide <- rollUp(esetNorm,featureDataColumnName= c("peptide","ptm")) # For progenesis
	# esetPeptide <- rollUp(esetNorm,featureDataColumnName= c("peptide"))  # For DIA-NN

  # esetPeptide_norm <- rollUp(esetNorm_only, featureDataColumnName= c("peptide"))
	# fdr filter
	# replace qValues by rollUp level qValues ()

	esetPeptide <- addIdQvalues(esetPeptide)
	# esetPeptide_norm <- addIdQvalues(esetPeptide_norm)
	if(fileType == "ScaffoldTMT"){
		fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | (fData(esetPeptide)$nbPeptides <  userOptions$minNbPeptidesPerProt)
	}else{
		# update filter to exclude peptide level hight qValues
		fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | (fData(esetPeptide)$idQValue > userOptions$fdrCutoff) | (fData(esetPeptide)$nbPeptides <  userOptions$minNbPeptidesPerProt)
		# fData(esetPeptide_norm)$isFiltered <- fData(esetPeptide_norm)$isFiltered | (fData(esetPeptide_norm)$idQValue > userOptions$fdrCutoff) | (fData(esetPeptide_norm)$nbPeptides <  userOptions$minNbPeptidesPerProt)
	}

	if(userOptions$proteinQuant){
		cat("INFO: ROLL-UP PROTEIN LEVEL\n")
		esetProtein <- rollUp(esetPeptide,featureDataColumnName= c("proteinName"))
		# esetProtein_norm <- rollUp(esetPeptide_norm,featureDataColumnName= c("proteinName"))
		# create one for the Row values

		esetProtein <- addIdQvalues(esetProtein)
		# esetProtein_norm <- addIdQvalues(esetProtein_norm)

		if(fileType == "ScaffoldTMT"){
			fData(esetProtein)$isFiltered <- fData(esetProtein)$isFiltered | isDecoy(fData(esetProtein)$proteinName) | (fData(esetProtein)$nbPeptides <  userOptions$minNbPeptidesPerProt)
		}else{
			fData(esetProtein)$isFiltered <- fData(esetProtein)$isFiltered | (fData(esetProtein)$idQValue > userOptions$fdrCutoff) | isDecoy(fData(esetProtein)$proteinName) | (fData(esetProtein)$nbPeptides <  userOptions$minNbPeptidesPerProt)
		}
		sqaProtein <- safeQuantAnalysis(esetProtein, method=statMethod, fcThrs=userOptions$ratioCutOff)
	}

	fData(esetPeptide)$isFiltered <- fData(esetPeptide)$isFiltered | isDecoy(fData(esetPeptide)$proteinName)
	sqaPeptide <- safeQuantAnalysis(esetPeptide, method=statMethod, fcThrs=userOptions$ratioCutOff)
	fData(esetNorm)$isFiltered <- fData(esetNorm)$isFiltered | isDecoy(fData(esetNorm)$proteinName) | (fData(esetNorm)$nbPeptides <  userOptions$minNbPeptidesPerProt)

	if(userOptions$top3 & userOptions$proteinQuant){
		cat("INFO: ROLL-UP TOP3\n")
		esetTop3 <-  rollUp(esetPeptide,featureDataColumnName= c("proteinName"), method="top3")
	}
}

### IBAQ
if(userOptions$iBAQ & userOptions$proteinQuant){
	cat("INFO: CALCULATING IBAQ VALUES\n")
	if(exists("proteinDB")){
		esetIBAQ <-  getIBAQEset(sqaProtein$eset, proteinDB=proteinDB)
	}else{
		cat("ERROR: proteinDB NOT FOUND NO iBAQ VALUES CALCULATED\n")
	}
}

### EXPRESSION ANALYSIS END

# all globalEnvironment
rm(addIdQvalues, addScaffoldPTMFAnnotations, createExpDesign, createExpressionDataset)

# I replace the 0 and 1 of the original store NA_IMP_CNT to the the actually values with NA when there was an imputated value
# Note: In the futere maybe this values should come after the original setting.
# I should also change the name but for know I will leave it like this
# if(!(exists("sqaPeptide"))){}

if(exists("sqaProtein")){
  expset <- exprs(sqaProtein$eset)
  if("rt" %in% statMethod || "quantile" %in% statMethod){
    ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("^NA_IMP_CNT", names(fData(sqaProtein$eset)))]
    ext_add_imp_col_repl <- sub("NA_IMP_CNT\\.", '', ext_add_imp_col)
    if ("file" %in% colnames(expDesign)){
      ext_add_imp_col_repl <- expDesign[ext_add_imp_col_repl, 'file']
      ext_add_imp_col_repl <- sub("Quantity\\.", "", ext_add_imp_col_repl)
    }
    if (userOptions$log2_i && userOptions$SRawDataAnalysis){
      if (fileType != "MaxQuantProteinGroup"){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Quantity", sep = ".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm.Quantity", sep = ".")
      }
    }else if (userOptions$log2_i){
      ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm.Quantity", sep = ".")
    }else if (userOptions$SRawDataAnalysis){
      if (fileType != "MaxQuantProteinGroup"){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Log2.Quantity", sep = ".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm_Log2.Quantity", sep = ".")
      }
    }else{
      ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm_Log2.Quantity", sep = ".")
    }
    ext_add_imp_col2 <- names(fData(sqaProtein$eset))[grepl("^NA_IMP_CNT", names(fData(sqaProtein$eset)))]
    na_cnt_info <- fData(sqaProtein$eset)[,ext_add_imp_col2]

    if (fileType == "DiaNNProteinGroup" || fileType == "MaxQuantProteinGroup"){
      na_cnt_info[na_cnt_info == 1] <- -1
      na_cnt_info[na_cnt_info == 0] <- NA
      na_cnt_info[is.na(na_cnt_info)] <- expset[is.na(na_cnt_info)]
      na_cnt_info[na_cnt_info == -1] <- NA
      fData(sqaProtein$eset)[, ext_add_imp_col_repl] <- na_cnt_info
    }
  }else{
    ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("NA_IMP_CNT", names(fData(sqaProtein$eset)))]
    ext_add_imp_col_repl <- sub("NA_IMP_CNT\\.", "", ext_add_imp_col)
    if ("file" %in% colnames(expDesign)){
      ext_add_imp_col_repl <- expDesign[ext_add_imp_col_repl, 'file']
      ext_add_imp_col_repl <- sub("Quantity\\.", "", ext_add_imp_col_repl)
    }
    if (userOptions$log2_i){
      if(fileType != "MaxQuantProteinGroup" && userOptions$SRawDataAnalysis){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Quantity", sep=".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm.Quantity", sep=".")
      }
    }else{
      if (fileType != "MaxQuantProteinGroup" && userOptions$SRawDataAnalysis){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Log2.Quantity", sep = ".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm_Log2.Quantity", sep = ".")
      }
    }
    ext_add_imp_col2 <- names(fData(sqaProtein$eset))[grepl("^NA_IMP_CNT", names(fData(sqaProtein$eset)))]
    na_cnt_info <- fData(sqaProtein$eset)[,ext_add_imp_col2]

    # if(fileType == "DiaNNProteinGroup" && fileType != "MaxQuantProteinGroup" && fileType != "SpectronautProteinGroup"){
    if(fileType == "DiaNNProteinGroup" || fileType == "ProgenesisPeptide"){
      na_cnt_info[na_cnt_info == 1] <- -1
      na_cnt_info[na_cnt_info == 0] <- NA
      na_cnt_info[is.na(na_cnt_info)] <- expset[is.na(na_cnt_info)]
      na_cnt_info[na_cnt_info == -1] <- NA
      # if (!(userOptions$log2_i)){
      #   na_cnt_info <- 2^na_cnt_info
      # }
      # Note: this should be activate if we want the no normalize - no imputated values
      # Note: needed for the future plots

      if (!(userOptions$SRawDataAnalysis)){
        gf <- pData(sqaProtein$eset)$globalNormFactors
        na_cnt_info2 <- sweep(na_cnt_info, 2, gf, FUN="/")

      }
      # if(!(userOptions$log2_i) && fileType != "ProgenesisPeptide"){
      #   # na_cnt_info <- log2(na_cnt_info)
      #   if(!(userOptions$SRawDataAnalysis)){
      #     na_cnt_info2 <- log2(na_cnt_info2)
      #   }
      # }
      fData(sqaProtein$eset)[,ext_add_imp_col_repl] <- na_cnt_info
      if (!(userOptions$SRawDataAnalysis)){
        colnames(na_cnt_info2) <- str_replace(colnames(na_cnt_info2), "NA_IMP_CNT.", "NA_Nonnormalize.")
        fData(sqaProtein$eset)[,str_replace(colnames(na_cnt_info2), "NA_IMP_CNT.", "NA_Nonnormalize.")] <- na_cnt_info2
      }
    }else if (fileType == "DIANN_Peptide"){
      na_cnt_info[na_cnt_info > 0] <- -1
      na_cnt_info[na_cnt_info == 0] <- NA
      na_cnt_info[is.na(na_cnt_info)] <- expset[is.na(na_cnt_info)]
      na_cnt_info[na_cnt_info == -1] <- NA
      if (!(userOptions$log2_i)){
        na_cnt_info <- 2^na_cnt_info
      }
      # Note: this should be activate if we want the no normalize - no imputated values
      # Note: needed for the future plots

      if (!(userOptions$SRawDataAnalysis)){
        gf <- pData(sqaProtein$eset)$globalNormFactors
        na_cnt_info2 <- sweep(na_cnt_info, 2, gf, FUN="/")

      }
      if(!(userOptions$log2_i)){
        na_cnt_info <- log2(na_cnt_info)
        na_cnt_info2 <- log2(na_cnt_info2)
      }
      colnames(na_cnt_info2) <- str_replace(colnames(na_cnt_info2), "NA_IMP_CNT.", "NA_Nonnormalize.")
      fData(sqaProtein$eset)[,ext_add_imp_col_repl] <- na_cnt_info
      fData(sqaProtein$eset)[,str_replace(colnames(na_cnt_info2), "NA_IMP_CNT.", "NA_Nonnormalize.")] <- na_cnt_info2
    }else if (fileType != "MaxQuantProteinGroup" && fileType != "SpectronautProteinGroup"){
      na_cnt_info[na_cnt_info == fData(sqaProtein$eset)[, c('nbPeptides')]] <- -1
      na_cnt_info[na_cnt_info < fData(sqaProtein$eset)[,c('nbPeptides')] & na_cnt_info >= 0] <- NA
      na_cnt_info[is.na(na_cnt_info)] <- expset[is.na(na_cnt_info)]
      na_cnt_info[na_cnt_info == -1] <- NA
      # if(!(userOptions$log2_i)){
      #   na_cnt_info <- 2^na_cnt_info
      # }

      # Note: this should be activate if we want the no normalize - no imputated values
      # Note: needed for the future plots

      # if(!(userOptions$SRawDataAnalysis)){
      #   print ("in the peptide leve")
      #   gf <- pData(sqaProtein$eset)$globalNormFactors
      #   na_cnt_info <- sweep(na_cnt_info, 2, gf , FUN="/")
      # }

      if (!(userOptions$log2_i)){
        na_cnt_info <- log2(na_cnt_info)
      }
      fData(sqaProtein$eset)[,ext_add_imp_col_repl] <- na_cnt_info

    }else if (fileType == "SpectronautProteinGroup"){
      na_cnt_info <- expset
      ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("PG.NrOfPrecursorsMeasured", names(fData(sqaProtein$eset)))]
      ext_add_imp_col_repl <- sub("\\.PG\\.NrOfPrecursorsMeasured", "\\.Intensity", ext_add_imp_col)

      if(!(userOptions$SRawDataAnalysis)){
        gf <- pData(sqaProtein$eset)$globalNormFactors
        na_cnt_info <- sweep(na_cnt_info, 2, gf , FUN="/")
      }
      fData(sqaProtein$eset)[,ext_add_imp_col_repl] <- na_cnt_info
    }else{
      na_cnt_info[na_cnt_info == 1] <- -1
      na_cnt_info[na_cnt_info == 0] <- NA
      na_cnt_info[is.na(na_cnt_info)] <- expset[is.na(na_cnt_info)]
      na_cnt_info[na_cnt_info == -1] <- NA
      # if(!(userOptions$log2_i)){
      #   na_cnt_info <- 2^na_cnt_info
      # }

      # Note: this should be activate if we want the no normalize - no imputated values
      # Note: needed for the future plots

      # if(!(userOptions$SRawDataAnalysis)){
      #   print ("in the peptide leve")
      #   gf <- pData(sqaProtein$eset)$globalNormFactors
      #   na_cnt_info <- sweep(na_cnt_info, 2, gf , FUN="/")
      # }

      # if (!(userOptions$log2_i)){
      #   na_cnt_info <- log2(na_cnt_info)
      # }
      fData(sqaProtein$eset)[,ext_add_imp_col_repl] <- na_cnt_info
    }
  }
}

# Temporary Solution for the Spectronaut files
# Solution is base the assumption that no NA values can been found in the file
if (fileType != "SpectronautProteinGroup"){
  # The sqaProtein_noImp is needed to plot the selected proteins plot before the imputated values
  sqaProtein_noImp <- sqaProtein
  ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("\\.Quantity$", names(fData(sqaProtein$eset)))]
  expr_noImp <- as.matrix(fData(sqaProtein$eset)[,ext_add_imp_col])
  colnames(expr_noImp) <- colnames(exprs(sqaProtein_noImp$eset))
  exprs(sqaProtein_noImp$eset) <- expr_noImp
}
if (length(names(fData(sqaProtein$eset))[grepl("^NA_Nonnormalize\\.", names(fData(sqaProtein$eset)))]) > 1 |
    length(names(fData(sqaProtein$eset))[grepl("^Intensity.", names(fData(sqaProtein$eset)))]) > 1){

  # The sqaProtein_noImp is needed to plot the selected proteins plot before the imputated values
  sqaProtein_noNorm <- sqaProtein
  if (fileType == "DiaNNProteinGroup" || fileType == "DIANN_Peptide"){
    # ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("^NA_Nonnormalize\\.", names(fData(sqaProtein$eset)))]
    # expr_noImp <- as.matrix(fData(sqaProtein$eset)[,ext_add_imp_col])
    # expr_noImp <- expr_noImp[,c(paste("Intensity.", pData(sqaProtein$eset)$c_Name, sep = ""))]
    ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("^NA_Nonnormalize\\.", names(fData(sqaProtein$eset)))]
    expr_noImp <- as.matrix(fData(sqaProtein$eset)[,ext_add_imp_col])
    expr_noImp <- expr_noImp[,c(paste("NA_Nonnormalize.", pData(sqaProtein$eset)$f_pos, sep = ""))]
  }else if(fileType == "MaxQuantProteinGroup"){

    ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("^Intensity.", names(fData(sqaProtein$eset)))]
    expr_noImp <- as.matrix(fData(sqaProtein$eset)[,ext_add_imp_col])

    # expr_noImp <- expr_noImp[,c(paste("Intensity ", pData(sqaProtein$eset)$file, sep = ""))]
  }
  colnames(expr_noImp) <- colnames(exprs(sqaProtein_noNorm$eset))
  # exprs(sqaProtein_noNorm$eset) <- 2^expr_noImp

  exprs(sqaProtein_noNorm$eset) <- expr_noImp
}else if (length(names(fData(sqaProtein$eset))[grepl("Intensity", names(fData(sqaProtein$eset)))]) > 1) {
  sqaProtein_noNorm <- sqaProtein
  ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("\\.Intensity", names(fData(sqaProtein$eset)))]
  expr_noImp <- as.matrix(fData(sqaProtein$eset)[,ext_add_imp_col])

  colnames(expr_noImp) <- colnames(exprs(sqaProtein_noNorm$eset))
  # exprs(sqaProtein_noNorm$eset) <- 2^expr_noImp

  exprs(sqaProtein_noNorm$eset) <- expr_noImp
}

# get the pca plot - possibilite to make it an optional
pca_flag <- TRUE
# The below code is causing problem with the output plots
if (pca_flag){
  pc <- pca(sqaProtein$eset, method = "nipals")
  df <- merge(pData(sqaProtein$eset), scores(pc), by = 0)
  if ("o_pos" %in% colnames(expDesign)){
    df$o_pos <- as.integer(df$o_pos)
    df <- df[order(-df$isControl, as.numeric(df$o_pos)),]
    rownames(df) <- df$Row.names
  }else{
    file_o <- rownames(pData(sqaProtein$eset))
    rownames(df) <- df$Row.names
    df <- df[file_o,]
  }

  pData(sqaProtein$eset) <- df
}


############################################################### EXPORTS ###############################################################
cat("INFO: PREPARING EXPORTS","\n")

#### SET WORKING DIR



if(!file.exists(userOptions$outputDir)) dir.create(userOptions$outputDir)
if(userOptions$verbose) cat("INFO: CREATED DIRECTORY",  userOptions$outputDir,"\n")
#### SET WORKING DIR

##I/O: set export file paths
userOptions$pdfFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,".pdf",sep=""))
userOptions$peptideReportFilePath <- file.path(userOptions$outputDir, paste0(userOptions$resultsFileLabel,"_PEPTIDE.",userOptions$sSheetExtension))
userOptions$proteinReportFilePath <- file.path(userOptions$outputDir, paste0(userOptions$resultsFileLabel,"_PROTEIN.",userOptions$sSheetExtension))
userOptions$paramsFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,"_SQ_PARAMS.TXT",sep=""))
userOptions$rDataFilePath <- file.path(userOptions$outputDir, paste(userOptions$resultsFileLabel,"_SQ.rData",sep=""))

############################################################### GRAPHICS ###############################################################
if (!(userOptions$noPDF)){

  # plot protein or peptide level results
  if(exists("sqaProtein")){
  	sqaDisp <- sqaProtein
  	lab <- "Protein"
  }else{
  	sqaDisp <- sqaPeptide
  	lab <- "Peptide"
  }

  # Change the column names

  ### only disp. a subset for some plots
  rowSelEset <- 1:nrow(eset) %in% sample(nrow(eset),min(c(2000,nrow(eset))) ,replace=F)
  rowSelSqaDisp <- 1:nrow(sqaDisp$eset) %in% sample(nrow(sqaDisp$eset),min(c(2000,nrow(sqaDisp$eset))) ,replace=F)
  pdf(userOptions$pdfFile)
  parDefault <- par()
  CONDITIONCOLORS <- .getConditionColors(esetNorm)
  ### EXPDESIGN PLOT
  # Note: visualize the information pData(esetNorm)
  plotExpDesign(esetNorm, version=VERSION)
  ### EXPDESIGN PLOT END
  ### IDENTIFICATION PLOTS
  if(userOptions$verbose) cat("INFO: IDENTIFICATION PLOTS \n")
  par(mfrow=c(2,2))

  #.idOverviewPlots()
  #@ NOT CRAN COMPATIBLE

  if(exists("sqaPeptide")){
    sPe <- sqaPeptide
  }else{
    sPe <- 0
  }
  if(exists("sqaProtein")){
    sPr <- sqaProtein
  }else{
    sPr <- 0
  }

  .idOverviewPlots(userOptions=userOptions
  		,esetNorm=esetNorm
  		,fileType=fileType
  		,sqaPeptide= sPe # HACK to pass check
  		,sqaProtein= sPr # HACK to pass check
  )
  if(fileType %in% c("ProgenesisFeature","ProgenesisPeptide")){
  	par(mfrow=c(3,2))
  	.idPlots(eset, selection=c(1,3), main="Feature Level", qvalueThrs=userOptions$fdrCutoff, userOptions=userOptions)
  	if(exists("sqaPeptide")) .idPlots(sqaPeptide$eset, selection=c(1,3), main="Peptide Level", qvalueThrs=userOptions$fdrCutoff)
  	if(exists("sqaProtein")) .idPlots(sqaProtein$eset, selection=c(1,3), main="Protein Level", qvalueThrs=userOptions$fdrCutoff)
  }
  par(parDefault)
  ### IDENTIFICATIONS PLOTS END
  ### QUANT. QC PLOTS
  if(userOptions$verbose) cat("INFO: QUANT QC. PLOTS \n")

  ### MASS ERROR
  par(parDefault)
  #if("pMassError" %in% names(fData(eset))){
  if(fileType %in% c("ProgenesisFeature","ProgenesisPeptide")){
  	par(mfrow=c(2,1), mar=c(4.5,6.1,4.1,6.1))
    # Note: Visualize the data fData(eset)
  	plotPrecMassErrorDistrib(eset, pMassTolWindow=userOptions$precursorMassFilter)

  	plotPrecMassErrorVsScore(eset[rowSelEset,], pMassTolWindow=userOptions$precursorMassFilter)
  	par(parDefault)
  }

  ### Violin plots ###

  # if (fileType %in% c("DiaNNProteinGroup", "SpectronautProteinGroup", "DIANN_Peptide")){
    par(mfrow=c(2,1), mar=c(4.5,6.1,4.1,6.1))

    if(userOptions$SRawDataAnalysis){
      p1 <- plot_violin(sqaProtein$eset, fType = fileType) + ggtitle("Violin plots of the Non-normalize values") + xlab("File")
    }else{
      p1 <- plot_violin(sqaProtein$eset, fType = fileType) + ggtitle("Violin plots of the normalize values") + xlab("File")
    }

    if (length(names(fData(sqaProtein$eset))[grepl("^NA_Nonnormalize\\.", names(fData(sqaProtein$eset)))]) > 1 |
        length(names(fData(sqaProtein$eset))[grepl("^Intensity", names(fData(sqaProtein$eset)))]) > 1){
      p2 <-plot_violin(sqaProtein_noNorm$eset, fType = fileType) + ggtitle("Violin plots of the Non-normalize values") + xlab("File")
      grid.arrange(p2, p1, nrow = 2)
    }else if(length(names(fData(sqaProtein$eset))[grepl("Intensity", names(fData(sqaProtein$eset)))]) > 1){
      p2 <-plot_violin(sqaProtein_noNorm$eset, fType = fileType) + ggtitle("Violin plots of the Non-normalize values") + xlab("File")
      grid.arrange(p2, p1, nrow = 2)
    }else{
      print(p1)
    }
  # }


  par(parDefault)

  ### PCA plot
  if (pca_flag){
    if(userOptions$verbose) cat("INFO: PCA plot \n")
    par(parDefault)
    plot(ggplot(pData(sqaDisp$eset), aes(PC1, PC2, shape=isControl, color = condition)) +
           geom_point(size = 3) +
           xlab("PC1") +
           ylab("PC2") +
           theme_bw())
  }

  layout(rbind(c(1), c(2,2)))
  ### missing values
  par( mar=c(6.5,5.1,2.5,3.1))
  missinValueBarplot(eset)

  ### total intensity sum
  # barplotMSSignal(eset)
  rm(barplotMSSignal)
  par( mar=c(6.5,5.1,2.5,3.1))
  cvBoxplot(sqaDisp$eset)

  par(parDefault)

  ### CORRELATION PLOTS
  ### COR OR PAIRS PLOT. IF FEWER THAN X SAMPLES

  if(ncol(sqaDisp$eset) < 8){
  	pairsAnnot(log10(exprs(sqaDisp$eset))[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered ,],textCol=as.character(CONDITIONCOLORS[pData(sqaDisp$eset)$condition,]))
  }else{
    d_c <- log10(exprs(sqaDisp$eset))[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered,]
    colnames(d_c) <- str_replace(str_replace(expDesign$file, "Quantity.", ""), ".raw", "")

  	#.correlationPlot(log10(exprs(sqaDisp$eset))[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered,], labels=as.character(unique(pData(sqaDisp$eset)$condition)), textCol=as.character(CONDITIONCOLORS[pData(sqaDisp$eset)$condition,]))
  	.correlationPlot(d_c, labels=as.character(unique(pData(sqaDisp$eset)$condition)), textCol=as.character(CONDITIONCOLORS[pData(sqaDisp$eset)$condition,]))

  	}

  ### COR OR PAIRS PLOT. IF FEWER THAN X CONDITIONS
  if(length(unique(pData(sqaDisp$eset)$condition)) < 8){
  	pairsAnnot(log10(getSignalPerCondition(sqaDisp$eset[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered,]))[,as.character(unique(pData(sqaDisp$eset)$condition)) ],textCol=as.character(CONDITIONCOLORS[as.character(unique(pData(sqaDisp$eset)$condition)),]))
  }else{
  	.correlationPlot(log10(getSignalPerCondition(sqaDisp$eset[rowSelSqaDisp & !fData(sqaDisp$eset)$isFiltered,]))[,as.character(unique(pData(sqaDisp$eset)$condition)) ],textCol=as.character(CONDITIONCOLORS[as.character(unique(pData(sqaDisp$eset)$condition)),]))
  }

  par(parDefault)

  ### TMT calibration mix
  if(exists("intAdjObj")){

  	# plot adjusted ratios vs org ratio
  	# boxplot noise fraction
  	#if(ncol(sqaDisp$ratio) > 1) par(mfrow=c(2,2))
  	boxplot(intAdjObj$noiseFraction*100, border=ifelse(intAdjObj$selectedPairs,"blue","black")
  			, ylab="Noise Fraction (%)",xlab="Calibration Mix Pair", cex.axis=1.5,cex.lab=1.5)

  	if(ncol(sqaDisp$ratio) > 1) par(mfrow=c(2,2))
  	plotAdjustedVsNonAdjustedRatio(sqaDisp$ratio,sqaDisp$unAdjustedRatio)
  	par(parDefault)

  }



  ### QUANT. QC PLOTS END

  par(parDefault)
  if(userOptions$verbose) cat("INFO: HEAT MAP \n")
  #hClustHeatMapOld(sqaDisp$eset[(1:nrow(sqaDisp$eset) %in% sample(nrow(sqaDisp$eset),min(c(nrow(sqaDisp$eset),10000)))) &  !fData(sqaDisp$eset)$isFiltered,],main= paste(lab,"Level"))

  # set smaller range if TMT
  breaks=seq(-2,2,length=20)
  if(fileType == "ScaffoldTMT" ){
    breaks=seq(-1.5,1.5,length=20)
  }
  hClustHeatMap(sqaDisp$eset[(1:nrow(sqaDisp$eset) %in% sample(nrow(sqaDisp$eset),min(c(nrow(sqaDisp$eset),10000)))) &  !fData(sqaDisp$eset)$isFiltered,]
                ,main= paste(lab,"Level")
                ,breaks = breaks)

  ### QUANT. STAT. PLOTS

  ### VAILD FEATURES VS. pValue/qValue
  if(userOptions$verbose) cat("INFO: QUANT RES. PLOTS \n")

  par(mfrow=c(1,2))
  plotNbValidDeFeaturesPerFDR(sqaDisp,
  		upRegulated=F
  		,log2RatioCufOff=log2(userOptions$ratioCutOff)
  		,pvalRange=c(0,0.15)
  		,pvalCutOff=userOptions$deFdrCutoff
  		,isLegend=T
  		,isAdjusted=T
  		,ylab=paste(lab, "Counts")
  		,main="DOWN REGULATION"
  )

  plotNbValidDeFeaturesPerFDR(sqaDisp,
  		upRegulated=T
  		,log2RatioCufOff=log2(userOptions$ratioCutOff)
  		,pvalRange=c(0,0.15)
  		,pvalCutOff=userOptions$deFdrCutoff
  		,isLegend=F
  		,isAdjusted=T
  		,ylab=paste(lab, "Counts")
  		,main="UP REGULATION"
  )

  par(parDefault)

  # plotVolcano(sqaDisp
  # 		, main=paste(lab,"Level")
  # 		, ratioThrs= userOptions$ratioCutOff
  # 		, pValueThreshold= userOptions$deFdrCutoff
  # 		, adjusted = T)

  # Note Check sqaDisp
  plotAllGGVolcanoes(sqaDisp
                     ,log2RatioThrs =userOptions$ratioCutOff %>% log2
                     ,pValueThrs= userOptions$deFdrCutoff
                     ,ylab = "log10 Adj. P-Value"
                     ,title = paste(lab,"Level")
                     ,textSize = 15
  )


  ###### Here is were need to change the plot legend for ebayes

  if(userOptions$eBayes){

  	par(mfrow=c(1,2))
  	plotNbValidDeFeaturesPerFDR(sqaDisp,
  			upRegulated=F
  			,log2RatioCufOff=log2(userOptions$ratioCutOff)
  			,pvalRange=c(0,0.15)
  			,pvalCutOff=userOptions$deFdrCutoff
  			,isLegend=T
  			,isAdjusted=F
  			,ylab=paste(lab, "Counts")
  			,main="DOWN REGULATION"
  	)

  	plotNbValidDeFeaturesPerFDR(sqaDisp,
  			upRegulated=T
  			,log2RatioCufOff=log2(userOptions$ratioCutOff)
  			,pvalRange=c(0,0.15)
  			,pvalCutOff=userOptions$deFdrCutoff
  			,isLegend=F
  			,isAdjusted=F
  			,ylab=paste(lab, "Counts")
  			,main="UP REGULATION"
  	)
  	par(parDefault)

  	# plotVolcano(sqaDisp
  	# 		, main=paste(lab,"Level")
  	# 		, ratioThrs= userOptions$ratioCutOff
  	# 		, pValueThreshold= userOptions$deFdrCutoff
  	# 		, adjusted = F)

  	plotAllGGVolcanoes(sqaDisp
  	                   ,log2RatioThrs =userOptions$ratioCutOff %>% log2
  	                   ,pValueThrs= userOptions$deFdrCutoff
  	                   ,ylab = "log10 P-Value"
  	                   ,title = paste(lab,"Level")
  	                   ,textSize = 15
  	                   ,isAdjusted = F
  	)

  	par(mfrow=c(2,2))
  	if(nrow(CONDITIONCOLORS) > 4) par(mfrow=c(3,3))
  	.allpValueHist(sqaDisp)
  	plotQValueVsPValue(sqaDisp, lim=c(0,1))
  	par(parDefault)
  }

  ### QUANT. STAT. PLOTS END

  par(parDefault)


  ### SOME ADDITIONAL QC PLOTS

  if(userOptions$addQC){

  #	if(exists("sqaPeptide")){
  #		plotXYDensity(fData(sqaPeptide$eset)$retentionTime,fData(sqaPeptide$eset)$pMassError, disp=c("")
  #			, xlab="Retention time (min)"
  #			, ylab="Precursor Mass Error (ppm)"
  #			, cex.axis=1.5
  #			, cex.lab=1.5)

  	if( all(c("retentionTime","pMassError")  %in% names(fData(eset)) )){
  		plotXYDensity(fData(eset)$retentionTime,fData(eset)$pMassError, disp=c("")
  			, xlab="Retention time (min)"
  			, ylab="Precursor Mass Error (ppm)"
  			, cex.axis=1.5
  			, cex.lab=1.5)

  		abline(h=c(userOptions$precursorMassFilter[1],0,userOptions$precursorMassFilter[2]),lty=2, lwd=2)

  		# rt vs signal
  		sel <- 1:nrow(esetNorm) %in% sample(nrow(esetNorm),min(c(4000,nrow(esetNorm))) ,replace=F) & (!(fData(esetNorm)$proteinName %in% names(CALIBMIXRATIOS)))
  		plotRTNormSummary(esetNorm[sel,])

  		par(mfrow=c(2,2))
  		plotRTNorm(getRTNormFactors(esetNorm[sel,], minFeaturesPerBin=100),esetNorm[sel,])

  		par(parDefault)

  	}

  	par(mfrow=c(2,2))
  	#all ma plots
  	for(s in colnames(exprs(esetNorm))){
  		sel <- 1:nrow(esetNorm) %in% sample(nrow(esetNorm),min(c(4000,nrow(esetNorm))) ,replace=F) & (!(fData(esetNorm)$proteinName %in% names(CALIBMIXRATIOS)))
  		maPlotSQ(esetNorm[sel,],sample=s)
  	}

  	par(parDefault)

  }


  if (userOptions$selectedProteinsList != ""){
    par(parDefault)
    plotSelectedProteins(sqaDisp, userOptions$selectedProteinsList)
    par(parDefault)
    plotSelectedProteins(sqaProtein_noImp, userOptions$selectedProteinsList, title_n = "Protein groups distribution - No imputated values")
  }

  #
  # par(parDefault)
  # plotSelectedProteins2(sqaDisp)

  rm(esetNorm)
  cat("INFO: CREATED FILE ", userOptions$pdfFile,"\n")

  graphics.off()
}
############################### GRAPHICS END

############################################################### SPREADSHEET EXPORT ###############################################################

if(exists("sqaPeptide")){

  # if("idScore" %in% names(fData(sqaPeptide$eset))){
  #   selFDataCol <- c("peptide","proteinName", "ac","geneName", "proteinDescription", "idScore","idQValue"
  #                    ,"retentionTime",	"ptm", "nbPtmsPerPeptide",	"nbRolledFeatures" )
  # }else{
  #   selFDataCol <- c("peptide","proteinName", "ac","geneName", "proteinDescription", "Global.PG.Q.Value","idQValue"
  #                    ,"retentionTime",	"ptm", "nbPtmsPerPeptide",	"nbRolledFeatures", "missCl", "Razor_genes", "Grouped",
  #                    "Protein.Ids")
  #
  # }

  # Why originally it was idScore
  if("Razor_gene" %in% names(fData(sqaPeptide$eset))){
    selFDataCol <- c("peptide","proteinName", "ac","geneName", "proteinDescription", "idScore","idQValue"
                     ,"retentionTime",	"ptm", "nbPtmsPerPeptide",	"nbRolledFeatures" )
  }else{
    selFDataCol <- c("peptide","proteinName", "ac","geneName", "proteinDescription", "Global.PG.Q.Value","idQValue"
                     ,"retentionTime",	"ptm", "nbPtmsPerPeptide",	"nbRolledFeatures", "missCl", "Razor_genes", "Grouped",
                     "Protein.Ids")

  }


	selFDataCol <-	selFDataCol[selFDataCol %in% names(fData(sqaPeptide$eset))]

	### add modif coord
	if("motifX" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"motifX","modifCoord")
	}

	### add allAccessions
	if("allAccessions" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"allAccessions")
	}

	### add ptmPeptide
	if("ptmPeptide" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"ptmPeptide")
	}

	### add ptmLocProb
	if("ptmLocProb" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"ptmLocProb")
	}

	### add ptmLocMascotConfidence
	if("ptmLocMascotConfidence" %in% names(fData(sqaPeptide$eset))){
		selFDataCol <- c(selFDataCol,"ptmLocMascotConfidence")
	}

	# add sqa object data
	cv <- sqaPeptide$cv
	names(cv) <- paste(names(cv), "Imp_cv", sep=".")
	ratio <- sqaPeptide$ratio
	if(ncol(ratio) > 0 ) names(ratio) <- paste(names(ratio), "log2ratio", sep="_")
	pValue <- sqaPeptide$pValue
	if(ncol(pValue) > 0 ) names(pValue) <- paste(names(pValue), "pValue", sep="_")
	log10_pValue <- -log10(sqaPeptide$pValue)     # Added: 26072021
	if(ncol(log10_pValue) > 0 ) names(log10_pValue) <- paste(" -log10_pValue",names(log10_pValue),sep="_") # Added: 26072021
	qValue <- sqaPeptide$qValue
	if(ncol(qValue) > 0 )  names(qValue) <- paste(names(qValue), "qValue",sep="_")
	log10_qValue <- -log10(sqaPeptide$qValue)         # Added: 26072021
	if(ncol(log10_qValue) > 0 )  names(log10_qValue) <- paste(" -log10_qValue",names(log10_qValue),sep="_")  # Added: 26072021

	if (userOptions$medianInfo){
    medianSignalDf <- log2(getSignalPerCondition(sqaPeptide$eset))
	  #names(medianSignalDf) <- paste("medianInt",names(medianSignalDf),sep="_")
    if (userOptions$log2_i){
      medianSignalDf <- 2^medianSignalDf
      names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_medianInt", sep="_")
    }else{
      names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_Log2_medianInt", sep="_")
    }
	}else{
	  medianSignalDf <- log2(getSignalPerCondition(sqaPeptide$eset, method='mean'))
	  #names(medianSignalDf) <- paste("medianInt",names(medianSignalDf),sep="_")
	  if (userOptions$log2_i){
	    medianSignalDf <- 2^medianSignalDf
	    names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_meanInt", sep="_")
	  }else{
	    names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_Log2_meanInt", sep="_")
	  }
	}
	expset <- exprs(sqaPeptide$eset)
	expset <- log2(expset)

	if ("file" %in% colnames(expDesign)){
	  colnames(expset) <- expDesign$file
	}else{
	  colnames(expset) <- rownames(expDesign)
	}

	if (userOptions$log2_i){
	  expset <- 2^expset
	  colnames(expset) <- sub("Quantity\\.", "", colnames(expset))
	  colnames(expset) <- paste(colnames(expset), "Imp.Quantity", sep=".")
	}else{
	  colnames(expset) <- sub("Quantity\\.", "", colnames(expset))
	  colnames(expset) <- paste(colnames(expset), "Imp_Log2.Quantity", sep=".")
	}

	expdDesign_file_name_cnv <- gsub("-", '.', sub('Quantity\\.', '', expDesign$file))
	ext_add_col <- c()
	# Alternative solution
	if(length(names(fData(sqaPeptide$eset))[grepl("^NrOfPrecursors\\.",names(fData(sqaPeptide$eset)))]) > 0){

	  for ( name_i in expdDesign_file_name_cnv){
	    ext_add_col <- c(ext_add_col, paste('NrOfPrecursors', name_i, sep = "."))
	  }
	}

	selFDataCol <- c(selFDataCol, ext_add_col)

  if ("rt" %in% statMethod || "quantile" %in% statMethod){
    ext_add_imp_col <- names(fData(sqaPeptide$eset))[grepl("^NA_IMP_CNT",names(fData(sqaPeptide$eset)))]
    ext_add_imp_col_repl <- sub('NA_IMP_CNT\\.', '', ext_add_imp_col)
    if ("file" %in% colnames(expDesign)){
      ext_add_imp_col_repl <- expDesign[ext_add_imp_col_repl, 'file']
      ext_add_imp_col_repl <- sub("Quantity\\.", "", ext_add_imp_col_repl)
    }
    if (userOptions$log2_i && userOptions$SRawDataAnalysis){
      if (fileType != "MaxQuantProteinGroup"){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Quantity", sep=".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm.Quantity", sep=".")
      }
    }else if (userOptions$log2_i){
      ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm.Quantity", sep=".")
    }else if (userOptions$SRawDataAnalysis){
      if (fileType != "MaxQuantProteinGroup" && userOptions$SRawDataAnalysis){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Log2.Quantity", sep=".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm_Log2.Quantity", sep=".")
      }
    }else{
      ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm_Log2.Quantity", sep=".")
    }
    selFDataCol <- c(selFDataCol, ext_add_imp_col)
    ext_add_imp_col2 <- names(fData(sqaPeptide$eset))[grepl("^NA_IMP_CNT",names(fData(sqaPeptide$eset)))]

    # I replace the 0 and 1 of the original store NA_IMP_CNT to the the actually values with NA when there was an imputated value
    # Note: In the futere maybe this values should come after the original setting.
    # I should also change the name but for know I will leave it like this
    na_cnt_info <- fData(sqaPeptide$eset)[,ext_add_imp_col2]
    na_cnt_info[na_cnt_info >= 1 ] <- -1
    na_cnt_info[na_cnt_info == 0 ] <- NA
    na_cnt_info[is.na(na_cnt_info)] <-expset[is.na(na_cnt_info)]
    na_cnt_info[na_cnt_info == -1 ] <- NA
  }else{
    ext_add_imp_col <- names(fData(sqaPeptide$eset))[grepl("^NA_IMP_CNT",names(fData(sqaPeptide$eset)))]
    ext_add_imp_col_repl <- sub('NA_IMP_CNT\\.', '', ext_add_imp_col)
    if ("file" %in% colnames(expDesign)){
      ext_add_imp_col_repl <- expDesign[ext_add_imp_col_repl, 'file']
      ext_add_imp_col_repl <- sub("Quantity\\.", "", ext_add_imp_col_repl)
    }
    if (userOptions$log2_i){
      if (fileType != "MaxQuantProteinGroup" && userOptions$SRawDataAnalysis){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Quantity", sep=".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm.Quantity", sep=".")
      }
    }else{
      if (fileType != "MaxQuantProteinGroup" && userOptions$SRawDataAnalysis){
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Log2.Quantity", sep=".")
      }else{
        ext_add_imp_col_repl <- paste(ext_add_imp_col_repl, "Norm_Log2.Quantity", sep=".")
      }
    }
    selFDataCol <- c(selFDataCol, ext_add_imp_col)
    ext_add_imp_col2 <- names(fData(sqaPeptide$eset))[grepl("^NA_IMP_CNT",names(fData(sqaPeptide$eset)))]

    # I replace the 0 and 1 of the original store NA_IMP_CNT to the the actually values with NA when there was an imputated value
    # Note: In the futere maybe this values should come after the original setting.
    # I should also change the name but for know I will leave it like this
    na_cnt_info <- fData(sqaPeptide$eset)[,ext_add_imp_col2]
    na_cnt_info[na_cnt_info >= 1 ] <- -1
    na_cnt_info[na_cnt_info == 0 ] <- NA
    na_cnt_info[is.na(na_cnt_info)] <-expset[is.na(na_cnt_info)]
    na_cnt_info[na_cnt_info == -1 ] <- NA
    if (!(userOptions$log2_i)){
      na_cnt_info <- 2^na_cnt_info
    }

    # gf <- pData(sqaPeptide$eset)$globalNormFactors
    # na_cnt_info <- sweep(na_cnt_info, 2, gf , FUN="/")
    if (!(userOptions$log2_i)){
      na_cnt_info <- log2(na_cnt_info)
    }
  }


	fData(sqaPeptide$eset)[,ext_add_imp_col] <- na_cnt_info

	out <- cbind(
			fData(sqaPeptide$eset)[,selFDataCol]
			, expset
			, medianSignalDf
			, cv
			, ratio
			# , fracNAFeatures = getNAFraction(sqaPeptide$eset,method=c("cond","count"))
			, pValue
			, log10_pValue
			, qValue
			, log10_qValue
			# , FTestPValue = sqaPeptide$FPValue
			# , FTestQValue = sqaPeptide$FQValue
			)[!fData(sqaPeptide$eset)$isFiltered,]

	colnames(out)[colnames(out) %in% ext_add_imp_col] <- ext_add_imp_col_repl

	### add unadjusted ratios if TMT ratio correction
	if(userOptions$TAdjustRatios){
		unadjPeptideRatios <- sqaPeptide$unAdjustedRatio[!fData(sqaPeptide$eset)$isFiltered,]
		names(unadjPeptideRatios) <- paste(names(sqaPeptide$ratio), "log2_unadjRatio", sep=".")
		out <- cbind(out,unadjPeptideRatios)
	}

	### paired expDesign ratio export
	if("subject" %in% names(pData(sqaPeptide$eset))){
		allRatios <- getRatios(sqaPeptide$eset,method="paired")[!fData(sqaPeptide$eset)$isFiltered,]
		names(allRatios) <- paste(names(allRatios), "log2_pairedRatio", sep=".")
		out <- cbind(out,allRatios)
	}
	if (fileType == "DIANN_Peptide"){
	  out <- out[,!(colnames(out) %in% c("ac", "idQValue", "missCl", "allAccessions"))]
	  new_order <- c("peptide", "proteinName", "Protein.Ids", "geneName", "Razor_genes", "proteinDescription", "Global.PG.Q.Value", "retentionTime", "ptm", "nbPtmsPerPeptide")
	  new_order <- c(new_order, names(out)[!(names(out) %in% new_order)])
	  out <- out[, new_order]
	}

	write.table(out
			, file=userOptions$peptideReportFilePath
			, sep=userOptions$sSheetExportDelimiter
			, row.names=F
	)
	if (fileType == "DIANN_Peptide"){
	  peptide_df <- out
	}

	cat("INFO: CREATED FILE ", userOptions$peptideReportFilePath,"\n")
}

if(exists("sqaProtein")){

  if ("idScore" %in% colnames(fData(sqaProtein$eset))){
    selFDataCol <- c("proteinName","ac","geneName","proteinDescription","idScore")
  }else if ("Global.PG.Q.Value" %in% colnames(fData(sqaProtein$eset))){
    if("RazorG" %in% colnames(fData(sqaProtein$eset))){
      selFDataCol <- c("proteinName","ac", "RazorP", "geneName", "RazorG", "proteinDescription","Global.PG.Q.Value")
    }else{
      selFDataCol <- c("proteinName","ac","geneName","proteinDescription","Global.PG.Q.Value")
    }
  }else if("Protein.Q.Value" %in% colnames(fData(sqaProtein$eset))){
    selFDataCol <- c("proteinName","ac","geneName","proteinDescription","Protein.Q.Value")
  }
	#selFDataCol <- c("proteinName","ac","geneName","proteinDescription","idScore","idQValue")
	if("nbPeptides" %in% names(fData(sqaProtein$eset))){
	  selFDataCol <- c(selFDataCol, 'nbPeptides')
	}else if(fileType == "DiaNNProteinGroup"){
	  col_info <- names(fData(sqaProtein$eset))[grepl("NrOfPrecursorsMeasured\\.",names(fData(sqaProtein$eset)))]
	  fData(sqaProtein$eset)[,"nbPeptides"] <- apply(fData(sqaProtein$eset)[,col_info], 1, max)
	  selFDataCol <- c(selFDataCol, 'nbPeptides')
	}

  if(fileType == "MaxQuantProteinGroup"){
    col_info <- names(fData(sqaProtein$eset))[grepl("razor_unique_peptides", names(fData(sqaProtein$eset)))]
    selFDataCol <- c(selFDataCol, col_info)
    col_info <- names(fData(sqaProtein$eset))[grepl("^Intensity.", names(fData(sqaProtein$eset)))]
    if(!(userOptions$log2_i)){
      fData(sqaProtein$eset)[,col_info] <- log2(fData(sqaProtein$eset)[, col_info])
    }
    selFDataCol <- c(selFDataCol, col_info)
  }

	selFDataCol <- selFDataCol[selFDataCol %in%  names(fData(sqaProtein$eset))]

	### add allAccessions
	if("allAccessions" %in% names(fData(sqaProtein$eset))){
		selFDataCol <- c(selFDataCol,"allAccessions")
	}

	# add sqa object data

	cv <- sqaProtein$cv
	names(cv) <- paste(names(cv), "Imp_cv",sep=".")
	ratio <- sqaProtein$ratio
	if(ncol(ratio) > 0 ) names(ratio) <- paste(names(ratio), "log2ratio",sep=".")
	pValue <- sqaProtein$pValue
	if(ncol(pValue) > 0 ) names(pValue) <- paste(names(pValue), "pValue",sep=".")
	log10_pValue <- -log10(sqaProtein$pValue)      # Added: 26072021
	if(ncol(log10_pValue) > 0 ) names(log10_pValue) <- paste(" -log10_pValue",names(log10_pValue),sep=".")  # Added: 26072021
	qValue <- sqaProtein$qValue
	if(ncol(qValue) > 0 ) names(qValue) <- paste(names(qValue), "qValue",sep=".")
	log10_qValue <- -log10(sqaProtein$qValue)     # Added: 26072021
	if(ncol(log10_qValue) > 0 ) names(log10_qValue) <- paste(" -log10_qValue",names(log10_qValue),sep=".")  # Added: 26072021

	if (userOptions$medianInfo){

	  medianSignalDf <- log2(getSignalPerCondition(sqaProtein$eset))
	  if (userOptions$log2_i){
	    medianSignalDf <- 2^medianSignalDf
	    names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_medianInt",sep=".")
	  }else{
	    names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_Log2_medianInt",sep=".")
	  }

	}else{
	  medianSignalDf <- log2(getSignalPerCondition(sqaProtein$eset, method="mean"))
    if (userOptions$log2_i){
      medianSignalDf <- 2^medianSignalDf
      names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_meanInt",sep=".")
    }else{
      names(medianSignalDf) <- paste(names(medianSignalDf), "Imp_Log2_meanInt",sep=".")
    }
	}
	expset <- exprs(sqaProtein$eset)

	expset <- log2(expset)

	if ("file" %in% colnames(expDesign)){
	  colnames(expset) <- expDesign$file
	}else{
	  colnames(expset) <- rownames(expDesign)
	}
	if (userOptions$log2_i){
	  expset <- 2^expset
	  colnames(expset) <- sub("Quantity\\.", "", colnames(expset))
	  colnames(expset) <- paste(colnames(expset), "Imp.Quantity", sep=".")
	}else{
	  colnames(expset) <- sub("Quantity\\.", "", colnames(expset))
	  colnames(expset) <- paste(colnames(expset), "Imp_Log2.Quantity", sep=".")
	}

	if (fileType == "SpectronautProteinGroup"){
	  ext_add_col <- names(fData(sqaProtein$eset))[grepl("NrOfPrecursorsMeasured$",names(fData(sqaProtein$eset)))]
	  selFDataCol <- c(selFDataCol, ext_add_col)
	  ext_add_col <- names(fData(sqaProtein$eset))[grepl("\\.Intensity",names(fData(sqaProtein$eset)))]
	  if(!(userOptions$log2_i)){
	    fData(sqaProtein$eset)[, ext_add_col] <- log2(fData(sqaProtein$eset)[, ext_add_col])
	  }
	  selFDataCol <- c(selFDataCol, ext_add_col)
	}else if (fileType == "DiaNNProteinGroup"){
	  # ext_add_col <- names(fData(sqaProtein$eset))[grepl("^NrOfPrecursorsMeasured\\.",names(fData(sqaProtein$eset)))]
	  # ext_add_col_cnv <- sub('NrOfPrecursorsMeasured\\.', '', ext_add_col)
	  expdDesign_file_name_cnv <- gsub("-", '.', sub('Quantity\\.', '', expDesign$file))
	  # One option but the problem it doesn't give the correct order of the data according to the quantity values when the C is not the first samples
	  # for (name_i in ext_add_col_cnv){
	  #   print (name_i)
	  #   if (name_i %in% expdDesign_file_name_cnv){
	  #     next
	  #   }else{
	  #     ch_str <- paste('NrOfPrecursorsMeasured.', name_i, sep = "")
	  #     ext_add_col <- ext_add_col[ext_add_col != ch_str]
	  #   }
	  # }


	  ext_add_col <- c()
	  # Alternative solution
	  if(length(names(fData(sqaProtein$eset))[grepl("^NrOfPrecursorsMeasured\\.",names(fData(sqaProtein$eset)))]) > 0){

	    for ( name_i in expdDesign_file_name_cnv){
	      ext_add_col <- c(ext_add_col, paste('NrOfPrecursorsMeasured', name_i, sep = "."))
	    }
	  }else{
	    for ( name_i in expdDesign_file_name_cnv){
	      ext_add_col <- c(ext_add_col, paste(name_i, 'NrOfFullyPeptides', sep = "."))
	    }
	  }



	  selFDataCol <- c(selFDataCol, ext_add_col)

	  if(length(names(fData(sqaProtein$eset))[grepl("^NrOfRazorPrecursors\\.",names(fData(sqaProtein$eset)))]) > 0){
	    ex_col_v2 <- c()
	    for ( name_i in expdDesign_file_name_cnv){
	      ex_col_v2 <- c(ex_col_v2, paste("NrOfRazorPrecursors", name_i, sep = "."))
	    }
	    selFDataCol <- c(selFDataCol, ex_col_v2)
	  }

	}
  if(length(names(fData(sqaProtein$eset))[grepl("\\.Quantity$",names(fData(sqaProtein$eset)))]) > 0){
    ext_add_imp_col <- names(fData(sqaProtein$eset))[grepl("\\.Quantity$",names(fData(sqaProtein$eset)))]
    selFDataCol <- c(selFDataCol, ext_add_imp_col)
    if(!(userOptions$log2_i)){
      fData(sqaProtein$eset)[,ext_add_imp_col] <- log2(fData(sqaProtein$eset)[, ext_add_imp_col])
    }

  }




  if (userOptions$extra_info){
    out <- cbind(
      fData(sqaProtein$eset)[,selFDataCol]
      , expset
      , medianSignalDf
      , cv
      , ratio
      , fracNAFeatures = getNAFraction(sqaProtein$eset,method=c("cond","count"))
      , pValue
      , log10_pValue
      , qValue
      , log10_qValue
      , FTestPValue = sqaProtein$FPValue
      , FTestQValue = sqaProtein$FQValue
    )[!fData(sqaProtein$eset)$isFiltered,]


  }else{
    out <- cbind(
      fData(sqaProtein$eset)[,selFDataCol]
      , expset
      , medianSignalDf
      , cv
      , ratio
      #, fracNAFeatures = getNAFraction(sqaProtein$eset,method=c("cond","count"))
      , pValue
      , log10_pValue
      , qValue
      , log10_qValue
      #, FTestPValue = sqaProtein$FPValue
      #, FTestQValue = sqaProtein$FQValue
    )[!fData(sqaProtein$eset)$isFiltered,]
  }
	# if (!(exists("sqaPeptide"))){
	  colnames(out)[colnames(out) %in% ext_add_imp_col] <- ext_add_imp_col_repl
	# }


	# add median top3
	if(exists("esetTop3")){

		# medians
#		tmpOut <- getSignalPerCondition(esetTop3)
#		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
#		names(tmpOut) <- paste("medianInt_top3",names(tmpOut),sep="_")

		tmpOut <- exprs(esetTop3)
		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
		colnames(tmpOut) <- paste("top3",colnames(tmpOut),sep="_")

		out <- cbind(out,tmpOut)
	}

	# add iBAQ
	if(exists("esetIBAQ")){

		# medians
#		tmpOut <- getSignalPerCondition(esetIBAQ)
#		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
#		names(tmpOut) <- paste("medianInt_top3",names(tmpOut),sep="_")

		tmpOut <- exprs(esetIBAQ)
		if ("file" %in% colnames(expDesign)){
		  colnames(tmpOut) <- sub("Quantity\\.", "", expDesign$file)
		}else{
		  colnames(tmpOut) <- rownames(expDesign)
		}

		tmpOut <- tmpOut[match(rownames(out),rownames(tmpOut)), ]
		if (!(userOptions$log2_i)){
		  tmpOut <- log2(tmpOut)
		  colnames(tmpOut) <- paste(colnames(tmpOut), "iBAQ.Log2", sep=".")
		}else{
		  colnames(tmpOut) <- paste(colnames(tmpOut), "iBAQ", sep=".")
		}


		out <- cbind(out,tmpOut)

	}

	### add unadjusted ratios if TMT ratio correction
	if(userOptions$TAdjustRatios){
		unadjProteinRatios <- sqaProtein$unAdjustedRatio[!fData(sqaProtein$eset)$isFiltered,]
		names(unadjProteinRatios) <- paste(names(sqaProtein$ratio), "log2_unadjRatio", sep= ".")
		out <- cbind(out,unadjProteinRatios)
	}


	### paired expDesign ratio export
	if("subject" %in% names(pData(sqaProtein$eset))){
		allRatios <- getRatios(sqaProtein$eset,method="paired")[!fData(sqaProtein$eset)$isFiltered,]
		names(allRatios) <- paste(names(allRatios), "log2_pairedRatio",sep=".")
		out <- cbind(out,allRatios)
	}

	if (fileType == "DIANN_Peptide"){
    out <- proteinTable_reformat(userOptions$inputFile, out, peptide_df)
	}


	write.table(out
			, file=userOptions$proteinReportFilePath
			, sep=userOptions$sSheetExportDelimiter
			, row.names=F
	)

	cat("INFO: CREATED FILE ", userOptions$proteinReportFilePath,"\n")
}
rm(expDesign)
### SPREADSHEET EXPORT END

############################################################### PARAMS EXPORT ###############################################################
write.table(data.frame(
				param=row.names(data.frame(unlist(userOptions[names(userOptions)])))
				,value=as.vector((unlist(userOptions[names(userOptions)])))
		)
		,file=userOptions$paramsFilePath
		,sep="\t"
		,row.names=F
		,quote=F
)

cat("INFO: CREATED FILE ", userOptions$paramsFilePath,"\n")

### EXPORT PARAMS

############################################################### RDATA EXPORT ###############################################################

if(userOptions$isSaveRObject){
	save.image(file=userOptions$rDataFilePath)
	cat("INFO: CREATED FILE ", userOptions$rDataFilePath,"\n")
}
### EXPORT RDATA END
