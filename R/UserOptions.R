# TODO: Add comment
#
# Author: erikahrne, georgiaAngelidou
###############################################################################



### CMD OPTIONS
#' Command Line Option List
#' @export
option_list <- list(


### I/O
		make_option(c("-i", "--inputFile"), type="character", default="",
				help="I/O:  Input file: Progenesis (Feature,Protein or Peptide) .csv,
			or Scaffold Q+ (Raw Export, for TMT quant) .xls (REQUIRED)",
		),
		make_option(c("-o", "--outputDir"), type="character", default=NA,
				help=paste("I/O:  Results Output Directory ", red$bold("[default FOLDER OF INPUTFILE]"), sep = ""),
		),

		make_option(c("-l", "--resultsFileLabel"), type="character", default="SQ_Results",
				help=paste("I/O: results file directory ", red$bold("[default %default]"), sep = ""),
		),

		make_option(c("-f", "--fastaFile"), type="character", default="",
				help="I/O:  Protein DB .fasta file " %+% red$bold("[default ./]"),
		),

		make_option(c("-p", "--scaffoldPTMSpectrumReportFile"), type="character", default="",
				help="I/O:  Scaffold PTM Spectrum Report File " %+% red$bold("[default ./]"),
		),

		make_option(c("-d","--spreadsheetExportDelimiter"), type="integer", default=1,
		            help="I/O: Spreadsheet Export Delimiter 1) <tab> 2) <,> " %+% red$bold("[default %default]"),
    ),
		make_option(c("-s", "--selectedProteinsList"), type="character", default="",
		            help="I/O: List of Selected proteins to track their abundance through the different condition"),

		make_option(c("--noPDF"), action="store_true", default= FALSE,
		            help="O: When activated the PDF output is not created"),

### I/O END

# FILTER (--F)
		make_option(c("--FProteinAccessionSelection"), type="character", default=".",
				help="FILTER: --FP Filter features by Accession Regular Expression " %+% red$bold("[default %default]") %+% " (all features kept)",
				metavar="Protein Accession Reg. expr."),

		#### peptide analysis specfic
		make_option(c("--FModificationSelection"), type="character", default="",
				help="FILTER (LFQ PEP ONLY): --FM Only keep Peptides with modifications matching Regular Expression " %+% red$bold("[default %default]") %+% "
				(all features kept).",
				metavar="modification name Reg. expr."),

		make_option(c("--FFdrCutoff"), type="double", default=0.01,
				help="FILTER (LFQ ONLY): --FF Identification level False Discovery Rate Cutoff.  [0-1] " %+% red$bold("[default %default]"),
				metavar="Peptide/Protein FDR cutoff"),

#		make_option(c("--FCoefficientOfVarianceMax"), type="double", default=Inf,
#				help="FILTER: --FC Do not include features with C.V. above this threshold in statistical
#				test for differential expression [default %default]",
#				metavar="Coefficent of Variance cutoff"),

		#### peptide analysis specfic
		make_option(c("--FDeltaMassTolerancePrecursor"), type="character", default="AUTO SET",
				help="FILTER (LFQ PEP ONLY): --FD Precursor mass Error Range filter (ppm) " %+% red$bold("[default %default]") %+% ".
				Peptide imports ONLY",
				metavar="Mass Range [x,y]"),

		#### protein analysis specfic
		make_option(c("--FNumberOfPeptidesPerProteinMin"), type="integer", default=1,
				help="FILTER: --FN Only include those proteins with at least x identified peptides " %+% red$bold("[default %default]") %+% "
				Protein analysis ONLY.",
				metavar="Number of peptides"),

		#### peptide analysis specfic
		make_option(c("--FSitesPerPeptide"), type="integer", default=99999,
				help="FILTER: --FS Max Nb. Modifications Per Peptide " %+% red$bold("[default Inf]") %+% "
						Peptide analysis ONLY.",
				metavar="Max Number of PTM sites Per Petptide"),

		#### peptide analysis specfic
		make_option(c("--FLengthPeptide"), type="integer", default=1,
				help="FILTER: --FL Min Peptide Length (Nb. AA's) " %+% red$bold("[default 1] ") %+% "
						Peptide analysis ONLY.",
				metavar="Min Peptide Length (>=)"),

		####
		make_option(c("--FExclusivePeptides"), action="store_true", default=FALSE,
				help="FILTER: --FE Discard all peptides mapping to multiple protein entries " %+% red$bold("[default %default]") %+% "
			Note that by default all peptides are used for quantification and assigned to proteins using
			a Occam's Razor based algorithm.
				"),

		make_option(c("--FRatioCutOff"), type="double", default=1,
				help="FILTER: --FR Intensity ratio cut-off. " %+% red$bold("[default %default]"),
				metavar="Intensity ratio cutoff"),

    make_option(c("--extraOut"), action="store_true", default=FALSE,
            help="Provides extra information in the output File"),


# FILTER (--F) END

# TMT (--T)

		# correct tmt ratios
		make_option(c("--TAdjustRatios"), action="store_true", default=FALSE,
				help="TMT: --TA Adjust TMT ratios using calibration mix proteins " %+% red$bold("[default %default]")),

# TMT (--T) END

# STATISTICS (--S)


  # Use median instead of median
  make_option(c("--Median"), action="store_true", default=FALSE,
            help="The return values will correspond to the median values"
  ),

  make_option(c("--log2toI"), action="store_true", default=FALSE,
            help="transforms log2 values to their linear scale"
  ),


	make_option(c("--SAnchorProtein"), type="character", default=".",
			help="STATISTICS: --SA Normalize Intensities by selected protein(s) Regular Expression " %+% "
			" %+% red$bold("[default %default]") %+% " (use all proteins).",
			metavar="Protein Accession Reg. expr."),

  make_option(c("--SMissingValuesImutationMethod"), type="character", default="nDist",
            help="STATISTICS: --SM 'ppca', 'knn','gMin','lMin','gMean,'lMean', 'nDist', " %+% "
              " %+% red$bold("[default %default]") %+% " (use all proteins).",
            metavar=" ppca: probabilistic pca (+ gMin, if not enough data)
                            knn: k-nearest neighbour (+ gMin, if not enough data)
                            gMin: global minimum
                            lMin: local minimum
                            gMean: global mean
                            lMean: local mean
                            nDist: according to the normal Distribution
            "),


  make_option(c("--SNonPairWiseStatTest"), action="store_true", default=FALSE,
            help="STATISTICS: --SN non pairwise eBayes moderated t-statistic p-values.
              I.e. variance is pooled, per protein/peptide, across all runs of the study " %+% red$bold("[default %default]")),

  make_option(c("--SPvalueInclude"), action="store_true", default=FALSE,
            help="STATISTICS: --SP output eBayes moderated t-statistic p-values " %+% red$bold("[default %default]")),

	make_option(c("--SRawDataAnalysis"), action="store_true", default=FALSE,
			help="STATISTICS: --SR No data normalization " %+% red$bold("[default %default]")),

# STATISTICS (--S) END

# EXPERIMENTAL DESIGN (--E)

	make_option(c("--EXperimentalDesign"), type="character", default=NA,
			help='EXPERIMENTAL DESIGN: --EX "," seperated samples, ":" separated conditions
					Example: 1,2,3:4,5,6
					   condition1 (REF) : channel 1,2,3
					   condition2: channel 4,5,6
					Note: for 10-plex default is "1,4,7,10:2,5,8:3,6,9"
					' %+% red$bold('[default %default]')),

	make_option(c("--EProteinQuantOff"), action="store_false", default=TRUE,
			help='EXPERIMENTAL DESIGN: --EP Disable Protein Level Quantification ' %+% red$bold('[default %default]')),

	make_option(c("--ECorrelatedSamples "), action="store_true", default=FALSE,
			help='EXPERIMENTAL DESIGN: --EC Apply "paired" statistical tests ' %+% red$bold('[default %default]')),

# EXPERIMENTAL DESIGN (--E) END

# PDF-REPORT (--P)

	make_option(c("--PQvalueCutOff"), type="double", default=0.01,
			help="PDF-REPORT: --PQ Qvalue cut-off used for graphics.
			High-lighting features with a qval < specified value. [0-1] " %+% red$bold("[default %default]"),
			metavar="Differential expression qvalue cutOff"),

# ADDITIONAL-REPORTS (--A)
	make_option(c("--ARDataFile"), action="store_true", default=FALSE,
		help="ADDITIONAL-REPORTS: --AR Save R objects in 'label'.RData file " %+% red$bold("[default %default]")),

	make_option(c("--AIbaq"), action="store_true", default=FALSE,
			help="ADDITIONAL-REPORTS : --AI add iBAQ values to results spreadsheet. " %+% red$bold("[default %default]")),

	make_option(c("--ATop3"), action="store_true", default=FALSE,
			help="ADDITIONAL-REPORTS : --AT add Top3 values to results spreadsheet. " %+% red$bold("[default %default]")),

	make_option(c("--AQC"), action="store_true", default=FALSE,
			help="ADDITIONAL-REPORTS : --AQ adds additional QC plots to .pdf report " %+% red$bold("[default %default]")),



# ADDITIONAL-REPORTS (--A) END

# TEST (peptide analysis specific)
	make_option(c("-t", "--test"), action="store_true", default=FALSE,
			help="TEST: test option, include first 2000 entries only " %+% red$bold("[default %default]") %+% "
			Peptide analysis ONLY."),
# TEST END
	make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
			help="Print extra output " %+% red$bold("[default %default]"))
	)

#' Read User Specified Command Line Options
#' @param version Safequant version number
#' @return user options list
#' @import  optparse
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getUserOptions <- function(version=version){
  # alert <- combine_styles("bold", "red4", "bgCyan")
	epilogue <- "Examples:
	Progenesis LFQ Protein Quant:
	>Rscript safeQuant.R -i /path/to/peptide_measurment.csv

	Progenesis LFQ Protein Quant (QE):
	>Rscript safeQuant.R -i /path/to/peptide_measurment.csv --FL 7

	Progenesis LFQ Phospho Quant:
	>Rscript safeQuant.R -i /path/to/peptide_measurment.csv -f /path/to/proteins.fasta --FM phospho --FS 3 --EP

	MaxQuant Protein Quant with Experinental Design File:
	>Rscript safeQuant.R -i /path/to/proteinGroups.txt --EX experimentalDesignTemplate.txt

	"

	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	cmdOpt <- parse_args(OptionParser( prog=paste("SafeQuant",version), option_list=option_list, epilogue=epilogue))

	### CMD OPTIONS END

	### SET USER OPTIONS
	userOptions <- list()

### VERBOSE
	#VERBOSE: verbose
	userOptions$verbose <- cmdOpt$verbose
### VERBOSE	END

# I/O
	#I/O: progenesisFilePath
	userOptions$inputFile <- cmdOpt$inputFile
	if( userOptions$inputFile == "" | !file.exists(userOptions$inputFile)){
		cat("ERROR. Please specify input file.",userOptions$inputFile, "Not found!","\n")
		q(status=-1)
	}

	#I/O: resultsFileLabel
	userOptions$resultsFileLabel <- cmdOpt$resultsFileLabel

	#I/O: outputDir
	userOptions$outputDir <- cmdOpt$outputDir
	if(is.na(userOptions$outputDir)){ # see default
		userOptions$outputDir <- dirname(userOptions$inputFile)
	}
	if(!file.exists(userOptions$outputDir) & userOptions$outputDir != "" ){
		cat("ERROR. No such directory",userOptions$outputDir,"\n")
		q(status=-1)
	}else{
		userOptions$outputDir <- file.path(userOptions$outputDir, userOptions$resultsFileLabel)
	}

	#I/O: proteinFastaFile
	userOptions$proteinFastaFile <- NA
	if(nchar(cmdOpt$fastaFile) > 0 ){
		### check if file exists
		if(file.exists(cmdOpt$fastaFile)){
			userOptions$proteinFastaFile <- cmdOpt$fastaFile
		}else{
			cat("ERROR. File does not exist",cmdOpt$fastaFile,"\n")
			q(status=-1)
		}
	}

	#I/O: scaffoldPTMSpectrumReportFile
	userOptions$scaffoldPTMSpectrumReportFile <- NA
	if(nchar(cmdOpt$scaffoldPTMSpectrumReportFile) > 0 ){
		### check if file exists
		if(file.exists(cmdOpt$scaffoldPTMSpectrumReportFile)){
			userOptions$scaffoldPTMSpectrumReportFile <- cmdOpt$scaffoldPTMSpectrumReportFile
		}else{
			cat("ERROR. File does not exist",cmdOpt$scaffoldPTMSpectrumReportFile,"\n")
			q(status=-1)
		}
	}

	#I/O: spreadsheetExportDelimiter
	if(cmdOpt$spreadsheetExportDelimiter == 1){
	  userOptions$sSheetExtension = "tsv"
	  userOptions$sSheetExportDelimiter = "\t"
	}else{
	  userOptions$sSheetExtension = "csv"
	  userOptions$sSheetExportDelimiter = ","
	}


	#I/O: Selected Proteins List
	userOptions$selectedProteinsList <- cmdOpt$selectedProteinsList
	if(userOptions$selectedProteinsList != "" && !file.exists(userOptions$selectedProteinsList)){
	  cat("ERROR. File does not exist:",userOptions$selectedProteinsListe,"\n")
	  q(status=-1)
	}


# I/O END

# FILTER (--F)

	#FILTER: selectedProteinName
	userOptions$selectedProteinName <- cmdOpt$FProteinAccessionSelection

	#FILTER: selectedModifName
	userOptions$selectedModifName <- cmdOpt$FModificationSelection

	#FILTER: fdrCutoff
	userOptions$fdrCutoff <- cmdOpt$FFdrCutoff
	if(is.na(userOptions$fdrCutoff) | userOptions$fdrCutoff <= 0 | userOptions$fdrCutoff > 1 ){
		cat("ERROR. fdrCutOff must be in the range [0-1]. You specified",userOptions$fdrCutoff,"\n")
		q(status=-1)
	}

	#FILTER: precursorMassFilter
	if(cmdOpt$FDeltaMassTolerancePrecursor == "AUTO SET" ){
		userOptions$precursorMassFilter <- NA
	}else{
		### set by user -> add lower and upper mass bound to vector
		userOptions$precursorMassFilter <- gsub("(\\[)","",cmdOpt$FDeltaMassTolerancePrecursor)
		userOptions$precursorMassFilter <- gsub("(\\])","",userOptions$precursorMassFilter)
		userOptions$precursorMassFilter <- sort(as.numeric(unlist(strsplit(userOptions$precursorMassFilter,","))))

		### check input format precursorMassFilter
		if((length(userOptions$precursorMassFilter) != 2) | sum(is.na(userOptions$precursorMassFilter)) > 0 ){
			cat("ERROR. Invalid FDeltaMassTolerancePrecursor", userOptions$minNbPeptidesPerProt, "\n")
			q(status=-1)
		}
	}

	#FILTER: cvCutOff
#	userOptions$cvCutOff <- cmdOpt$FCoefficientOfVarianceMax
#	if(is.na(userOptions$cvCutOff) | (userOptions$cvCutOff < 0)){
#		print(paste("ERROR. cvCutOff must be > 0. You specified ", userOptions$cvCutOff))
#		q(status=-1)
#	}

	#FILTER: minNbPeptidesPerProt
	userOptions$minNbPeptidesPerProt <- cmdOpt$FNumberOfPeptidesPerProteinMin
	if(is.na(userOptions$minNbPeptidesPerProt) | userOptions$minNbPeptidesPerProt < 0 ){
		print(paste("ERROR. FNumberOfPeptidesPerProteinMin must be >= 0. You specified ", userOptions$minNbPeptidesPerProt))
		q(status=-1)
	}

	#FILTER: maxNbPTMsPerPeptide
	userOptions$maxNbPtmsPerPeptide <- cmdOpt$FSitesPerPeptide
	if(is.na(userOptions$maxNbPtmsPerPeptide) | userOptions$maxNbPtmsPerPeptide < 0 ){
		print(paste("ERROR. FSitesPerPeptide must be >= 0. You specified ", userOptions$maxNbPtmsPerPeptide))
		q(status=-1)
	}

	#FILTER: maxNbPTMsPerPeptide
	userOptions$minPeptideLength <- cmdOpt$FLengthPeptide
	if(is.na(userOptions$minPeptideLength) | userOptions$minPeptideLength < 0 ){
		print(paste("ERROR. FLengthPeptide must be >= 0. You specified ", userOptions$minPeptideLength))
		q(status=-1)
	}

	#FILTER: FExclusivePeptides
	userOptions$FExclusivePeptides <- cmdOpt$FExclusivePeptides

	#FILTER: ratioCutOff
	userOptions$ratioCutOff <- cmdOpt$FRatioCutOff
	if(is.na(userOptions$ratioCutOff) | userOptions$ratioCutOff < 1){
		cat("ERROR. ratioCutoff must be > 1. You specified",userOptions$ratioCutOff,"\n")
		q(status=-1)
	}

	#FILTER: Exclude the PDF file as an output
	userOptions$noPDF <- cmdOpt$noPDF
	#FILTER: Output file
	userOptions$extra_info <- cmdOpt$extraOut

# FILTER (--F) END

# TMT
	userOptions$TAdjustRatios <- cmdOpt$TAdjustRatios


# TMT END


# STATISTICS

	#STATISTICS: median values
	userOptions$medianInfo <- cmdOpt$Median

	#STATISTICS: transform the i values to the log2
	userOptions$log2_i <- cmdOpt$log2toI

	#STATISTICS: normAC
	userOptions$normAC <- cmdOpt$SAnchorProtein

	#STATISTICS: SMissingValuesImutationMethod
	userOptions$SMissingValuesImutationMethod <- cmdOpt$SMissingValuesImutationMethod

	#STATISTICS: SNonPairWiseStatTest
	userOptions$SNonPairWiseStatTest <- cmdOpt$SNonPairWiseStatTest

	#STATISTICS: eBayes
	userOptions$eBayes <- cmdOpt$SPvalueInclude

	#STATISTICS: SRawDataAnalysis
	userOptions$SRawDataAnalysis <- cmdOpt$SRawDataAnalysis

# STATISTICS END

# EXPERIMENTAL DESIGN

	#EXPERIMENTAL DESIGN: EXperimentalDesign
	userOptions$expDesignTag <- cmdOpt$EXperimentalDesign

	userOptions$proteinQuant <- cmdOpt$EProteinQuant
	#userOptions$proteinQuant <- userOptions$selectedModifName != "."

	userOptions$ECorrelatedSamples <- cmdOpt$ECorrelatedSamples

# EXPERIMENTAL DESIGN END

# PDF-REPORT (--P)

	# PDF-REPORT: deFdrCutoff
	userOptions$deFdrCutoff <- cmdOpt$PQvalueCutOff
	if(is.na(userOptions$deFdrCutoff) | userOptions$deFdrCutoff <= 0 | userOptions$deFdrCutoff > 1 ){
		cat("ERROR. deFdrCutoff must be in the range [0-1]. You specified",userOptions$deFdrCutoff,"\n")
		q(status=-1)
	}

	# PDF-REPORT: PSelectedGraphics
#	userOptions$isDispExpDesign <- !regexpr("e",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isFdrPlots <- !regexpr("f",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isIntensityDistributionPlots <- !regexpr("i",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isVolcanoPlots <- !regexpr("v",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isHClustPlot <- !regexpr("h",cmdOpt$PSelectedGraphics) > -1
#	userOptions$isDeFdrPlot <- !regexpr("d",cmdOpt$PSelectedGraphics) > -1


# PDF-REPORT (--P) END

# TSV-REPORT (--T)
#
#	# TSV-REPORT: proteinFastaFile
#	userOptions$proteinFastaFile <- NA
#	if(nchar(cmdOpt$TFastaFile) > 0 ){
#		### check if file exists
#		if(file.exists(cmdOpt$TFastaFile)){
#			userOptions$proteinFastaFile <- cmdOpt$TFastaFile
#		}else{
#			cat("ERROR. File does not exist",cmdOpt$TFastaFile,"\n")
#			q(status=-1)
#		}
#	}
#
#	# TSV-REPORT: proteaseTarget	(deprecated)
#	userOptions$protease <- cmdOpt$TProtease

# TSV-REPORT (--T) END

# ADDITIONAL-REPORTS (--A)

	#ADDITIONAL-REPORTS iBaq
	userOptions$iBAQ <- cmdOpt$AIbaq

    #ADDITIONAL-REPORTS top3
	userOptions$top3 <- cmdOpt$ATop3

	#ADDITIONAL-REPORTS additional QC plots
	userOptions$addQC <- cmdOpt$AQC

	#ADDITIONAL-REPORTS rDataFile, isSaveRObject
	userOptions$isSaveRObject <- cmdOpt$ARDataFile
	#userOptions$rDataFile <- paste(userOptions$outputDir,userOptions$resultsFileLabel,".rData",sep="")


# ADDITIONAL-REPORTS (--A) END


# TEST

	### test run to define parameters (peptide analysis specific)
	userOptions$test <- cmdOpt$test

# TEST END
	return(userOptions)

}


#userInputTag <- "1,2,3:4,5,6" or experimental design file
# tag: 1,2:3:4,5,6
#condition isControl
#1 Condition 1      TRUE
#2 Condition 1      TRUE
#3 Condition 1     TRUE
#4 Condition 2     FALSE
#5 Condition 2     FALSE
#6 Condition 2     FALSE

# Experimental Design File Format should have the following 4 columns:
# Name: the names of the column where the quantify values are located
# Experiment: Condition Replication
# Groups: General Condition
# File: C for control files and T for files which should included in the analysis


#' Create experimental design data.frame from user input string
#' @param tag tag
#' @param expDesignDefault data.frame
#' @return data.frame describing experimental design
#' @export
#' @note  No note
#' @details  tag: 1,2:3:4,5,6
#'		condition isControl
#'	1 Condition 1 TRUE
#'	2 Condition 1 TRUE
#'	3 Condition 1 TRUE
#'	4 Condition 2 FALSE
#'	5 Condition 2 FALSE
#'	6 Condition 2 FALSE
#' @references NA
#' @examples print("No examples")
expDesignTagToExpDesign <- function(tag, expDesignDefault, ...){

  if (str_detect(tag, ".txt$")){
    expDesign <- expDesignTagToExpDesign_file_2(tag, ...)

    #(expDesign) <- rownames(expDesignDefault)[as.numeric(rownames(expDesign))]
  }else{
	sampleOrder <- as.numeric(unlist(strsplit(tag,"[\\,\\:]")))
	# make sure no duplicates, within range etc.
	if(is.na(sampleOrder[1])
			| (max(table(sampleOrder))>1)
			| (min(sampleOrder) < 1)
			| max(as.numeric(sampleOrder)) > nrow(expDesignDefault)
			){
		stop("ERROR: expDesignTagToExpDesign, INVALID EXPERIMENTAL DESIGN ",tag,"\n")

	}
	expDesign <- data.frame(row.names=sampleOrder,condition=rep(NA,length(sampleOrder)), isControl=rep(FALSE,length(sampleOrder))  )
	condNb <- 1
	for(cond in unlist(strsplit(tag,":"))){

		#cat(as.character(unlist(strsplit(cond,","))), paste("Condition",condNb) , "\n")
		expDesign[as.character(unlist(strsplit(cond,","))),]$condition <- paste("Condition",condNb,sep="")
		condNb <- condNb + 1
	}

	expDesign[ expDesign[,1] == "Condition1" ,]$isControl <- T
	expDesign[,1] <- as.factor(expDesign[,1]) ### has to be factor and not character
	### get original sample names

	rownames(expDesign) <- rownames(expDesignDefault)[as.numeric(rownames(expDesign))]
	#expDesignUser$condition <- expDesign[rownames(expDesignUser) ,]$condition

	# get original condition names, unless conditions have been split

	# if more conditions than originally use cond_1, cond_n labellinf @TODO can be done better, to avoid loosing org condition names
	if(length(unique(expDesign$condition)) > length(unique(expDesignDefault$condition))) return(expDesign)

	# check if conditions are split in new expDesign
	for(cond in unique(expDesign$condition)){
		runs <- rownames(expDesign)[expDesign$condition == cond]

		if(length(unique(expDesignDefault[runs,]$condition)) > 1){

			return(expDesign)
			#stop("SPLIT")
		}
	}

	### make sure a single condirion has not been split into two
	# I.e there should be just as many conditions before and anfter retreival of orginial condition names
	if(length(unique(	expDesign$condition)) == length(unique(	 expDesignDefault[rownames(expDesign),]$condition)) ){
		# get original condition names
		expDesign$condition <- expDesignDefault[rownames(expDesign),]$condition

	}
  }
	return(expDesign)

}


#' Create experimental design data.frame from user input string
#' @param tag tag
#' @param expDesignDefault data.frame
#' @return data.frame describing experimental design
#' @export
#' @note  No note
#' @details  Example Files for the different file formats
#' MaxQuant File:
#' Name	Experiment	Groups	Type
#' LFQ intensity 6892 - Calcium -1	6892 - Calcium -1	6892 - Calcium	C
#' LFQ intensity 6892 + Calcium -1	6892 + Calcium -1	6892 + Calcium	T
#' LFQ intensity 7648 - Calcium -1	7648 - Calcium -1	7648 - Calcium	T
#' LFQ intensity 7648 - Calcium -2	7648 - Calcium -2	7648 - Calcium	T
#'
#' @references NA
#' @examples print("No examples")
expDesignTagToExpDesign_file_2 <- function(t, ...){
  f_t <- list(...)
  f_n <- names(list(...))

  df_exp_design <- read.table(t, header = TRUE,  sep = "\t",  quote = "",  as.is = TRUE)

  #df_exp_design <- read.table(t,header = TRUE,  sep = "\t",  quote = "",  as.is = TRUE, fileEncoding = "UTF-16LE")
  expD <- df_exp_design[df_exp_design$Type %in% c('C', 'T'), c('Groups', 'Type', 'Name', 'Experiment')]

  names(expD) <- c('condition', 'isControl', 'file', 'c_Name')
  expD <- expD[
    with(expD, order(isControl),  decreasing = c(FALSE)),
    ]
  expD$isControl[expD$isControl == 'T'] <- FALSE
  expD$isControl[expD$isControl == 'C'] <- TRUE
  expD$isControl <- as.logical(expD$isControl)
  if ('file' %in% f_n){
    if (f_t$fileT == "MaxQuantProteinGroup"){
      res <- read.csv(f_t$file,allowEscapes=T, check.names=F,sep="\t")
      r <- colnames(res)
      ### get
      # Note 20210906: may need to ad method control here
      s_r <- str_replace_all(r[grepl("^LFQ",r)], "LFQ intensity ", "" )
      expD$f_pos <- match(expD$c_Name,s_r)

      expD$o_pos <- rownames(expD)
      rownames(expD) <- expD$f_pos

      expDesign <- expD[order(-expD$isControl, as.numeric(expD$o_pos)),]

    }else if(f_t$fileT == "SpectronautProteinGroup"){
      f_info <- gsub(".*\\.", "", f_t$file)
      if (f_info == "csv"){
        #res <- read.csv(f_t$file,allowEscapes=T, check.names=F,sep=",")
        res <- read.csv(f_t$file,allowEscapes=T, check.names=F,sep="\t", skip=1)
      }else{
        res <- read_excel(f_t$file)
      }

      r <- colnames(res)
      ### get
      # Note 20210906: may need to ad method control here
      if (sum(grepl(".*MS2Quantity$", r)) > 0){
        #s_r <- str_replace_all(r[grepl(".MS2Quantity",r)], "MS2Quantity", "")
        s_r <- r[grepl(".MS2Quantity",r)]
      }else if (sum(grepl(".*Quantity$", r)) > 0){
        #s_r <- str_replace_all(r[grepl(".Quantity",r)], "Quantity", "")
        s_r <- r[grepl(".Quantity",r)]
      }

      for (i in 1:length(rownames(expD))){
        pos <- sub("\\[(\\d+)\\].*", "\\1",expD[i, 'file'])
        expD[i, 'f_pos'] <- pos
      }

      expD$o_pos <- rownames(expD)
      rownames(expD) <- expD$f_pos
      expDesign <- expD[order(-expD$isControl, as.numeric(expD$o_pos)),]

    }else if (f_t$fileT %in% c("DiaNNProteinGroup", "DIANN_Peptide")){

      res <- read.table(file = f_t$file, sep = '\t', header = TRUE, quote = "\"")
      #res <- read.table(file = f_t$file)
      r <- colnames(res)
      if (sum(grepl("^Quantity\\.", r)) > 0){
        #s_r <- str_replace_all(r[grepl(".Quantity",r)], "Quantity", "")
        s_r <- r[grepl("^Quantity\\.",r)]
      }

      for (i in 1:length(rownames(expD))){
        if (gsub('-', '.', expD[i, 'file']) %in% s_r){
          expD[i, 'f_pos'] <- match(gsub('-', '.', expD[i, 'file']), s_r)
        }else{
          stop("ERROR: expDesignTagToExpDesign, INVALID COLUMN NAME IN THE EXPERIMENTAL FILE ",expD[i, 'file'],"\n")
        }

      }
      expD$o_pos <- rownames(expD)
      rownames(expD) <- expD$f_pos
      expDesign <- expD[order(-expD$isControl, as.numeric(expD$o_pos)),]
    }
  }else{
    expDesign <- expD[order(-expD$isControl),]

  }
  return(expDesign)

}

