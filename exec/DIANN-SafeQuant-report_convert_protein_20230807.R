# rm(list=setdiff(ls(), ""))
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# getwd()


suppressWarnings(suppressPackageStartupMessages(library("readxl")))
suppressWarnings(suppressPackageStartupMessages(library ("openxlsx")))
suppressWarnings(suppressPackageStartupMessages(library ("ggplot2")))
suppressWarnings(suppressPackageStartupMessages(library("dplyr")))
suppressWarnings(suppressPackageStartupMessages(library("hash")))
suppressWarnings(suppressPackageStartupMessages(library(ggrepel)))
suppressWarnings(suppressPackageStartupMessages(library(dplyr)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(magrittr)))
suppressWarnings(suppressPackageStartupMessages(library(optparse)))

option_list <- list(
  ### I/O
  make_option(c("-i", "--inputFile"), type="character", default="",
              help="I/O:  Input file: DIANN file report.tsv",
  ),
  make_option(c("-l", "--resultsFileLabel"), type="character", default="report_SQ_final_convert.tsv",
              help="I/O: results file name [default %default]",
  ),
  make_option(c("--QValuefilter"), type="character", default="Global",
              help="FILTER: --QV 'Global', 'Protein' [default %default]",
              metavar=" Global: Global.PG.Q.Value
                            Protein: Protein.Q.Value
            ")

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
getUserOptions <- function(){
  # alert <- combine_styles("bold", "red4", "bgCyan")
  epilogue <- "Examples:
	>Rscript DIANN-SafeQuant-report_convert_protein_20230807.R  -i /path/to/report.tsv
	>Rscript DIANN-SafeQuant-report_convert_protein_20230807.R  -i /path/to/report.tsv --QV Protein
	"
  cmdOpt <- parse_args(OptionParser( prog=paste("DIANN-SafeQuant-report_convert_protein"), option_list=option_list, epilogue=epilogue))

  userOptions <- list()
  userOptions$inputFile <- cmdOpt$inputFile
  if( userOptions$inputFile == "" | !file.exists(userOptions$inputFile)){
    cat("ERROR. Please specify input file.",userOptions$inputFile, "Not found!","\n")
    q(status=-1)
  }

  userOptions$resultsFileLabel <- cmdOpt$resultsFileLabel
  if (cmdOpt$QValuefilter == "Global"){
    userOptions$QValuefilter <- "Global.PG.Q.Value"
  }else if (cmdOpt$QValuefilter == "Protein"){
    userOptions$QValuefilter <- "Protein.Q.Value"
  }else{
    userOptions$QValuefilter <- "Global.PG.Q.Value"
  }


  return(userOptions)

}



getPeptides <- function(proteinSeq,proteaseRegExp="[KR]",nbMiscleavages=0){

  allAA <- as.vector(unlist(strsplit(proteinSeq,"")))

  ### get all full cleaved peptides without cleaved residue
  fcPeptides <- strsplit(proteinSeq,proteaseRegExp,perl=TRUE)[[1]]

  ### add cleaved residue
  matchPos <- gregexpr(proteaseRegExp,proteinSeq,perl=TRUE)[[1]]
  separator <- allAA[ matchPos ]
  ###if c-term peptide isn't tryptic
  if(length(separator) < length(fcPeptides)) separator <- c(separator,"")
  fcPeptides <- paste(fcPeptides,separator,sep="")
  fcPeptides <- fcPeptides[nchar(fcPeptides) > 0 ]

  ### if no mis-cleavages that's it
  if(nbMiscleavages == 0){
    return(fcPeptides)
  }

  allPeptides <- c()
  ### handle miscleavages
  for(i in 1:length(fcPeptides)){
    #cat(fcPeptides[i]," ",i," -------------\n")
    for(j in i:(i+nbMiscleavages)){
      if(j <= length(fcPeptides)){
        pept <- paste(fcPeptides[i:j],collapse="")
        #cat(j," ---",pept,"\n")
        allPeptides <- c(allPeptides,pept)
      }
    }
  }

  return(allPeptides[nchar(allPeptides) > 0])

}

sum_calc <- function(d){
  if (d %in% df2$Protein.Group){
    return(sum(df2[df2$Protein.Group == d, nrPrecursor]))
  }else{
    return("0")
  }

}
sum_calc2 <- function(d){
  if (d %in% all_prot_convert_df2$Protein.Group){
    return(sum(all_prot_convert_df2[all_prot_convert_df2$Protein.Group == d, "NrOfuniquePeptides"]))
  }else{
    return("0")
  }

}

sum_calc_pept <- function(d){
  if (d %in% df2$Protein.Group){
    g_sm_df <- df2[df2$Protein.Group == d,quantFiles]

    g_sm_df[g_sm_df > 1] <- 1
    g_sm_df[g_sm_df == 0] <- 0


    all_Files <- colSums(g_sm_df)

    return(sum(all_Files))
  }else{
    return("0")
  }

}

i_found <- function(d){
  if(d %in% df2$Protein.Group){
    return(max(df2[df2$Protein.Group == d, quantFiles]))
  }else{
    return("0")
  }

}

i_found2 <- function(d){
  if(d %in% all_prot_convert_df2$Protein.Group){
    return(max(all_prot_convert_df2[all_prot_convert_df2$Protein.Group == d, quantFiles]))
  }else{
    return("0")
  }

}

allSame <- function(x) length(unique(x)) == 1

randFilter <- 'Include'

userOptions <- getUserOptions()
cols <- c("File.Name", "Run", "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
          "PG.Quantity", "Modified.Sequence", "Stripped.Sequence", "Precursor.Id", "Precursor.Charge",
          "Global.PG.Q.Value", "PG.Q.Value", "Precursor.Quantity", "RT", "First.Protein.Description")
dt <- fread(userOptions$inputFile, select = cols)


dt$File.Name <- basename(dt$File.Name)
dt[["Precursor.Quantity"]] %<>% as.numeric()
dt %<>% magrittr::extract(order(Protein.Group, Run, -get("Precursor.Quantity")))

# dt <- dt[Global.PG.Q.Value <= 0.01,]
dt <- dt[dt[[userOptions$QValuefilter]] <= 0.01,]


uniq_PG <- sort(unique(dt$Protein.Group))



# Note: The Protein.Ids was used instead of the Protein.Group because it was observed that there possibilities
# where there are multible unique proteins for the Protein.Ids and only one option for the Protein.Group
# This is mostly for the contamination which it will be filter out in a later step and a few rare cases with the main sample.
# This it was one of the biggest issue from Beau dataset and for this reason there was an option that was disselect from the DIANN run and then all proteins mention for both Columns were present
# Which means for the bel,ow command it is not necessary to split by the Protein.Ids and used Protein.Group instead.
dt[, PG.Sum := sum(sort(get('Precursor.Quantity'))),by = c('Protein.Group',  'Run')]
dt[, NrOfPrecursorsMeasured := .N, by = c('Protein.Group', 'Run')]
dt[, NrOfuniquePeptides := length(unique(get('Precursor.Id'))), by = c('Protein.Group')]


cols <- c('Protein.Group', 'Protein.Names', 'Genes', 'Run', 'Protein.Ids', 'NrOfuniquePeptides',
         'First.Protein.Description', 'Global.PG.Q.Value', 'PG.Sum', 'NrOfPrecursorsMeasured')
dt_PG <- dt %>% magrittr::extract(, cols, with = FALSE)

dt_PG_unique <- dt_PG %>% unique()


precursors <- data.table::dcast(dt_PG_unique, Protein.Group +
                                  Protein.Ids + NrOfuniquePeptides ~ Run, value.var = 'NrOfPrecursorsMeasured')
sum_quantity <- data.table::dcast(dt_PG_unique, Protein.Group +
                                    Protein.Ids + NrOfuniquePeptides ~ Run, value.var = 'PG.Sum')


# Get the column names
precursors_names <- colnames(precursors)
first_names <- c("Protein.Group", "Protein.Ids", "NrOfuniquePeptides")
replace_names <- precursors_names[!(precursors_names %in% first_names)]
precursors_n <- paste("NrOfPrecursorsMeasured", replace_names, sep = ".")

precursors_names <- c(first_names, precursors_n)
colnames(precursors) <- precursors_names

quantity_n <- paste("Quantity", replace_names, sep = ".")
quantity_names <- c(first_names, quantity_n)
colnames(sum_quantity) <- quantity_names

all_prot_convert_df <- as.data.frame(merge(precursors, sum_quantity, by = c("Protein.Group", "Protein.Ids", "NrOfuniquePeptides")))
all_prot_convert_df2 <- as.data.frame(merge(precursors, sum_quantity, by = c("Protein.Group", "Protein.Ids", "NrOfuniquePeptides")))
all_prot_convert_df2[is.na(all_prot_convert_df2)] <- 0


uniq_File <- unique(dt$File.Name)

all_prot_convert_df$Protein.Group.orig <- all_prot_convert_df$Protein.Group
all_prot_convert_df <- separate_rows(all_prot_convert_df, Protein.Group, sep=";", convert= TRUE)


df2 <- all_prot_convert_df
df2[is.na(df2)] <- 0
names_col <- colnames(df2)
nrPrecursor <- names_col[grepl("NrOfPrecursorsMeasured", names_col)]
quantFiles <- names_col[grepl("Quantity", names_col)]


df_redandant <- df2[grepl(";", df2$Protein.Ids), ]
df_redandant$majorP <- ''
unq_Peptide <- unique(df_redandant$Protein.Ids)


# When I add the Gene column the below part is not function correctly
for (i in unq_Peptide){
  name_i <- strsplit(i, ";")
  sm_df <- data.frame(matrix("", ncol = 1, nrow = length(name_i[[1]])))
  colnames(sm_df) <- c("Protein")
  sm_df$Protein <- name_i[[1]]
  # sm_df$un_Peptide <- mapply(sum_calc, sm_df$Protein)
  sm_df$un_Peptide <- mapply(sum_calc2, sm_df$Protein)
  # sm_df$i_max <- mapply(i_found, sm_df$Protein)
  sm_df$i_max <- mapply(i_found2, sm_df$Protein)

  if (allSame(sm_df$un_Peptide)){
    # break
    if (unique(sm_df$un_Peptide == 0)){
       df_redandant[df_redandant$Protein.Ids == i, ]$majorP <- df_redandant[df_redandant$Protein.Ids == i,]$Protein.Group.orig
       df_redandant[df_redandant$Protein.Ids == i, ]$Protein.Group <- df_redandant[df_redandant$Protein.Ids == i,]$Protein.Group.orig
       df2[df2$Protein.Ids == i,]$Protein.Group <- df_redandant[df_redandant$Protein.Ids == i,]$Protein.Group.orig
      next
    }else if(allSame(sm_df$i_max)){
      df_redandant[df_redandant$Protein.Ids == i, ]$majorP <- df_redandant[df_redandant$Protein.Ids == i, ]$Protein.Group.orig
      replace_name_i <- df_redandant[df_redandant$Protein.Ids == i, ]$Protein.Group.orig[1]
      majorP <- df_redandant[df_redandant$Protein.Ids == i, ]$Protein.Group[1]
      g <- df2[df2$Protein.Ids == i,]
      notmajor <- sm_df[sm_df$Protein != majorP, "Protein"]
      df2 <- df2[!(df2$Protein.Group %in% notmajor & df2$Protein.Ids ==i),]
      df2[df2$Protein.Group == majorP & df2$Protein.Ids == i,]$Protein.Group <- replace_name_i
      df_redandant <- df_redandant[!(df_redandant$Protein.Group %in% notmajor & df_redandant$Protein.Ids ==i),]
      df_redandant[df_redandant$Protein.Group == majorP & df_redandant$Protein.Ids == i,]$Protein.Group <- replace_name_i
    }else{
      majorP <- sm_df[which.max(sm_df$i_max),'Protein']
      g <- df2[df2$Protein.Ids == i,]
      notmajor <- sm_df[sm_df$Protein != majorP, "Protein"]
      df2 <- df2[!(df2$Protein.Group %in% notmajor & df2$Protein.Ids ==i),]
      df_redandant[df_redandant$Protein.Ids == i, ]$majorP <- majorP
      df_redandant <- df_redandant[!(df_redandant$Protein.Group %in% notmajor & df_redandant$Protein.Ids ==i),]
    }
  }else{
    max_number_Peptides <- sm_df[which.max(sm_df$un_Peptide), 'un_Peptide']
    max_number_intensity <- sm_df[which.max(sm_df$un_Peptide), 'i_max']
    check_sm_df <- sm_df[sm_df$un_Peptide == max_number_Peptides & sm_df$i_max == max_number_intensity,]
    if (nrow(check_sm_df) == 1){
      majorP <- sm_df[which.max(sm_df$un_Peptide),'Protein']
      g <- df2[df2$Protein.Ids == i,]
      notmajor <- sm_df[sm_df$Protein != majorP,'Protein']
      df2 <- df2[!(df2$Protein.Group %in% notmajor & df2$Protein.Ids ==i),]
      df_redandant <- df_redandant[!(df_redandant$Protein.Group %in% notmajor & df_redandant$Protein.Ids ==i),]
      df_redandant[df_redandant$Protein.Ids == i, ]$majorP <- majorP
    }else{
      protein_group_orig <- unique(df_redandant[df_redandant$Protein.Ids == i, ]$Protein.Group.orig)
      replace_name_list <- check_sm_df$Protein
      replace_name_df <- paste(sort(replace_name_list), collapse= ";")
      origin_sort <- paste(sort(strsplit(protein_group_orig, ";")[[1]]), collapse = ";")
      if (replace_name_df == origin_sort){
        df2 <- df2[!(df2$Protein.Group %in% replace_name_list[-1] & df2$Protein.Ids ==i),]
        df2[df2$Protein.Ids == i,]$Protein.Group <- df2[df2$Protein.Ids == i, ]$Protein.Group.orig
        df_redandant <- df_redandant[!(df_redandant$Protein.Group %in% replace_name_list[-1] & df_redandant$Protein.Ids == i),]
        df_redandant[df_redandant$Protein.Ids == i,]$Protein.Group <- df_redandant[df_redandant$Protein.Ids == i, ]$Protein.Group.orig
        df_redandant[df_redandant$Protein.Ids == i, ]$majorP <- protein_group_orig
      }else{
        print("Looks like there is an error here please chek what is going on")
        print(i)
        break
      }
    }
  }
}
#

df2 <- unique(df2)
df_redandant <- unique(df_redandant)


df2_short <- df2[, c("Protein.Group", nrPrecursor, quantFiles)]

df_new <- aggregate(. ~ Protein.Group, df2_short, sum)
df_new$RazorG <- ""
df_new$RazorP <- ""
df_new$Genes <- ""
df_new$Protein.Names <- ""
df_new$Global.PG.Q.Value <- ""
df_new$First.Protein.Description <- ""
colNames_1 <- names_col[!(names_col %in% colnames(df_new))]
df_new[, colNames_1] <- ''
nrPrecursor_2 <- gsub('NrOfPrecursorsMeasured', 'NrOfRazorPrecursors', nrPrecursor)
df_new[nrPrecursor_2] <- ''

df_gene_info <- as.data.frame(unique(dt_PG_unique[, c("Protein.Group", "Genes", "First.Protein.Description", "Global.PG.Q.Value", "Protein.Ids", "Protein.Names")]))
colnames_sel <- c("Protein.Ids", "Genes", "Global.PG.Q.Value", "First.Protein.Description", "Protein.Names")
for (i in 1:nrow(df_new)){
  # Note: keep in mind that this it may cause problems
  if (df_new$Protein.Group[i] %in% df_redandant$majorP){
    df_new[i, colNames_1] <- df2[df2$Protein.Group == df_new$Protein.Group[i], colNames_1]
    if (nrow(df_redandant[df_redandant$majorP == df_new$Protein.Group[i], ]) == 1){
      df_new$Genes[i] <- df_gene_info[df_gene_info$Protein.Group == df_new$Protein.Group[i],]$Genes
      df_new$Protein.Names[i] <- df_gene_info[df_gene_info$Protein.Group == df_new$Protein.Group[i],]$Protein.Names
      df_new$Global.PG.Q.Value[i] <- df_gene_info[df_gene_info$Protein.Group == df_new$Protein.Group[i],]$Global.PG.Q.Value
      df_new$First.Protein.Description[i] <- df_gene_info[df_gene_info$Protein.Group == df_new$Protein.Group[i],]$First.Protein.Description
      df_new$RazorG[i] <- df_gene_info[df_gene_info$Protein.Ids == df_redandant[df_redandant$majorP == df_new$Protein.Group[i],]$Protein.Ids,]$Genes[1]
      df_new$RazorP[i] <- df_redandant[df_redandant$majorP == df_new$Protein.Group[i], ]$Protein.Ids
      razor_info <- df_redandant[df_redandant$majorP == df_new$Protein.Group[i], nrPrecursor]
      colnames(razor_info) <- nrPrecursor_2
      df_new[i, nrPrecursor_2] <- razor_info
  # break
    }else{
      # get all the gene names
      genesName_info <- df_redandant[df_redandant$majorP == df_new$Protein.Group[i], ]$Genes
      genesName_info <- paste(genesName_info, collapse=";")
      genesName_info <- unique(strsplit(genesName_info, ";")[[1]])
      genesName_info <- paste(genesName_info, collapse=";")
      df_new$RazorG[i] <- genesName_info
      # get all the protein ids
      proteinIds_info <- df_redandant[df_redandant$majorP == df_new$Protein.Group[i], ]$Protein.Ids
      proteinIds_info <- paste(proteinIds_info, collapse=";")
      proteinIds_info <- unique(strsplit(proteinIds_info, ";")[[1]])
      proteinIds_info <- paste(proteinIds_info, collapse=";")
      df_new$RazorP[i] <- proteinIds_info
      razor_info <- df_redandant[df_redandant$majorP == df_new$Protein.Group[i], nrPrecursor]
      colnames(razor_info) <- nrPrecursor_2
      razor_info <- colSums(razor_info)
      df_new[i, nrPrecursor_2] <- razor_info
# break
    }

  }else{
    colnames_sel <- c("Protein.Ids", "Protein.Names","Genes", "Global.PG.Q.Value", "First.Protein.Description")
    df_new[i, colnames_sel] <- df_gene_info[df_gene_info$Protein.Group == df_new$Protein.Group[i], colnames_sel]
    # break
  }

}
# df_new <- df_new[,c("Protein.Group", "Protein.Ids", "Protein.Names","RazorP", "Genes", "RazorG", "First.Protein.Description", "Global.PG.Q.Value", nrPrecursor, nrPrecursor_2, quantFiles)]
df_new <- df_new[,c("Protein.Group", "Protein.Ids", "Protein.Names","RazorP", "Genes", "RazorG", "First.Protein.Description", "Global.PG.Q.Value", nrPrecursor, nrPrecursor_2, quantFiles)]

write.table(df_new, file=paste(dirname(userOptions$inputFile),'/', userOptions$resultsFileLabel, sep = ""), quote=FALSE, sep='\t', col.names = NA)

cat("INFO: CREATED FILE ", paste(dirname(userOptions$inputFile),'/', userOptions$resultsFileLabel, sep = ""),"\n")
