
proteinTable_reformat <- function(filePath, df2, df3){
  df <- read.table(file = filePath, sep="\t", header = TRUE, quote="\"")

  df_sm <- df[,c("Protein.Ids", "Genes", "Razor_genes", "Global.PG.Q.Value", "missCl")]
  df_sm$razor <- ""

  df_sm[df_sm$Protein.Ids != "",]$razor <- 2

  # df_sm[grepl(";", df_sm$Protein.Ids),]$razor <- 2
  df_sm[df_sm$Protein.Ids == "",]$razor <- 1

  # df_sm[!(grepl(";", df_sm$Protein.Ids)),]$razor <- 1

  df3_sm <- df3[grepl('NrOfPrecursors', names(df3))]
  proteinNames <- df3$proteinName
  proteinIds <- df3$Protein.Ids
  df3_sm <- cbind(proteinNames, proteinIds, df3_sm)
  df3_sm$razor <- ""
  df3_sm[df3_sm$proteinIds !="",]$razor <- 2
  df3_sm[df3_sm$proteinIds == "",]$razor <- 1
  df3_sm_razor <- df3_sm[df3_sm$razor == 2,]
  df3_sm_Notrazor <- df3_sm[df3_sm$razor == 1,]
  df3_sm_razor <- aggregate(.~proteinNames, df3_sm_razor, function(x) sum(as.numeric(x)))
  names_ch <- names(df3_sm_razor)[grepl('NrOfPrecursors', names(df3_sm_razor))]
  names_ch_rpl <- str_replace(names_ch, "NrOfPrecursors", "NrOfRazorPeptides")
  names(df3_sm_razor)[names(df3_sm_razor) %in% names_ch] <- names_ch_rpl
  df3_sm_all <- aggregate(.~proteinNames, df3_sm, function(x) sum(as.numeric(x)))

  df_s <- df_sm %>%group_by(Genes) %>% arrange(desc(razor)) %>%  filter(row_number() == 1)
  # df_s <- df_s[df_s$razor == 2,]

  df2$Razor_genes <- ''
  df2$Razor_protein <- ''
  df2$Global.PG.Q.Value <- ''
  # df2$missCl <- ''
  df2[, names_ch_rpl] <- ''
  df2[, names_ch] <- ''
  for (i in 1:nrow(df2)){
    if (df2$geneName[i] %in% df_s$Genes){
      rGenes <- df_s[df_s$Genes == df2$geneName[i], ]$Razor_genes
      mGenes <- df2$geneName[i]
      r_gene_lst <- str_split(rGenes, ";")
      r_gene_lst <- r_gene_lst[[1]][!(grepl(mGenes, r_gene_lst[[1]]))]
      df2$Razor_genes[i] <- paste(r_gene_lst, collapse = ";")

      rIDs <- df_s[df_s$Genes == df2$geneName[i], ]$Protein.Ids
      mIDs <- df2$proteinName[i]
      r_ids_lst <- str_split(rIDs, ";")
      r_ids_lst <- r_ids_lst[[1]][!(grepl(mIDs, r_ids_lst[[1]]))]
      df2$Razor_protein[i] <- paste(r_ids_lst, collapse = ";")

      df2$Global.PG.Q.Value[i] <- df_s[df_s$Genes == df2$geneName[i], ]$Global.PG.Q.Value
      # df2$missCl[i] <- df_s[df_s$Genes == df2$geneName[i], ]$missCl
    }
    for (i2 in names_ch){
      df2[i, i2] <- df3_sm_all[df3_sm_all$proteinName == df2[i, 'proteinName'], i2]
    }
    if (df2[i, 'proteinName'] %in% df3_sm_razor$proteinNames){
      for (i2 in names_ch_rpl){
        df2[i, i2] <- df3_sm_razor[df3_sm_razor$proteinNames == df2[i, 'proteinName'], i2]
      }
    }
  }


  colnames(df2)[which(names(df2) == "Razor_genes")] <- "Other_geneNames"
  colnames(df2)[which(names(df2) == "Razor_protein")] <- "Other_proteins"

  df2 <- df2[,!(colnames(df2) %in% c("ac", "idScore", "nbPeptides", "allAccessions"))]

  neworder <- c("proteinName", "Other_proteins", "geneName", "Other_geneNames", "proteinDescription", "Global.PG.Q.Value", names_ch, names_ch_rpl)
  neworder <- c(neworder, names(df2)[!(names(df2) %in% neworder)])
  df2 <- df2[,neworder]
  names(df2)[names(df2) %in% names_ch] <- paste(sub('NrOfPrecursors\\.', '', names(df2)[names(df2) %in% names_ch]), "NrOfPrecursors", sep=".")
  names(df2)[names(df2) %in% names_ch_rpl] <- paste(sub('NrOfRazorPeptides\\.', '', names(df2)[names(df2) %in% names_ch_rpl]), "NrOfRazorPeptides", sep=".")

  return(df2)
}

# peptideTable_reformat <- function(df2){
#
#
#
#   df2$Razor_genes <- ''
#   df2$Razor_protein <- ''
#   df2$Global.PG.Q.Value <- ''
#   # df2$missCl <- ''
#   df2[, names_ch_rpl] <- ''
#   df2[, names_ch] <- ''
#   for (i in 1:nrow(df2)){
#
#       rGenes <- df_s[df_s$Genes == df2$geneName[i], ]$Razor_genes
#       mGenes <- df2$geneName[i]
#       r_gene_lst <- str_split(rGenes, ";")
#       r_gene_lst <- r_gene_lst[[1]][!(grepl(mGenes, r_gene_lst[[1]]))]
#       df2$Razor_genes[i] <- paste(r_gene_lst, collapse = ";")
#
#       rIDs <- df_s[df_s$Genes == df2$geneName[i], ]$Protein.Ids
#       mIDs <- df2$proteinName[i]
#       r_ids_lst <- str_split(rIDs, ";")
#       r_ids_lst <- r_ids_lst[[1]][!(grepl(mIDs, r_ids_lst[[1]]))]
#       df2$Razor_protein[i] <- paste(r_ids_lst, collapse = ";")
#
#       df2$Global.PG.Q.Value[i] <- df_s[df_s$Genes == df2$geneName[i], ]$Global.PG.Q.Value
#       # df2$missCl[i] <- df_s[df_s$Genes == df2$geneName[i], ]$missCl
#     }
#     for (i2 in names_ch){
#       df2[i, i2] <- df3_sm_all[df3_sm_all$proteinName == df2[i, 'proteinName'], i2]
#     }
#     if (df2[i, 'proteinName'] %in% df3_sm_razor$proteinNames){
#       for (i2 in names_ch_rpl){
#         df2[i, i2] <- df3_sm_razor[df3_sm_razor$proteinNames == df2[i, 'proteinName'], i2]
#       }
#     }
#   }
#
#
#   # print(names(df2)[names(df2) %in% names_ch])
#
#   colnames(df2)[which(names(df2) == "Razor_genes")] <- "Other_geneNames"
#   colnames(df2)[which(names(df2) == "Razor_protein")] <- "Other_proteins"
#
#   df2 <- df2[,!(colnames(df2) %in% c("ac", "idScore", "nbPeptides", "allAccessions"))]
#
#   neworder <- c("proteinName", "Other_proteins", "geneName", "Other_geneNames", "proteinDescription", "Global.PG.Q.Value", names_ch, names_ch_rpl)
#   neworder <- c(neworder, names(df2)[!(names(df2) %in% neworder)])
#   df2 <- df2[,neworder]
#   names(df2)[names(df2) %in% names_ch] <- paste(sub('NrOfPrecursors\\.', '', names(df2)[names(df2) %in% names_ch]), "NrOfPrecursors", sep=".")
#   names(df2)[names(df2) %in% names_ch_rpl] <- paste(sub('NrOfRazorPeptides\\.', '', names(df2)[names(df2) %in% names_ch_rpl]), "NrOfRazorPeptides", sep=".")
#
#   return(df2)
# }
