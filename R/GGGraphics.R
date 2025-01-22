
#' Plots volcano, data points colored by max cv of the 2 compared conditions
#' @param data data.frame
#' @param log2RatioThrs default log2(0.5)
#' @param pValueThrs default 0.01
#' @param thrsLineCol default "lightgrey"
#' @param defalut 2
#' @param xlab default "log2 ratio"
#' @param ylab default "-log10 pValue"
#' @param title default no title
#' @param xlim xlim
#' @param ylim ylim
#' @param abline c("none","both","ratio","pvalue")
#' @param topNlabels default 10, label top proteins/peptides ordered by p-value
#' @param textSize default 20
#' @return ggplot2 object
#' @import ggplot2 ggrepel
#' @export
#' @note  No note
#' @details data.frame input object should contain columns ("ratio","pValue","geneName","ac","cv", "description")
#' @references NA
#' @examples print("No examples")
ggVolcanoPlot = function(data=data
		, title=""
		, pValueThrs=0.05
		, log2RatioThrs=0.5849625
		, thrsLineCol = "lightgrey"
		, thrsLineLty = 2
		, xlab = "log2 ratio"
		, ylab = "-log10 pValue"
		, textSize = 20
		, xlim = range(data$ratio,na.rm=T)
		, ylim = range(-log10(data$pValue),na.rm=T)
		, abline = c("both")
		, topNlabels = 10
){

  # plotted data
	p =  ggplot(data,aes(x = ratio,y=-log10(pValue)
					,label=pValue
					,label2=geneName
					,label3=ac
					,lable4=description
					,color=intensity))

	# axis
	p = p + labs(x=xlab, y=ylab, title=title)
	p = p + theme_bw()
	p = p + scale_x_continuous(limits =xlim)
	p = p + scale_y_continuous(limits =ylim)

	# abline
	#	 pvalue thrs
	if(abline %in% c("pvalue","both") ) p =  p + geom_abline(intercept = -log10(pValueThrs),slope=0, lty=thrsLineLty, col=thrsLineCol)
	# ratio thrs
	if(abline %in% c("ratio","both") ) p =  p + geom_vline(xintercept=c(-log2RatioThrs,log2RatioThrs), lty=thrsLineLty, col=thrsLineCol)

	# point style
	p = p + geom_point()
	p = p + scale_colour_gradientn(colours=c("yellow","blue"), name="Max Int\n(rank)" ) +
	  geom_point(subset(data,naHighLightSel &  (pValue <= pValueThrs)  & (abs(ratio) >= log2RatioThrs) ), mapping= aes(shape = naCat), col="grey",size=3) +
	  scale_shape_manual(values = c(2, 6, 11),name="NA Features\n(>=50%) ")
	#p=plabs
	# theme
	p = p + theme(text = element_text(size=textSize)
			#, axis.text.x = element_text(angle=0, hjust=1)
			, legend.position="right"
			, legend.direction="vertical"
			, legend.title = element_text(size=textSize*0.8)
			, legend.text=element_text(size=textSize*0.6)

	)

	# add labels
	# geom_GeomTextRepel() has yet to be implemented in plotly (status at v. 4.5.6).
	# disp geneName above thrs
	labPvalueThrs = ifelse(topNlabels > 0,sort(data$pValue)[min(topNlabels,length(data$pValue))],0)
	dfLab = subset(data, (pValue <= min(labPvalueThrs,pValueThrs) ) & (abs(ratio) >= log2RatioThrs))
	dfLab = dfLab[order(dfLab$pValue,decreasing = F), ]
	dfLab = dfLab[1:min(10,nrow(dfLab)),]
	dfLab = subset(dfLab, !is.na(geneName)) # avoid warning

	if(nrow(dfLab) > 0){
		p = p + geom_text_repel(data= dfLab,aes(x =ratio,y=-log10(pValue),label=geneName ))
	}
	return(p)
}


#' Plots volcano of all condition comparisons
#' @param sqa SafeQuantAnalysis object
#' @param isAdjusted (T/F) plot adjusted pvalues
#' @param see ggVolcanoPlot
#' @return ggplot2 object
#' @import ggplot2 ggrepel
#' @export
#' @note  No note
#' @details data.frame input object should contain columns ("ratio","pValue","geneName","ac","cv", "description")
#' @references NA
#' @examples print("No examples")
plotAllGGVolcanoes = function(sqa, isAdjusted=T ,...){
  
  # need more than own data feature
  if(nrow(sqa$eset) <= 1 ) return()

  # plot all volcanoes
  ctrlCondition = pData(sqa$eset)$condition[pData(sqa$eset)$isControl][1] %>% as.character
  caseConditions = setdiff(pData(sqa$eset)$condition %>% unique, ctrlCondition)

  if(isAdjusted){
    allPValue = sqa$qValue
  }else{
    allPValue = sqa$pValue
  }

  # avoid crash if no caseConditions
  if(length(caseConditions) > 0){
    xlim  = range(sqa$ratio, na.rm=T)
    ylim = range(abs(log10(allPValue)), na.rm=T)
  }

  # highlight features with >50% NA Int
  aboveNAThrs = getNAFraction(sqa$eset,method=c("cond","count")) >= 0.5

  for(cond in caseConditions){

    # create naCat factors
    aboveNAThrsLoc = aboveNAThrs[,c(ctrlCondition,cond)]
    aboveNAThrsCount = rowSums(aboveNAThrsLoc)
    naCat = ifelse(aboveNAThrsCount == 2,"both","none" )
    naCat[(naCat == "none") & aboveNAThrsLoc[,1]  ]  = ctrlCondition
    naCat[(naCat == "none") & aboveNAThrsLoc[,2]  ]  = cond
    naCat = factor(naCat, levels = c(ctrlCondition,cond,"both","none"))

    # compile df
    #cv = apply(sqa$cv[, c(ctrlCondition,cond) ],1,max, na.rm=T)*100
    #cv = apply(getNAFraction(sqa$eset,"cond")[, c(ctrlCondition,cond)],1, max, na.rm=T)*100
    intensity = apply(getSignalPerCondition(sqa$eset)[, c(ctrlCondition,cond)],1, max, na.rm=T) %>% rank
    ggDf = data.frame(ratio = sqa$ratio[,cond]
                      , pValue=allPValue[,cond]
                      , geneName = fData(sqa$eset)$geneName
                      , ac=fData(sqa$eset)$ac
                      , intensity = intensity
                      , description=fData(sqa$eset)$proteinDescription
                      , naHighLightSel = (aboveNAThrsCount > 0)
                      , naCat = naCat

                      )
    #plot
    plot(ggVolcanoPlot(data=ggDf, xlab = paste("log2", cond,"/",ctrlCondition ),xlim=xlim,ylim=ylim,  ... ))

  }
}

#' Plots the trends of selected proteins from all condition
#' @param sqa SafeQuantAnalysis object
#' @param isAdjusted (T/F) plot adjusted pvalues
#' @param proteinList a list of selevted gene names e.g. c("cbbTC", "cbbL2")
#' @param see ggVolcanoPlot
#' @return ggplot2 object
#' @import ggplot2 ggrepel
#' @export
#' @note  No note
#' @details data.frame input object should contain columns ("ratio","pValue","geneName","ac","cv", "description")
#' @references NA
#' @examples print("No examples")
plotSelectedProteins <- function(sqa, proteinList , isAdjusted=T , ...){
  fun_arg <- as.list(match.call())
  if("title_n" %in% names(fun_arg)){
    title_name <- fun_arg$title_n
  }else{
    title_name <- "Protein groups distribution - after imputation"
   exprs(sqa$eset) <- log2(exprs(sqa$eset))
  }


  t_df_line <- cbind(exprs(sqa$eset),fData(sqa$eset)[,c("proteinName", "geneName")])
  t_df_line_melt = melt(t_df_line)
  t_df_line_melt$variable <- pData(sqa$eset)[t_df_line_melt$variable, 'c_Name']
  t_exp_type=as.character(t_df_line_melt$variable)
  t_exp_type=substr(t_exp_type,1,nchar(t_exp_type)-2)
  t_df_line_melt=cbind(t_df_line_melt,t_exp_type)

  df_sel <- read.csv(proteinList,header = F,  sep = ',', row.names = NULL, quote = "", as.is = TRUE)
  df_sel <- df_sel[,1]

  t_df_line_sel = t_df_line_melt[t_df_line_melt$geneName %in% df_sel ,]

  # t_df_line_sel_melt = melt(t_df_line_sel)
  # t_exp_type=as.character(t_df_line_sel_melt$variable)
  # t_exp_type=substr(t_exp_type,1,nchar(t_exp_type)-2)
  # t_df_line_sel_melt=cbind(t_df_line_sel_melt,t_exp_type)
  t_df_line_sel_melt <- t_df_line_sel
  #
  #
  # 'add boxplots!'
  p2_line_raw = ggplot(data=t_df_line_melt, aes(x=variable, y=value, group=geneName)) + geom_tile(aes(fill = t_exp_type), width = 0.90, height = Inf, alpha = 0.8) +
    geom_line( size=1, alpha=0.1)+  geom_point( size=1, alpha=0.1)+ scale_fill_manual(values = qualitative_hcl(length(levels(pData(sqa$eset)$condition)),palette = "Pastel1")) +
    labs(title= title_name,subtitle= "log2 intensities ", x="Sample", y = "Intensity")+
    geom_point(data=subset(t_df_line_sel_melt), size=4, aes(color=geneName),alpha = 0.9)+
    geom_line(data=subset(t_df_line_sel_melt), size=2,aes(color=geneName),alpha = 0.7)+
    #geom_boxplot(data=t_df_line_melt, aes(variable, value),width=0.3, fill="white", size=1)  +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90))
  plot(p2_line_raw)
}

#' Plots the trends of selected proteins from all condition
#' @param sqa SafeQuantAnalysis object
#' @param isAdjusted (T/F) plot adjusted pvalues
#' @param see ggVolcanoPlot
#' @return ggplot2 object
#' @import ggplot2 ggrepel
#' @export
#' @note  No note
#' @details data.frame input object should contain columns ("ratio","pValue","geneName","ac","cv", "description")
#' @references NA
#' @examples print("No examples")
plotSelectedProteins2 <- function(sqa, isAdjusted=T ,...){

  t_df_line <- cbind(exprs(sqa$eset),fData(sqa$eset)[,c("proteinName", "geneName")])
  t_df_line_melt = melt(t_df_line)
  t_df_line_melt$variable <- pData(sqa$eset)[t_df_line_melt$variable, 'c_Name']
  t_exp_type=as.character(t_df_line_melt$variable)
  t_exp_type=substr(t_exp_type,1,nchar(t_exp_type)-2)
  t_df_line_melt=cbind(t_df_line_melt,t_exp_type)

  df_sel <- read.csv('D:\\proteomics\\SQ_v2.4.1_testRuns\\DIANN\\R_script-selection.csv',header = F,  sep = ',', row.names = NULL, quote = "", as.is = TRUE)
  df_sel <- df_sel[,1]

  t_df_line_sel = t_df_line_melt[t_df_line_melt$geneName %in% df_sel ,]

  # t_df_line_sel_melt = melt(t_df_line_sel)
  # t_exp_type=as.character(t_df_line_sel_melt$variable)
  # t_exp_type=substr(t_exp_type,1,nchar(t_exp_type)-2)
  # t_df_line_sel_melt=cbind(t_df_line_sel_melt,t_exp_type)
  t_df_line_sel_melt <- t_df_line_sel
  #
  #
  # 'add boxplots!'
  p2_line_raw = ggplot(data=t_df_line_melt, aes(x=variable, y=value, group=geneName)) +
    facet_grid(~ t_exp_type, scale = 'free_x') +
    geom_line( size=1, alpha=0.1)+  geom_point( size=1, alpha=0.1)+
    labs(title="Protein groups ",subtitle= "log2 intensities - filtered of con/sit/rev + half valid", x="Sample", y = "Intensity")+
    geom_point(data=subset(t_df_line_sel_melt), size=4, aes(color=geneName),alpha = 0.9)+
    geom_line(data=subset(t_df_line_sel_melt), size=2,aes(color=geneName),alpha = 0.7)+
    #geom_boxplot(data=t_df_line_melt, aes(variable, value),width=0.3, fill="white", size=1)  +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90))
  plot(p2_line_raw)
}

#' Plots violin for all the conditions
#' @param eset SafeQuantAnalysis object
#' @param fType defines from which software the data originally were analyzed 
#' @param see ggVolcanoPlot
#' @return ggplot2 object
#' @import ggplot2 ggrepel
#' @export
#' @note  No note
#' @details data.frame input object should contain columns ("ratio","pValue","geneName","ac","cv", "description")
#' @references NA
#' @examples print("No examples")
plot_violin <- function(eset, fType,... ){
  if (fType %in% c("DIANN_Peptide")){
    data_df <- as.data.frame(exprs(eset))
  }else{
    data_df <- as.data.frame(log2(exprs(eset)))
  }
  
  data_df$ID <- rownames(data_df)
  data_df <- melt(data_df, id = "ID", variable.name = "File", value.name = "Quantity")
  data_df$Group <- rep(pData(eset)$condition, each = length(rownames(exprs(eset))))
  if ("file" %in% colnames(pData(eset))){
    data_df$file2 <- str_replace(str_replace(str_replace(rep(pData(eset)$file, each = length(rownames(exprs(eset)))), "Quantity.", ""),".raw", ""), ".PG.Quantity", "")
    # Note 20221212: Temporary solution for the x-axis order problem
    return(ggplot(data_df, aes(x=file2, y = Quantity, fill = Group)) +
             geom_violin() +
             geom_boxplot(width = 0.1) +
             scale_x_discrete(limits = str_replace(str_replace(unique(data_df$file2), "Quantity.", ""), ".raw", "")) +
             theme_bw() +
             theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  }else{
    data_df$file2 <- data_df$File
    # Note 20221212: Temporary solution for the x-axis order problem
    return(ggplot(data_df, aes(x=file2, y = Quantity, fill = Group)) +
             geom_violin() +
             geom_boxplot(width = 0.1) +
             # scale_x_discrete(limits = data_df$file2) +
             theme_bw() +
             theme(axis.text.x = element_text(angle = 90, hjust = 1)))
  }
  

}
