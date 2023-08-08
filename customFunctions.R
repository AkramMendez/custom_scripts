################################
# Enrichment analysis, Differential expression untreated CellType1 (KD vs Control)
#########################################################################################

myplotEnrichment <- function(pathway,stats,ticksSize ){
  
  p <- plotEnrichment(pathway = pathway, stats = stats,ticksSize =ticksSize ) + 
    scale_x_continuous(expand = expansion(mult = c(0,0)))
    #+
  #labs(title = paste(i,title, sep = "\n"), y="Enrichment Score (ES)", x="Gene rank")
  
  ymin.coord= ggplot_build(p)$layout$panel_params[[1]]$y.range[1]
  
  p <- p + scale_y_continuous(expand = expansion(mult = c(0.01,0.01)), limits = c(-1.0,1.0)) 
  
  tmp <- stats
  rnk <- rank(-tmp)
  ord <- order(rnk)
  statsAdj <- tmp[ord]
  #statsAdj <- sign(statsAdj) * (abs(statsAdj)^0.5)
  #statsAdj <- statsAdj/max(abs(statsAdj))
  
  tmp1.stats <- as.data.frame(statsAdj) %>% 
    mutate(rank=row_number()) %>% 
    mutate(inter=cut_interval(rank,n=50)) %>% 
    group_by(inter) %>% 
    summarise(logFC=mean(statsAdj)) %>% 
    dplyr::select(inter,logFC) %>% 
    ungroup() %>% 
    distinct() %>% 
    separate(col = "inter", into = c("xmin","xmax"),sep=",") %>% 
    mutate(xmin=as.numeric(gsub("(\\[|\\()","",xmin,perl = T)),
           xmax=as.numeric(gsub("(\\]|\\))","",xmax,perl = T)))%>% 
    mutate(ymin=ymin.coord)
  
  message("ymin=",max(tmp1.stats$ymin))
  
  p <- p +  geom_rect(data=tmp1.stats,aes(xmin=xmin,xmax=xmax,ymin=-0.95,ymax=-0.85, fill=logFC)) + 
    scale_fill_gradient2(low="navy",mid = "white", high = "firebrick",midpoint = 0) +
    theme_prism(base_size = 10, base_line_size = 0.5, border = T) +
    theme(legend.position = "none")
  
  return(p)
  
  
}

plotGSEAindiv <- function(rankedList, termsList,terms,ontologyName,title,outdir,outdirName){
  # Plot individually GSEA results for significant terms  
  for(i in terms){
    
    p1 <- myplotEnrichment(pathway = termsList[[i]], stats = rankedList, ticksSize = 0.2) + 
      #scale_x_continuous(expand = expansion(mult = c(0,0.01))) +
      labs(title = paste(i,title, sep = "\n"), y="Enrichment Score (ES)", x="Gene rank")
    
    # p1 <- p1  + geom_tile(data=rankedList %>%
    #               as.data.frame() %>%
    #               dplyr::rename(log2FoldChange=".") %>%
    #               mutate(rank=row_number()),
    #             aes(x=rank, y=-0.075,fill=cut(log2FoldChange,seq(1,length(rankedList),by=500))), 
    #             height=0.05,interpolate = T,inherit.aes = F, alpha=0.95) + 
    #   theme(legend.position = "none") +
    #   scale_fill_manual(values = coolwarm(15000))
    # 
    #p <- p1 / p2
    
    path <- paste0(outdir,outdirName,"/",ontologyName,"_",gsub("\\s+|\n","_",title))
    
    file <- paste0(path,"/",i,"_",gsub("\\s+|\n","_",title),".pdf")
    
    if(!dir.exists(paths = path))
    {
      message("Creating dir:",path)
      dir.create(path,recursive = TRUE)
      message("Saving plot to",file)
      ggsave(file, plot = p1,width = 8, height = 6, dpi = 300)
      
    }else{
      message("Saving plot to",file)
      ggsave(file, plot = p1,width = 8, height = 6, dpi = 300)
    }
    
  }
  
  message("Done")
}



######################################################################
# Plotting GSEA results
#####################################################################

plotGSEAbarplot <- function(fgseaResult,title){
  
  ntop <- 10 
  p <- fgseaResult %>%
    mutate(col=if_else(NES > 0,"#D85847","#5C7BE4")) %>%
    dplyr::group_by(type=sign(NES)) %>% arrange(dplyr::desc(NES)) %>% filter(row_number() <= ntop) %>% ungroup() %>%
    mutate(label=pathway) %>% 
    #mutate(label=if_else(grepl("HOX",leadingEdge, perl = T),glue("<b style='color: blue'>{pathway}</b>"),str_trunc(pathway,70,"right"))) %>% 
    rowwise() %>% 
    #mutate(label=if_else(any(hox %in% leadingEdge),glue("<b style='color: blue'>{pathway}</b>"),str_trunc(pathway,70,"right"))) %>% 
    ungroup() %>% 
    #ggplot(aes(NES,fct_reorder(str_wrap(pathway,25),desc(desc(NES*-log10(padj)))),fill=-log10(padj))) + 
    #ggplot(aes(NES,fct_reorder(str_trunc(pathway,70,"right"),desc(desc(NES))),fill=I(col)),text=leadingEdge) + #fill=fct_relevel(col,c("up","down")))) #desc(desc(NES*-log10(padj))
    ggplot(aes(NES,fct_reorder(str_trunc(label,70,"right"),dplyr::desc(dplyr::desc(NES))),fill=I(col))) + #fill=fct_relevel(col,c("up","down")))) #desc(desc(NES*-log10(padj))
    geom_col() +
    #geom_text(aes(label=size, hjust=if_else(NES < 0,0,1)), position = position_dodge(0.5)) +
    #scale_fill_gradient(low = "blue", high = "salmon") +
    scale_fill_manual(values = c("#5C7BE4","#D85847"), labels=c("Down","Up")) +
    coord_fixed(ratio = 0.5) + 
    theme_prism(base_size = 10, base_line_size = 0.5, border = T) +
    theme(legend.title = element_text(), 
          axis.text.y = element_markdown()) + 
    labs(y="",x="NES",title=title, fill="")
  
  return(p)
}

################
# Enrichment analysis, clusterProfiler, enrichR
#######################################################################################################################################################################################

# mydotplot<-function(df,title, showCategory,x){
#   Count<-x
#   df@result$GeneRatio<-sapply(df@result$GeneRatio,function(x){eval(parse(text=x))})
#   df %>% slot("result") %>% arrange(pvalue) %>% filter(Description!="NA") %>% filter(pvalue<0.05) %>% head(showCategory) %>% ggplot(aes(reorder(str_wrap(Description,70),Count),Count,color=pvalue,size=GeneRatio)) + 
#     geom_point() + coord_flip() + scale_color_gradient(low="red", high="blue") + 
#     theme_light()+ 
#     theme(panel.background = element_blank(), panel.grid = element_blank(), aspect.ratio = 1.5, axis.text = element_text(face="bold", family = "sans",color = "black")) +
#     #scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
#     labs(title=title) + guides(color = guide_colorbar(reverse = TRUE)) + xlab("") + ylab("Count") -> p
#   return(p)
# }


mydotplot<-function(df,title, showCategory,x){
  Count<-x
  df@result$GeneRatio<-sapply(df@result$GeneRatio,function(x){eval(parse(text=x))})
  df %>% slot("result") %>% 
    arrange(pvalue) %>% 
    filter(Description!="NA") %>% 
    filter(pvalue<0.05) %>% 
    head(showCategory) %>%
    mutate(label=str_trunc(Description,70)) %>% 
    #rowwise() %>% 
    #mutate(label=if_else(any(str_detect(geneID,c(hox))), glue("<b style='color: blue'>{label}</b>"),label)) %>% 
    #ungroup() %>%
    ggplot(aes(reorder(label,Count),Count,color=pvalue,size=GeneRatio)) + 
    geom_point() + 
    coord_flip() + 
    scale_color_gradient(low="red", high="blue") + 
    #theme_light(base_rect_size = 0.75) +
    theme_prism(base_size = 10,base_line_size = 0.5) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),aspect.ratio = 1.5, legend.title = element_text()) + #aspect.ratio = 1.5 #axis.text.y = element_markdown()
    #scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    labs(title=title) + 
    guides(color = guide_colorbar(reverse = TRUE)) + 
    xlab("") + 
    ylab("Count") -> p
  return(p)
}

mybarplot<-function(df,title, showCategory,x){
  Count<-x
  df@result$GeneRatio<-sapply(df@result$GeneRatio,function(x){eval(parse(text=x))})
  df %>% slot("result") %>% 
    arrange(pvalue) %>% 
    filter(Description!="NA") %>% 
    filter(pvalue<0.05) %>% 
    head(showCategory) %>%
    mutate(label=str_trunc(Description,50,side = "right")) %>%
    #rowwise() %>% 
    #mutate(label=if_else(any(str_detect(geneID,c(hox))), glue("<b style='color: blue'>{label}</b>"),label)) %>% 
    #ungroup() %>%
    ggplot(aes(reorder(label,-log10(pvalue)),-log10(pvalue),fill=-log10(pvalue))) + #fill=-log10(pvalue)
    #geom_col(fill="#08306B", width = 0.85) +
    geom_col(width = 0.85) +
    coord_fixed(ratio = 0.333) +
    coord_flip() +
    #scale_fill_gradient(low="navy", high="salmon") +
    scale_fill_gradientn(colours = c("navy",rev(pals:::brewer.ylorbr(15)[8:15]))) +
    geom_text(aes(label=paste("",Count)), color="white",position = position_dodge(0.9), hjust=1.25, size=3.5, show.legend=TRUE) +
    #scale_color_gradient(low="red", high="blue") + 
    #theme_light(base_rect_size = 0.75) +
    theme_prism(base_size = 10,base_line_size = 0.5) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),aspect.ratio = 0.333, legend.position = "none") +
          #legend.title = element_text()) + #aspect.ratio = 1.5 #axis.text.y = element_markdown()
    #scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    labs(x="",y=bquote(-log[10]~"p-value"),fill=bquote(-log[10]~"p-value"), title=title) -> p
  return(p)
}

myplotEnrich <- function(df,numChar, showTerms){
  df %>% 
    mutate(Count=as.numeric(str_extract(Overlap,"(.*)(?=\\/.*)"))) %>% 
    arrange(P.value) %>% 
    head(showTerms) %>% 
    #ggplot(aes(fct_reorder(str_trunc(Term,numChar,"right"),Count),Count)) + #fill=-log10(P.value)
    ggplot(aes(fct_reorder(str_trunc(Term,numChar,"right"),-log10(P.value)),-log10(P.value),fill=-log10(P.value))) + #fill=-log10(P.value)
    #geom_col(fill="#08306B", width = 0.85) +
    geom_col( width = 0.85) +
    coord_flip() +
    #scale_fill_gradient(low="navy", high="salmon") + #low="#0f0877", high="#FE9929"
    scale_fill_gradientn(colours = c("navy",rev(pals:::brewer.ylorbr(15)[8:15]))) + #pals:::brewer.ylorbr(15)[5:15]
    geom_text(aes(label=paste("",Count)), color="white",position = position_dodge(0.9), hjust=1.25, size=3.5, show.legend=TRUE) +
    #scale_color_gradient(low="red", high="blue") + 
    #theme_light(base_rect_size = 0.75) +
    theme_prism(base_size = 10,base_line_size = 0.5) +
    theme(panel.background = element_blank(), panel.grid = element_blank(),aspect.ratio = 0.3333, legend.position = "none") +
    #legend.title = element_text()) + #aspect.ratio = 1.5 #axis.text.y = element_markdown()
    #scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
    #guides(color = guide_colorbar(reverse = TRUE)) + 
    labs(x="",y=bquote(-log[10]~"p-value"),fill=bquote(-log[10]~"p-value"), title=title) -> p
  
  return(p)
}

runEnrichmentAll <- function(geneList,title,path,orgDb){
  
  #path <- paste0("~/MondalLab/NCC_RNAseq/diffExpAnalysis_nccProject/plots/compiledFigures/")
  message("Enrichment GOBP")
  goBP <- enrichGO(gene = geneList ,OrgDb = orgDb, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
  message("Enrichment GOMF")
  goMF <- enrichGO(gene = geneList ,OrgDb = orgDb, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
  message("Enrichment GOCC")
  goCC <- enrichGO(gene = geneList ,OrgDb = orgDb, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
  message("Enrichment WP")
  enrichrWP <- enrichr(geneList,"WikiPathway_2021_Human")
  message("Enrichment REACTOME")
  enrichrReact <- enrichr(geneList,"Reactome_2022")
  message("Enrichment KEGG")
  enrichrKEGG <- enrichr(geneList,"KEGG_2021_Human")
  message("Enrichment MsigDB Hallmark")
  enrichrMsigdbH <- enrichr(geneList,"MSigDB_Hallmark_2020")
  message("Enrichment Bioplanet")
  enrichrBioplanet <- enrichr(geneList,"BioPlanet_2019")
  message("Enrichment ENCODE_Histone_Modifications_2015")
  enrichrENCODE_Histone_Modifications_2015 <- enrichr(geneList,"ENCODE_Histone_Modifications_2015")
  message("Enrichment ChEA_2016")
  enrichrChEA_2016 <- enrichr(geneList,"ChEA_2016")
  message("Enrichment ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
  enrichrENCODE_ChEA_TFs <- enrichr(geneList,"ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
  message("ENCODE_TF_ChIP-seq_2015")
  enrichrENCODE_TF_ChIP <- enrichr(geneList,"ENCODE_TF_ChIP-seq_2015")
  message("Enrichment Elsevier_Pathway_Collection")
  enrichrElsevierPath <- enrichr(geneList,"Elsevier_Pathway_Collection")
  message("Enrichment CORUM")
  enrichrCORUM <- enrichr(geneList,"CORUM")
  message("Enrichment ESCAPE")
  enrichrESCAPE <- enrichr(geneList,"ESCAPE")
  
  
  
  message("Saving enrichment analysis tables")
  
  list(
    "GOBioProc"=goBP@result %>% filter(pvalue < 0.05), #%>% rowwise() %>% mutate(HOXrelated=if_else(any(str_detect(geneID,hox)),TRUE,FALSE)) %>% ungroup(),
    "GOMolFun"=goMF@result %>% filter(pvalue < 0.05) , #%>% rowwise() %>% mutate(HOXrelated=if_else(any(str_detect(geneID,hox)),TRUE,FALSE)) %>% ungroup(),
    "GOCellComp"=goCC@result %>% filter(pvalue < 0.05) ,#%>% rowwise() %>% mutate(HOXrelated=if_else(any(str_detect(geneID,hox)),TRUE,FALSE)) %>% ungroup(),
    "EnrichR_WikiPath"=enrichrWP$WikiPathway_2021_Human %>% filter(P.value < 0.05),# %>% rowwise() %>% mutate(HOXrelated=if_else(any(str_detect(Genes,hox)),TRUE,FALSE)) %>% ungroup(),
    "EnrichR_Reactome"=enrichrReact$Reactome_2022 %>% filter(P.value < 0.05), # %>% rowwise() %>% mutate(HOXrelated=if_else(any(str_detect(Genes,hox)),TRUE,FALSE)) %>% ungroup(),
    "EnrichR_KEGG"=enrichrKEGG$KEGG_2021_Human %>% filter(P.value < 0.05), # %>% rowwise() %>% mutate(HOXrelated=if_else(any(str_detect(Genes,hox)),TRUE,FALSE)) %>% ungroup(),
    "EnrichR_MsigdbH"=enrichrMsigdbH$MSigDB_Hallmark_2020 %>% filter(P.value < 0.05), #%>% rowwise() %>% mutate(HOXrelated=if_else(any(str_detect(Genes,hox)),TRUE,FALSE)) %>% ungroup(),
    "EnrichR_Bioplanet"=enrichrBioplanet$BioPlanet_2019 %>% filter(P.value < 0.05),
    "EnrichR_CORUM"=enrichrCORUM$CORUM  %>% filter(P.value < 0.05),
    "EnrichR_ChEA"= enrichrChEA_2016$ChEA_2016 %>% filter(P.value < 0.05),
    "EnrichR_ENCODE_Histone"=enrichrENCODE_Histone_Modifications_2015$ENCODE_Histone_Modifications_2015 %>% filter(P.value < 0.05),
    "EnrichR_ENCODE_ChEA_TFs"=enrichrENCODE_ChEA_TFs$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X` %>% filter(P.value < 0.05),
    "EnrichR_ENCODE_TF_ChIP"=enrichrENCODE_TF_ChIP$`ENCODE_TF_ChIP-seq_2015` %>% filter(P.value < 0.05),
    "EnrichR_ESCAPE"=enrichrESCAPE$ESCAPE %>% filter(P.value < 0.05),
    "EnrichR_ElsevPaths"=enrichrElsevierPath$Elsevier_Pathway_Collection %>% filter(P.value < 0.05)
  ) %>% openxlsx::write.xlsx(.,paste0(path,"/enrichmentAnalysis_tables_",gsub("\\s+","_",title),".xlsx"), overwrite = T)
  
                             
                             pGOBP <- mybarplot(goBP , showCategory=10, x="Count", title=paste0("GO:Biological Process \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=7), legend.title = element_text())
                             
                             pGOMF <- mybarplot(goMF , showCategory=10, x="Count", title=paste0("GO:Molecular Function \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=7), legend.title = element_text())
                             
                             pGOCC <- mybarplot(goCC , showCategory=10, x="Count", title=paste0("GO:Cellular Component \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=7), legend.title = element_text())
                             
                             
                             pWP <- myplotEnrich(enrichrWP$WikiPathway_2021_Human %>% filter(P.value<0.05),numChar = 50, showTerms = 10)  + labs(title=paste0("Enrichr: WikiPathways \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             pReactome <- myplotEnrich(enrichrReact$Reactome_2022 %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: Reactome \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             pKEGG <- myplotEnrich(enrichrKEGG$KEGG_2021_Human %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: KEGG \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=7), legend.title = element_text())
                             
                             pMsigdbH <- myplotEnrich(enrichrMsigdbH$MSigDB_Hallmark_2020 %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: Msigdb Hallmark Gene Sets \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             pBioplanet <- myplotEnrich(enrichrBioplanet$BioPlanet_2019 %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: BioPlanet \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             pENCOHist <- myplotEnrich(enrichrENCODE_Histone_Modifications_2015$ENCODE_Histone_Modifications_2015 %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: ENCODE Histone Modifications \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=7), legend.title = element_text())
                             
                             pChea <-  myplotEnrich(  enrichrChEA_2016$ChEA_2016 %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: ChEA \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             pENCOCheaTFs <-  myplotEnrich(  enrichrENCODE_ChEA_TFs$`ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X` %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: ENCODE ChEA Consensus TFs from ChIP \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             
                             pENCOTFs <-  myplotEnrich( enrichrENCODE_TF_ChIP$`ENCODE_TF_ChIP-seq_2015` %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: ENCODE TFs from ChIP-seq \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             
                             pElsePath <-  myplotEnrich( enrichrElsevierPath$Elsevier_Pathway_Collection %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: Elsevier Pathway Collection \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             
                             pCORUM <-  myplotEnrich( enrichrCORUM$CORUM %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: CORUM \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             pESCAPE <-  myplotEnrich( enrichrESCAPE$ESCAPE %>% filter(P.value<0.05),numChar = 50,showTerms = 10) + labs(title=paste0("Enrichr: ESCAPE, Embryonic Stem Cell Atlas \n",title)) + theme_prism(base_size = 10,base_line_size = 0.5) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + theme(aspect.ratio = 0.5,axis.text.y = element_text(size=8), legend.title = element_text())
                             
                             
                             pEnrich <- (pGOBP | pWP  | pReactome ) / ( pMsigdbH | pBioplanet |  pKEGG ) / (pCORUM | pESCAPE | pChea ) / (pENCOHist | pENCOCheaTFs | pENCOTFs) / (pElsePath | plot_spacer() | plot_spacer())
                               #(pGOBP | pGOMF | pGOCC  | plot_spacer() ) / (pWP  | pReactome | pKEGG) / ( pMsigdbH | pBioplanet | plot_spacer())
                             
                             message("Saving enrichment analysis plots")
                             
                             ggsave(paste0(path,"/compilation_enrichmentAnalysis_",gsub("\\s+","_",title),".pdf"), plot = pEnrich, width = 24, height = 24)
                             #ggsave(paste0(path,"/compilation_enrichmentAnalysis_",gsub("\\s+","_",title),".pdf"), plot = pEnrich, width = 24, height = 24)
                             
                             message("Saving enrichment analysis objects")
                             
                             save(goBP,goMF,goCC,enrichrWP,enrichrReact,enrichrKEGG,enrichrMsigdbH,enrichrBioplanet,enrichrESCAPE,enrichrCORUM,enrichrElsevierPath,enrichrChEA_2016,enrichrENCODE_ChEA_TFs,enrichrENCODE_Histone_Modifications_2015,enrichrENCODE_TF_ChIP,file=paste0(path,"/compilation_enrichmentAnalysis_",gsub("\\s+","_",title),".Rdata"))
                             
                             gc()
                                    
                                    message("Done.")
                                    
}

###############################
# Relative log expression
myplotRLE <- function(plotRLEres,title){

  plotRLEres %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(names_to = "sample", values_to = "rle",-gene) %>% 
    ggplot(aes(sample,rle)) + geom_boxplot(position = position_dodge2(width = 0.5)) +
    theme_prism(base_line_size = 0.5,axis_text_angle = 45) +
    labs(title=title, x="",y="Relative Log Expression")
  
}



###############################
# Plot gProfiler results barplot

myGostPlot <- function(gostResult, top,title){
    
    
    df <- gostResult$result %>% 
    filter(p_value < 0.05) %>%
    arrange(p_value) %>% 
    select(term_name, p_value, intersection_size,source) %>% 
    group_by(source) %>% 
    arrange(p_value) %>% 
    slice_head(n=top) %>% 
    ungroup() %>%
    mutate(term_name=str_trunc(term_name,70,"right")) %>%
    mutate(source=if_else(source=="WP","WikiPath",if_else(source=="REAC","REACTOME", source)))

    df$source <- factor(df$source, levels=c("GO:BP","WikiPath","REACTOME","KEGG","GO:MF","GO:CC","CORUM"))
    
    p <-  df %>%
    ggplot(aes(x=intersection_size,y=fct_reorder(term_name,intersection_size), fill=-log10(p_value), group=source)) +
    geom_col(width=0.75,position = position_dodge2(width = 0.5, preserve = "total")) +
    scale_x_continuous(expand=expansion(mult=c(0,0.05))) +
    facet_grid(rows=vars(source), scales="free",as.table = T,space="free") + 
    #facet_wrap(rows=vars(source), scales="free",as.table = F,space="free") + 
    theme_prism(base_size = 10, base_rect_size = 0.333, border = T) +
    theme(legend.title = element_text(), strip.text.y=element_text(angle=0)) +
    scale_fill_gradient(low = "navy", high = "salmon") +
    labs(title=title,y="",x="Gene Count", fill=bquote(-log[10]~"p-value"))

 return(p)
}

mygostEnrich <- function(geneList,organism, outdir, filename){
    message("Running gPRofiler enrichment analysis")
  res1 <- gost(query = geneList,organism = organism,user_threshold = 0.05,significant = TRUE ,correction_method="fdr", exclude_iea=T,evcodes = TRUE,sources = c("GO:BP","GO:MF","GO:CC","KEGG","REAC","WP","CORUM"))
  
  message("Saving xlsx results to",outdir)
  list(
    "Enrich_gProfiler"=res1$result %>% filter(p_value < 0.05) %>% arrange(desc(precision), desc(recall)) %>% select(ID=term_id,Description=term_name,Count=intersection_size,p_value, term_size, geneID=intersection, source)
    ) %>% openxlsx::write.xlsx(paste0(outdir,"/enrichmentAnalysis_gProfiler_",filename,".xlsx"), overwrite = T)
  
  message("Saving gProfiler object")
  
  save("res1", file = paste0(outdir,"/enrichmentAnalysis_gProfiler_",filename,".Rda"))
  
  return(res1)
}

read_all_sheets = function(xlsxFile, ...) {
  sheet_names = openxlsx::getSheetNames(xlsxFile)
  sheet_list = lapply(sheet_names, function(sn){openxlsx::read.xlsx(xlsxFile, sheet=sn, ...)})
  names(sheet_list) <- sheet_names
  return(sheet_list)
}

