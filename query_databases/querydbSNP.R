#Sol Ochoa
#!/usr/bin/env Rscript
#query dbSNP for a rsids list
library(reutils)
library(pbapply)

#input: lista de rsids
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Necesito un archivo con rsids que buscar", call.=FALSE)
}else if (length(args)==1){
	query=read.table(args[1])
	if(ncol(query)>1){stop("El archivo debe contener sólo una columna de rsids",call.=FALSE)}
	fetcher=efetch(query,"snp", rettype = "docset")
#fetcher=pblapply(seq(500,nrow(query),499),function(x) efetch(query[(x-499):x],"snp", rettype = "docset"))#si son muchas, hazlo por lotes
	writer=content(fetcher)
	#writer=sapply(fetcher,content)#escribe así, si lo hiciste por lote
	writeLines(writer,"dbSNP.docset")
}