library("XenofilteR", lib="/domus/h1/amendez/private/myR_packages")

out.path <- "/aln/aln_xenofilteR_mouse/" # Change to this path to filter mouse reads
#out.path <- "/aln/aln_xenofilteR_human/" # Defaul path to filter human reads.

path_human <- "/aln/grch38"
path_mouse <- "/aln/mm10"

samples_human <- list.files(path = path_human, pattern = ".bam$", full.names = TRUE)
samples_mouse <- list.files(path = path_mouse, pattern = ".bam$", full.names = TRUE)

for(i in 1:length(samples_human)){
    #cat("Filtering samples:", samples_human[i]," and ",samples_mouse[i],"\n")
    sample.list <- data.frame(samples_mouse[i], samples_human[i]) # Change to this line when filtering mouse reads
    #sample.list <- data.frame(samples_human[i], samples_mouse[i]) # Default when filtering human reads
    #sample.list <- data.frame(samples_human, samples_mouse)
    bp.param <- SnowParam(workers = 4, type = "SOCK")
    #out.name <- gsub(".bam","",samples_human[i])
    name <- gsub("(/aln/grch38/|_humanAligned.sortedByCoord.out.bam)","",samples_human[i], perl=T)
    outdir <- paste0(out.path,name)

    XenofilteR(sample.list = sample.list, destination.folder = outdir, bp.param = bp.param, output.names = NULL, MM_threshold=4, Unmapped_penalty = 8)

}



