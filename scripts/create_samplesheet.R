#!/bin/bash R

library(optparse)
option_list <- list(
  make_option(opt_str = c("-i","--in_dir"),
              default = NULL,
              help = "Input directory fastq files",
              metavar = "character"),
  make_option(opt_str = c("-o","--out_dir"),
              type = "character",
              default = NULL,
              help = "output for samplesheet",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)

args <- parse_args(opt_parser)


fastq_dir <- args$in_dir
outdir    <- args$out_dir

#files = list.files(fastq_dir)
files = list.files(fastq_dir,pattern="\\.fastq.gz$")


#samples = unique(stringr::str_remove(files,"_R1.fastq.gz|_R2.fastq.gz|.fastq.gz|_[0-9]+_S.*_L.*_R.*\\.fastq\\.gz$|_L001_R.*_001.*\\.fastq\\.gz$"))
samples = unique(stringr::str_remove(files,"_R1_001.fastq.gz|_R2_001.fastq.gz|_R1.fastq.gz|_R2.fastq.gz|.fastq.gz|_L001_R.*_001.*\\.fastq\\.gz$"))


sample    = list()
fastq_1   = list()
fastq_2   = list()
long_read = list()
for (sam in samples){
    sample[[sam]]  = sam
    fastq_1[[sam]] = files[grepl(paste0("^",sam,".*_R1.*\\.fastq\\.gz$"),files)] 
    fastq_2[[sam]] = files[grepl(paste0("^",sam,".*_R2.*\\.fastq\\.gz$"),files)]
    long_read[[sam]] = files[grepl(paste0(sam,".fastq.gz"),files)]
    # If there is no read change character(0) to "No_Read"
    fastq_1   <- lapply(fastq_1, function(x) if(identical(x,character(0))) "No_Read_R1" else x)
    fastq_2   <- lapply(fastq_2, function(x) if(identical(x,character(0))) "No_Read_R2" else x)
    long_read <- lapply(long_read, function(x) if(identical(x,character(0))) "No_Read" else x)
}

test <- lapply(long_read, function(x) if(identical(x,character(0))) "No_Read" else x)
data_dir <- fastq_dir
df <- data.frame(sample_id   = unlist(sample),
                 fastq_1   = paste0(fastq_dir,unlist(fastq_1)),
                 fastq_2   = paste0(fastq_dir,unlist(fastq_2)),
                 long_read = paste0(fastq_dir,unlist(long_read)),
                 all_fastq_files = paste0(fastq_dir,'[!Undetermined_]*.fastq.gz'),
                 data_dir = paste0(data_dir))

write.table(df,paste0(outdir,'/',"samplesheet.csv"),row.names=FALSE,sep=",",quote=FALSE)
