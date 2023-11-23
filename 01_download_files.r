pkgTest <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, ask = F)
    }
    require(x)
  }
}

pkgTest("glue")
pkgTest("R.utils")
pkgTest("curl")

target_info <- read.csv("links.csv")
head(target_info)


############ Download fastq files ##########################

download.fastq.files <- function(row_info){
  dir.create("raw_fastq_files", showWarnings = F)
  
  sample_files <- c(glue("raw_fastq_files/{row_info$Sample}_R1.fastq"),glue("raw_fastq_files/{row_info$Sample}_R2.fastq"))
  
  # download R1
  dest_file <- sample_files[1]
  curl_download(row_info$R1, destfile = dest_file, quiet = FALSE)
  
  gunzip(filename = dest_file, destname = glue("{dest_file}.gz"))
  
  # download R2
  dest_file <- sample_files[2]
  curl_download(row_info$R2, destfile = dest_file, quiet = FALSE)
  
  gunzip(filename = dest_file, destname = glue("{dest_file}.gz"))
  
  sample_files <- c(glue("raw_fastq_files/{row_info$Sample}_R1.fastq.gz"),glue("raw_fastq_files/{row_info$Sample}_R2.fastq.gz"))
  
  return(sample_files)
}

sample_files <- lapply(1:nrow(target_info), function(i){
  row_info <- target_info[i,]
  sample_files <- download.fastq.files(row_info)
})

names(sample_files) <- target_info$Sample

dir.create("rds_data", showWarnings = F)
dir.create("rds_data/raw_fastq", showWarnings = F)

saveRDS(sample_files, file = "rds_data/raw_fastq/sample_files.rds")
