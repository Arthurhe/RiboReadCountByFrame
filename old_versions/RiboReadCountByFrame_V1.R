#!/usr/bin/env Rscript

# R.version 3.6.1

# Install packages if not already done
# BiocManager::install("GenomicRanges")
# install.packages("dplyr")
# install.packages("stringr")
# install.packages("optparse")

# Load packages
library(GenomicRanges)
library(dplyr)
library(stringr)
library(optparse)


# Find path of running Rscript
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

message(script.basename)


# Use optparse package to pass argument to Rscrip from command line (with flag)
option_list = list(make_option(c("-f", "--file"), type="character", default=NULL, 
                               help="read sam file name", metavar="character"), 
                   make_option(c("-g", "--gtf"), type="character", default=NULL, 
                               help="gtf file name", metavar="character"),
                   make_option(c("-o", "--out"), type="character", default="out.txt", 
                               help="output file name [default= %default]", metavar="character"))
                   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)



# ---------------------------
# Step 0: Run on terminal

# Create read bed file from sam file
sh_name = paste0(script.basename, "/Codes/sam_to_bed_strand_v2.sh")
read_basename = sub(".sam", "", basename(opt$file))
pre_bed_name = paste0(script.basename, "/Temp/", read_basename, "_pre_bedtools_intersect.tsv")

system(paste(sh_name, opt$file, pre_bed_name))

# Create gtf bed file and filter exon for pre bedtools intersect in linux
repre_gtf_bed = paste0(script.basename, "/Temp/", read_basename, "_repre_gtf2bed.bed")

# need a lot of \ to escape quotation mark 
# uniq again, becasue orfID can have many same exon coordinate > falsely increase read count
system(paste("grep \'exon\\s\'", opt$gtf, "| awk \'{print $1 \"\t\" $4 \"\t\" $5 \"\t\" $7}\' | uniq >", repre_gtf_bed))
 

# Run bedtools intersect, -s: match strand, -wa: if A intersect B, keep A
post_bed_name = paste0(script.basename, "/Temp/", read_basename, "_post_bedtools_intersect.tsv")

system(paste("bedtools intersect -a", pre_bed_name, "-b", repre_gtf_bed, "-wa>", post_bed_name))



# -----------------------------
# Step 1: Read in "repre.valid.ORF.*gtf" and prepare for intersect

gtf = read.table(opt$gtf, sep="\t", stringsAsFactors=FALSE)

# Remove duplicated rows, if input is from multiple gtf file concatenated from multiple samples
gtf = distinct(gtf)

# Select the columns we need and filter only "CDS" 
gtf_selected = gtf %>%
    select(V1, V3, V4, V5, V8, V9, V7) %>%
    filter(V3 == "CDS")


# extract orf_id
gtf_selected$V9 = str_match(gtf_selected$V9,"\\s(ENST[0-9]+.*:chr.+);\\stranscript_id")[,2]

colnames(gtf_selected) = c("chr", "CDS", "str", "end", "frame", "orf_id", "strand_col")

# remove CDS column
gtf_bed = gtf_selected[,-2]

# change to appropiate column data type
gtf_bed$frame = as.numeric(gtf_bed$frame)
gtf_bed$orf_id = as.character(gtf_bed$orf_id)

gtf_bed = arrange(gtf_bed, orf_id)



# -------------------------
# Step 1.2: Create a merged ORF table for annotation and a dataframe for intersect

# find identical elements in a vector, return a list of positions(vector), need to sort the vector first
find_identical_element = function(x){
    
    listv = list()
    
    # first index
    first_idx_v = which(!duplicated(x))
    last_idx_v = c((first_idx_v - 1)[-1], length(x))
    
    for (i in 1:length(first_idx_v)){
        
        listv[i] = list(first_idx_v[i]:last_idx_v[i])
    }
    return(listv)
}


# Create a string from a set of CDSs (also frame, strand) from each orfID, then make into a list
index = find_identical_element(gtf_bed$orf_id)

mORF_list = list()

for (i in 1:length(index)){
    
    idx_vec = unlist(index[i])
    
    str = sapply(as.character(gtf_bed$str[idx_vec]), paste0, ":")
    end = sapply(as.character(gtf_bed$end[idx_vec]), paste0, ":")
    frame = sapply(as.character(gtf_bed$frame[idx_vec]), paste0, "|")
    
    chr = paste0(unique(gtf_bed$chr[idx_vec]), ";")
    strand = unique(gtf_bed$strand[idx_vec])
    
    mORF = paste(c(rbind(str, end, frame)), collapse="")
    mORF = paste0(mORF, chr, strand)
    
    # make a list: name is orf_id, value is string of CDSs
    orf_id = gtf_bed$orf_id[idx_vec[1]]
    
    mORF_list[orf_id] = list(mORF)
    }



# Read in gene annotation file
anno_path = paste0(script.basename, "/Codes/hg38_t2g_R35.Rds")
anno = readRDS(anno_path)



# ------------------------
# Step 1.3: Create Annotation table for merged_ORF to geneID

# Create mORF_table by adding column "geneID" "tx_id" "merged_orfID" > for annotation
tx_id_vec = str_extract(gtf_bed$orf_id, "ENST[0-9]+")

idx = match(tx_id_vec, anno$target_id)

mORF_table = gtf_bed %>%
            mutate(tx_id = tx_id_vec, ens_gene = anno$ens_gene[idx], ext_gene = anno$ext_gene[idx], 
                   merged_orfID = as.character(mORF_list[gtf_bed$orf_id]))


# Use distinct() to merge same merged_ORF (CDSs string)
mORF_table = mORF_table %>%
    select(tx_id, ens_gene, ext_gene, merged_orfID) %>%  # Note: merged_ORF can match many tx_id, here just show one
    distinct(merged_orfID, .keep_all = TRUE) %>%
    arrange(ext_gene)


# create a vector of number to add to geneID, eg. A2M.1, A2M.2, A2M.3, AAAS.1 ...
gene_count = mORF_table %>% count(ext_gene)

gene_nvec = c()
for (i in 1:nrow(gene_count)){ 
    
    # create a vector for each geneID, starting from 1
    v = 1:gene_count$n[i]
    
    gene_nvec = append(gene_nvec, v)
}


# Add number to geneID
mORF_table = mORF_table %>%
    mutate(ens_gene = paste0(ens_gene, ".", gene_nvec), ext_gene = paste0(ext_gene, ".", gene_nvec))

# Replace output NA.$number with NA
NA_idx = str_which(mORF_table$ext_gene, "NA.*")
mORF_table$ens_gene[NA_idx] = NA


# ---------------------
# Step 1.4: Create a merged ORF dataframe for intersect

# Transform string into data.frame: chr, str, end, frame, merge_orfID, strand 

mORF_unique = unique(mORF_list)

mORF_df = data.frame()

# Make a df for each mORF, then join the df together
for (x in mORF_unique){
    
    str_end = unlist(strsplit(x, "\\|"))
    n = length(str_end)
    
    pos = unlist(lapply(str_end[-n], strsplit, ":"))

    mat = matrix(pos, ncol=3, byrow=TRUE)
    
    chr_strand = unlist(strsplit(str_end[n], ";"))
    
    nr = nrow(mat)
    df = data.frame(chr = rep(chr_strand[1], nr), str = mat[,1], end = mat[,2], frame = mat[,3],
                    merge_orfID = rep(x, nr), strand_col = rep(chr_strand[2], nr), stringsAsFactors=FALSE)

    df$str = as.numeric(df$str)
    df$end = as.numeric(df$end)
    df$frame = as.numeric(df$frame)
    
    mORF_df = bind_rows(mORF_df, df)
    
}


# -----------------------------
# Step 2: Read in "FC2_Index5_*offsetCorrect_post_bedtools_intersect.bed"

# Read in reads bed file 
reads = read.table(post_bed_name, sep="\t", stringsAsFactors=FALSE)
colnames(reads) = c("chr", "start", "end", "strand_col")

head(reads, 20)

# Number of reads
nrow(reads)


# Arthur's bedtools_intersect function #modified

bedtools_intersect=function(bed1,bed2,overlap=T,count=F,maxgap=0,minoverlap=0,ignore_strand=T){
    suppressMessages(require(GenomicRanges))
    if(ignore_strand){
        bed1 <- GRanges(seqnames = bed1[,1],
                        ranges = IRanges(start = bed1[,2],
                                         end = bed1[,3]))
        bed2 <- GRanges(seqnames = bed2[,1],
                        ranges = IRanges(start = bed2[,2],
                                         end = bed2[,3]))        
    }else{
        bed1 <- GRanges(seqnames = bed1[,1],
                        ranges = IRanges(start = bed1[,2],
                                         end = bed1[,3]),
                        strand=bed1[,'strand_col'])
        bed2 <- GRanges(seqnames = bed2[,1],
                        ranges = IRanges(start = bed2[,2],
                                         end = bed2[,3]),
                        strand=bed2[,'strand_col'])          
    }
    if(overlap){
        overlap=findOverlaps(bed1,bed2,ignore.strand = ignore_strand,maxgap=maxgap,minoverlap=minoverlap)
        overlap=data.frame(overlap)
        colnames(overlap)=c("bed1_idx","bed2_idx")
    }
    if(count){
        count=countOverlaps(bed1,bed2,ignore.strand = ignore_strand,maxgap=maxgap,minoverlap=minoverlap)
    }
    return(list(overlap=overlap,count=count))
}


# -----------------------
# Step 3: Intersect reads and mergeORF
# query (reads) is the first argument, subejct(ORF) the second

list = bedtools_intersect(reads, mORF_df, ignore_strand=FALSE)

overlap = list$overlap


# Use dplyr package, do all steps by piping: 
    # Add read_position, cds_start, cds_frame by index; 
    # Add remainder column by calculating (read_pos - cds_start - cds frame)/3 and output remainder
    # Add tx_id column by index

result = overlap %>%
    mutate(read_pos = reads[overlap$bed1_idx, 2], cds_start = mORF_df[overlap$bed2_idx, 2], 
           cds_frame = mORF_df[overlap$bed2_idx, 4], remainder = (read_pos - cds_start - cds_frame) %% 3, 
           orfID = mORF_df[overlap$bed2_idx, 5])

# Select result that are in frame
result_inframe = result %>%
    filter(remainder == 0)


# Number of in-frame reads 
paste0(" in_frame reads: ", nrow(result_inframe))
paste0(" total reads: ", nrow(result))


result_f1 = result %>% 
    filter(remainder==0) %>%
    count(orfID)

colnames(result_f1) = c("merged_orfID", "f1Num")
    
result_f2 = result %>% 
    filter(remainder==1) %>%
    count(orfID)

colnames(result_f2) = c("merged_orfID", "f2Num")
    
result_f3 = result %>% 
    filter(remainder==2) %>%
    count(orfID)

colnames(result_f3) = c("merged_orfID", "f3Num")

# join the three tables together
result_f1_f2_f3 = full_join(full_join(result_f1, result_f2, by='merged_orfID'), result_f3, by='merged_orfID')

# Substitue NA to 0
result_f1_f2_f3[is.na(result_f1_f2_f3)] = 0

# Order by f1Num
result_f1_f2_f3 = arrange(result_f1_f2_f3, desc(f1Num))

# Add gene IDs
idx1 = match(result_f1_f2_f3$merged_orfID, mORF_table$merged_orfID)

result_f1_f2_f3 = result_f1_f2_f3 %>%
                mutate(ens_gene = mORF_table$ens_gene[idx1], ext_gene = mORF_table$ext_gene[idx1])

#Re-order column, easier to see f1Num
result_f1_f2_f3 = result_f1_f2_f3[,c(2,3,4,5,6,1)]

head(result_f1_f2_f3,20)

write.table(result_f1_f2_f3, opt$out, sep="\t", quote = FALSE, row.names= FALSE)









