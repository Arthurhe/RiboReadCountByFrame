#updated support for mm9 and mm10, some code optimization by Arthur Jan 27th, 2020

#start timer
message("------------------")
cat("setting things up ... ")
ptm <- proc.time()

#check package availability 
if(!all(c("GenomicRanges","dplyr","stringr","optparse","data.table") %in% installed.packages())){stop("not all required package installed!\nRequired libraries: GenomicRanges, dplyr, stringr, optparse, data.table")}

suppressMessages(require(GenomicRanges))
suppressMessages(require(dplyr))
suppressMessages(require(stringr))
suppressMessages(require(optparse))
suppressMessages(require(data.table))

# Find path of running Rscript
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

# load functions
source(paste0(script.basename,"/Codes/functions.R"))

# Use optparse package to pass argument to Rscrip from command line (with flag)
option_list = list(make_option(c("-f", "--file"), type="character", default=NULL, 
                               help="read sam file name", metavar="character"), 
                   make_option(c("-g", "--gtf"), type="character", default=NULL, 
                               help="gtf file name", metavar="character"),
                   make_option(c("-o", "--out"), type="character", default=NULL, 
                               help="output file prefix, will output two files: prefix.tsv & prefix.gtf", metavar="character"),
                   make_option(c("-r", "--reference"), type="character", default="hg38", 
                               help="hg38, mm10 or mm9 [default= %default]", metavar="character"))
                   
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if(!opt$reference %in% c("hg38","mm9","mm10")){stop("-r/--reference must be hg38, mm9 or mm10")}
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2)) 
cat(paste0(t[1],"h-",t[2],"m-",round(t[3]),"s\n"))


# ---------------------------
# Step 0: Run on terminal
cat("processing input gtf ... ")

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

gtf = fread(opt$gtf, sep="\t", stringsAsFactors=F,header=F,data.table=F)

# Remove duplicated rows, if input is from multiple gtf file concatenated from multiple samples
gtf = distinct(gtf)
gtf=gtf[gtf$V3 == "CDS",]

# Select the columns we need and filter only "CDS" 
gtf_bed = gtf[,c(1, 4, 5, 8, 9, 7)]
colnames(gtf_bed) = c("chr", "str", "end", "frame", "info", "strand_col")

# extract orf_id
gtf_bed$orf_id = gsub("\"; transcript_id","",str_extract(gtf_bed$info,"(ENS[MUS]*T[0-9]+.*:chr.+; transcript_id)"))

#reorder
gtf_bed=gtf_bed[,c("chr", "str", "end", "frame", "orf_id", "strand_col")]

# change to appropiate column data type
gtf_bed$frame = as.numeric(gtf_bed$frame)

gtf_bed = arrange(gtf_bed, orf_id)

# -------------------------
# Step 1.2: Create a merged ORF table for annotation and a dataframe for intersect

# return list of positions for each element
index = split(seq_along(gtf_bed$orf_id), gtf_bed$orf_id)

mORF_vec = rep("",length(index))
for (i in 1:length(index)){
    
    idx_vec = index[[i]]
    
    str = sapply(as.character(gtf_bed$str[idx_vec]), paste0, ":")
    end = sapply(as.character(gtf_bed$end[idx_vec]), paste0, ":")
    frame = sapply(as.character(gtf_bed$frame[idx_vec]), paste0, "|")
    
    chr = paste0(gtf_bed$chr[idx_vec[1]], ";")
    strand = gtf_bed$strand[idx_vec[1]]
    
    mORF = paste(c(rbind(str, end, frame)), collapse="")
    mORF = paste0(mORF, chr, strand)
    
    # name is first old orf_id, value is string of CDSs
    mORF_vec[i] = mORF
    names(mORF_vec)[i]=gtf_bed$orf_id[idx_vec[1]]
    }

t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2)) 
cat(paste0(t[1],"h-",t[2],"m-",round(t[3]),"s\n"))

# ------------------------
# Step 1.3: Create Annotation table for merged_ORF to geneID

# Read in gene annotation file
cat("comparing gtf with reference annotation ... ")
anno_path = paste0(script.basename, "/Codes/",opt$reference,"_t2g_R35.Rds")
anno = readRDS(anno_path)

# Create mORF_table by adding column "geneID" "tx_id" "merged_orfID" > for annotation
tx_id_vec = str_extract(gtf_bed$orf_id, "ENS[MUS]*T[0-9]+")

idx = match(tx_id_vec, anno$target_id)

mORF_table = gtf_bed %>%
            mutate(tx_id = tx_id_vec, ens_gene = anno$ens_gene[idx], ext_gene = anno$ext_gene[idx], 
                   merged_orfID = as.character(mORF_vec[gtf_bed$orf_id]))


# Use distinct() to merge same merged_ORF (CDSs string)
mORF_table = mORF_table %>%
    select(tx_id, ens_gene, ext_gene, merged_orfID) %>%  # Note: merged_ORF can match many tx_id, here just show one
    distinct(merged_orfID, .keep_all = TRUE) %>%
    arrange(ext_gene)

# create a vector of number to add to geneID, eg. A2M.1, A2M.2, A2M.3, AAAS.1 ...
gene_nvec=occurrance_count_4_each_element(mORF_table$ext_gene)


#identify NAs
ext_na_idx=is.na(mORF_table$ext_gene)
ens_na_idx=is.na(mORF_table$ens_gene)

# Add number to geneID
mORF_table = mORF_table %>%
    mutate(ens_gene = paste0(ens_gene, ".", gene_nvec), ext_gene = paste0(ext_gene, ".", gene_nvec))

# Replace output NA.$number with NA
mORF_table$ens_gene[ens_na_idx] = "NA"
mORF_table$ext_gene[ext_na_idx] = "NA"


# ---------------------
# Step 1.4: Create a merged ORF dataframe for intersect

# Transform string into data.frame: chr, str, end, frame, merge_orfID, strand 

# ---------------------
# Step 1.4: Create a merged ORF dataframe for intersect

# Transform string into data.frame: chr, str, end, frame, merge_orfID, strand 

mORF_unique = mORF_vec[!duplicated(mORF_vec)]

CDS_blocks=gtf_bed
CDS_blocks=CDS_blocks[CDS_blocks$orf_id %in% names(mORF_unique),]
CDS_blocks$full_id=mORF_vec[match(CDS_blocks$orf_id,names(mORF_vec))]
CDS_blocks$simp_id=mORF_table$ext_gene[match(CDS_blocks$full_id,mORF_table$merged_orfID)]
#CDS_blocks$simp_id_count=occurrance_count_4_each_element(CDS_blocks$simp_id)
#CDS_blocks$exon_id=paste0(CDS_blocks$simp_id,".",CDS_blocks$simp_id_count)
CDS_blocks$gene_id=gsub("\\.[0-9]+$","",CDS_blocks$simp_id)

mORF_df=CDS_blocks[,c("chr","str","end","frame","simp_id","strand_col")]

t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2)) 
cat(paste0(t[1],"h-",t[2],"m-",round(t[3]),"s\n"))

# -----------------------------
# Step 2: Read in "FC2_Index5_*offsetCorrect_post_bedtools_intersect.bed"
cat("counting read numbers in ORF... ")

# Read in reads bed file 
reads = fread(post_bed_name, sep="\t", stringsAsFactors=F,data.table=F,header=F)
colnames(reads) = c("chr", "start", "end", "strand_col")
system(paste0("rm ",post_bed_name))

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

t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2)) 
cat(paste0(t[1],"h-",t[2],"m-",round(t[3]),"s\n"))

# ---------------------------
cat("writing output ... ")
# ---------------------------
#output gtf
gtf_CDS=data.frame(V1=CDS_blocks$chr,
                   V2="ArthurRiboPipeline",
                   V3="CDS",
                   V4=CDS_blocks$str,
                   V5=CDS_blocks$end,
                   V6=".",
                   V7=CDS_blocks$strand_col,
                   V8=".",
                   V9=paste0("gene_id \"",CDS_blocks$gene_id,"\"; transcript_id \"",CDS_blocks$simp_id,"\";"),
                   stringsAsFactors=F
                   )

fwrite(gtf_CDS,file = paste0(opt$out,".gtf"),quote = F,sep="\t",row.names = F,col.names = F)

# ---------------------------
#output tsv
fwrite(result_f1_f2_f3, paste0(opt$out,".tsv"), sep="\t", quote = FALSE, row.names= FALSE)
t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2)) 
cat(paste0(t[1],"h-",t[2],"m-",round(t[3]),"s\n"))

# Number of in-frame reads 
message(paste0("total reads: ", nrow(result)))
message(paste0("in-frame reads: ", nrow(result_inframe), " (", round(nrow(result_inframe)/nrow(result),3)*100), "% of total)")

