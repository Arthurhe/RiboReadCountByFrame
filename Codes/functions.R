
# human readable time (t from proc.time())
second_to_humanReadableTime=function(t){
  #change second to a vector of hour,min,second
  h=floor(t/3600)
  t=t-h*3600
  m=floor(t/60)
  t=t-m*60
  s=t
  t=c(h,m,s)
  return(t)
}

# Arthur's bedtools_intersect function #modified #Grange is 0 base
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

# Arthur's bedtools_merge function 
bedtools_merge_by_gene=function(bed,gene_chr_tbl,gene_col=5,strand_col=4,ignore_strand=F){
    suppressMessages(require(GenomicRanges))
    if(ignore_strand){
        bed <- GRanges(seqnames = bed[,gene_col],
                       ranges = IRanges(start = bed[,2],
                                        end = bed[,3]))   
    }else{
        bed <- GRanges(seqnames = bed[,gene_col],
                       ranges = IRanges(start = bed[,2],
                                        end = bed[,3]),
                       strand=bed[,strand_col])      
    }
    merged_out=data.frame(reduce(bed))
    merged_out$gene_id=merged_out$seqnames
    merged_out$chr=gene_chr_tbl$chr[match(merged_out$gene_id,gene_chr_tbl$gene_id)]
    merged_out=merged_out[,c("chr","start","end","strand","gene_id")]
    return(merged_out)
}

#return occurance of each element
# occurrance_count_4_each_element(c("A","A","B","C","C"))
# -> 1 2 1 1 2
occurrance_count_4_each_element=function(vec){
    posi_list=split(seq_along(vec), vec)
    out=rep(0,length(vec))
    for(i in 1:length(posi_list)){
        out[posi_list[[i]]]=1:length(posi_list[[i]])
    }
    return(out)
}