#!/bin/bash
# Create a bed file from offset_correct sam file, remove header and contain four columns: seqname, start position, end position,strand
# Tony 12/09/2019

FILENAME=$1
OUT_FILENAME=$2


# Remove sam file header
grep -v "^@" $FILENAME | \

# Substitute "16" and "0" to strand "-" and "+", respectively
# Select column "chr", "start" and create 3rd column "start+1" 
awk '{gsub("16","-",$2); gsub("0","+",$2); print $3 "\t" $4-3 "\t" $4+3 "\t" $2}' > $OUT_FILENAME  #expand the reads to make sure inclusion -3 so the frame doesn't change


