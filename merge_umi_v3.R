# hermes library merging umis within 700 bp windows and within 1 bp of each other
# g avecilla
# 3 Oct 2019

# this code takes bed files (output of bam2bed_v5.sh)
# Rscript from command line
# pass like so
# r/intel/3.4.2
# Rscript merge_umi.R /path/NAME_sorted.bed chromosomenumber /path/NAME_chr_umi.bed
# this version takes only one chromosome at a time, which will hopefully speed things u

#clear global env
rm(list=ls())


#get arguments from command line
args = commandArgs(TRUE)

#set variables
bed_files = args[1]
chr = args[2]
output = args[3]

library(tidyverse)
library(stringdist)

# function to get NAME_sorted.bed file and make the umi into a separate column
# x is the file path
# chromosome "17" is mitochondria
read_bed = function(x, c = chr) {
  data = read_tsv(x, col_names = F) %>%
    mutate(chromosome = case_when(X1== 'NC_001133.9' ~ 1,
                                  X1== 'NC_001134.8' ~ 2,
                                  X1== 'NC_001135.5' ~ 3,
                                  X1== 'NC_001136.10' ~ 4,
                                  X1== 'NC_001137.3' ~ 5,
                                  X1== 'NC_001138.5' ~ 6,
                                  X1== 'NC_001139.9' ~ 7,
                                  X1== 'NC_001140.6' ~ 8,
                                  X1== 'NC_001141.2' ~ 9,
                                  X1== 'NC_001142.9' ~ 10,
                                  X1== 'NC_001143.9' ~ 11,
                                  X1== 'NC_001144.5' ~ 12,
                                  X1== 'NC_001145.3' ~ 13,
                                  X1== 'NC_001146.8' ~ 14,
                                  X1== 'NC_001147.6' ~ 15,
                                  X1== 'NC_001148.4' ~ 16,
                                  X1 == 'NC_001224.1' ~ 17)) %>%
    rename(start = X2, end = X3, name = X4, mapq = X5, strand = X6) %>%
    separate(name, into = c("rname", "t"), sep = "_") %>%
    separate(t, into = "umi", sep = "/") %>%
    select(chromosome, start, end, rname, mapq, strand, umi) %>%
    filter(chromosome == c)
}

# function to merge umis within 700 bp windows
# allows for difference of 1 base within umi
# x is the bed data

roll_window = function(x){
  for(i in unique(x$chromosome)) {
    x_sub = x %>% filter(chromosome == i)
    min = min(x_sub$start)
    max = max(x_sub$end)
    for(j in min:max) {
      x_sub1 = x_sub %>% filter(start >= j & start <= j+700)
      if(length(unique(x_sub1$umi)) < nrow(x_sub1)) { #if there are fewer umis then there are reads, some reads must share umis, so pass to merge umi
        x_sub = x_sub %>% 
          filter(start < j | start > j+700) #%>%
          bind_rows(merge_umis(x_sub1))
      } else if(any(stringdistmatrix(x_sub1$umi, method = "hamming") == 1)) { #here we check if even though there is a "unique" umi for each read, none of them should be merged
        x_sub = x_sub %>% 
          filter(start < j | start > j+700) %>%
          bind_rows(merge_umis(x_sub1))
      } else { #if neither of the above conditions are met, bind the data as is
        x_sub = x_sub
      }
    }
    x = x %>% filter(chromosome != i) %>% bind_rows(x_sub)
  }
  return(x)
}

roll_window = function(x_sub){
  min = min(x_sub$start)
  max = max(x_sub$end)
  for(j in min:max) {
    x_sub1 = x_sub %>% filter(start >= j & start <= j+700)
    if(length(unique(x_sub1$umi)) < nrow(x_sub1)) { #if there are fewer umis then there are reads, some reads must share umis, so pass to merge umi
      x_sub = x_sub %>% 
        filter(start < j | start > j+700) %>%
        bind_rows(merge_umis(x_sub1))
    } else if(any(stringdistmatrix(x_sub1$umi, method = "hamming") == 1)) { #here we check if even though there is a "unique" umi for each read, none of them should be merged
      x_sub = x_sub %>% 
        filter(start < j | start > j+700) %>%
        bind_rows(merge_umis(x_sub1))
    } else { #if neither of the above conditions are met, bind the data as is
        x_sub = x_sub
      }
    }
  return(x_sub)
}


# this function chooses a representative read (the one with the highest quality)
# for reads with the same umi
merge_umis = function(data) {
  merged = NULL
  #first, merge change umis with hamming distance of 1 to have the same sequence
  if(any(stringdistmatrix(data$umi, method = "hamming") == 1)) {
    x = close_distance(data)
  } else {
    x = data
  }
  for(i in 1:length(unique(x$umi))) {
    if(length(which(x$umi == unique(x$umi)[i])) > 1) { #merge identical umis
      x_sub = x %>% filter(umi == unique(x$umi)[i])
      merged = bind_rows(merged, x_sub[which(x_sub$mapq == max(x_sub$mapq))[1],])
    } else { #keep umis only corresponding to one read
      merged = bind_rows(merged, x[which(x$umi == unique(x$umi)[i]),])
    }
  }
  return(merged)
}


# this function finds umis that have a hamming distance from one to each other
# and changes them all to be the same (based on which is more frequent)

###close distance needs help####
close_distance = function(x) {
  while(any(stringdistmatrix(x$umi, method = "hamming") == 1)) {
    close = x$umi[which(as.matrix(stringdistmatrix(x$umi, method = "hamming")) == 1, arr.ind = T)[1,1]]
    inds1 = which(stringdist(close, x$umi, method = "hamming") == 1)
    inds_ori = which(x$umi == close)
    if(length(inds_ori) > length(inds1)) {
      x$umi[inds1] = close
    } else {
      x$umi[inds_ori] = x$umi[inds1[1]]
    }
  }
  return(x)
}

# get data
bed_data = read_bed(bed_files, chr)

# merge umis
merged_bed = roll_window(bed_data)

# write output
write_tsv(merged_bed %>% arrange(start), file.path(output))



