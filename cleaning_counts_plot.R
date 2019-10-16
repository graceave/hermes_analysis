# hermes library stats
# g avecilla
# 3 Oct 2019

# libraries before and after cleaning
# this takes the output of clean_reads.sh (run on the command line)
# specificially, NAME_cleaning_counts.txt
# it outputs a plot with number of reads per library for each library

library(tidyverse)
dir= "/Volumes/GoogleDrive/My Drive/Gresham Lab_Grace/france_satay/library_prep/france_saturation_scripts/grace_adaptations/v5_code"
setwd(dir)

# function to get NAME_cleaning_counts.txt file and parse the info
# x is the file path
get_library_info = function(x) {
  file = read_lines(x)
  df = tibble(sample=file[1], reads_original=as.numeric(file[3]), 
              reads_preinsert=as.numeric(file[5]), 
              reads_cleanf=as.numeric(file[7]), reads_cleanr=as.numeric(file[9]))
}
counts_files = paste0(dir, '/cleaning_counts/', 
                     list.files(path  = './cleaning_counts/', pattern='*cleaning_counts.txt'))
counts_data = map(counts_files, get_library_info)
counts = do.call(rbind, counts_data) %>% 
  select(-reads_cleanr) %>%
  pivot_longer(cols = starts_with("reads"), names_to = "step", values_to = "reads")

# plot
ggplot(counts, aes(sample, reads, fill = step)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 16) +
  xlab("Sample") +
  ylab("Reads") +
  scale_y_continuous(breaks = c(1e6,3e6,5e6,7e6,9e6,11e6)) +
  scale_fill_brewer(palette = "Set3", name = "Library cleaning step", 
                      labels = c("Original", "With hermes TIR", "With TIR & without plasmid"))


