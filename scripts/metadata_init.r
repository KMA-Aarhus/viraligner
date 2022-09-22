library(tidyverse)
library(readr)
library(stringr)


#' author: "Tine Sneibjerg Ebsen"

#############################################################################################################################
# This script takes the long metadata directly from pappenheim (run on workstations) and converts the table to wide format. # 
#############################################################################################################################



#############
# Read data #
#############

# Parse arguments:
args = commandArgs(trailingOnly=TRUE)
write("These are the args:", stderr())
write(paste(args), stderr())
write("", stderr())


main_batch_id = args[1] 
file_prefix = args[2]
file_out = args[3]

###############################
# Read the long metadata file #
###############################

df_wide = read_delim(paste0(file_prefix,"/",main_batch_id, "_all.tsv"), delim = "\t", col_types = cols(.default = col_character())) %>%
        rename(batch_id = `#batch_id`) %>% 
    mutate_at(vars(batch_id), as.character)  %>% 
  mutate(totalMissing=0)

for(i in 1:length(df_wide$upload_consensus_file)){
  consensus = read_file(paste0(file_prefix,"/",df_wide$upload_consensus_file[i]))
  totalMissing = str_count(consensus, pattern = "N")
  df_wide$totalMissing[i] <- totalMissing
}

# Pivot to wide format
df_wide = df_wide %>%
  select(-starts_with("unnamed:_")) %>%
    mutate(totalMissing_interpreted = case_when(is.na(totalMissing) ~ 197233,
                                                   TRUE ~ as.numeric(totalMissing)),
           sample_name_prefix = str_sub(sample_id, 1, 2), # first two characters
           sample_name_suffix = str_sub(sample_id, 3), # the rest
           full_name = paste0(batch_id, "_", sample_id),
           sample_name_prefix_converted = recode(sample_name_prefix,
                                                 "87" = paste0("L"),
                                                 "88" = paste0("R"),
                                                 "89" = paste0("I"),
                                                 "90" = paste0("V"),
                                                 "96" = paste0("P"))) %>% 
    
    rowwise() %>% 
    mutate(ya_sample_id = paste0(sample_name_prefix_converted, sample_name_suffix)) %>% 
    ungroup() %>% 
    
    
    select(-sample_name_prefix, -sample_name_suffix, -sample_name_prefix_converted) %>% 
    select(full_name, batch_id, sample_id, ya_sample_id, type, everything()) 




# Read barcode_alignment.tsv which contains the number of reads for each sample
barcode_alignment = read_delim(paste0(file_prefix,"/",main_batch_id,"_barcode_alignment.tsv"), delim = "\t") %>% 
    select(barcode_base = barcode, ba_type = type, ba_reads = target_unclassified, started) %>% 
    mutate(barcode= str_replace(barcode_base,"barcode","NB"))

sum_reads = barcode_alignment %>% pull(ba_reads) %>% sum
unclassified_reads = barcode_alignment %>% filter(ba_type == "na") %>% pull(ba_reads)

# Join barcode_alignment onto df_wide, and add information about the number of reads.
df_wide_ba = df_wide %>% 
    left_join(barcode_alignment) 

df_wide_ba = df_wide_ba %>% 
    mutate(batch_sum_reads = sum_reads,
           batch_unclassified_reads = unclassified_reads, 
           batch_unclassified_reads_prop = unclassified_reads/batch_sum_reads)

# Read final_summary.txt, which contains information about time.
terminal_summary = read_delim(paste0(file_prefix,"/",main_batch_id, "_final_summary.txt"), delim = "=", col_names = c("name", "value")) %>% 
    pivot_wider(names_from = name, values_from = value) %>% 

    mutate_at(vars(started, acquisition_stopped, processing_stopped), lubridate::ymd_hms) %>% 
    # Calculate run time
    #Vi
    identity




# Incorporate time information
# Parse hours sequencing
start_time = terminal_summary %>% pull(started)
end_time = terminal_summary %>% pull(acquisition_stopped) %>% head(1)
#sequencing_time = end_time-start_time
hours_sequencing = lubridate::interval(start_time, end_time) / lubridate::hours(1)

# add to the main df
df_wide_ba_hours= df_wide_ba %>% 
    mutate(batch_hours_sequencing = hours_sequencing,
           batch_million_reads_per_day = (batch_sum_reads/1000000)/(batch_hours_sequencing/24))



##################################
# Check that the controls are OK #
##################################

# Make a summary that can be used for checking the presence of controls.
type_summary = df_wide_ba %>% group_by(very_long_batch_id, type) %>%
    summarize(n = length(type), min_totalMissing = min(totalMissing_interpreted), max_totalMissing = max(totalMissing_interpreted)) # If there is more than one positive or negative control, we need to know the min and max-values of them.


error_messages = c() # Empty vector to collect error_messages.

# TODO: Test this thoroughly!

# Positive control presence
if(type_summary %>% filter(type == "positive_control") %>% nrow < 1) { 
    error_messages = c(error_messages, "No positive control present.")
} else if (type_summary %>% filter(type == "positive_control") %>% pull(min_totalMissing) >= 3000) { # Positive control coverage
    error_messages = c(error_messages, "The positive control totalMissing is above the critical threshold.") 
}


# Negative control presence
if(type_summary %>% filter(type == "negative_control") %>% nrow < 1) { 
    error_messages = c(error_messages, "No negative control present.")
} else if (type_summary %>% filter(type == "negative_control") %>% pull(max_totalMissing) <= 29903/2) { # Negative control coverage
    error_messages = c(error_messages, "The negative control totalMissing is below the critical threshold.")
}


if(length(error_messages) > 0) {
    error_flag = TRUE
    write("During processing of the incoming samples, the following control problems were encountered:", stderr())
    write(paste("Error:", error_messages), stderr())
} else {
    error_flag = FALSE
    error_messages = "QC OK"
    write("Info: During processing of the incoming samples, no control problems were found :)", stderr())
}

df_wide_ba_hours_stamp = df_wide_ba_hours %>% 
    mutate(batch_control_stamp = case_when(error_flag ~ "failed",
                                           TRUE ~ "passed"),
           batch_control_error_messages = paste(error_messages, collapse = ", "))


# Write the init file out.
df_wide_ba_hours_stamp %>% write_tsv(file_out)

