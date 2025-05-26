library(readxl)
library(argparse)
library(tidyr)
library(dplyr)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-ex", "--ex", 
                    help="Path to excel file")

args = parser$parse_args()

excel_file <- args$ex

gt <- read_excel(excel_file, sheet = 1)
metadata <- read_excel(excel_file, sheet = 2)
markers <- read_excel(excel_file, sheet = 3)
loci_performance <- read_excel(excel_file, sheet = 4)
asv_table <- read_excel(excel_file, sheet = 5)
asv_seqs <- read_excel(excel_file, sheet = 6)
asv_seqs_masked <- read_excel(excel_file, sheet = 7)
vcf_like <- read_excel(excel_file, sheet = 8)
discarded_loci_gt <- read_excel(excel_file, sheet = 9)
discarded_loci_loci_performance <- read_excel(excel_file, sheet = 10)
discarded_samples_gt <- read_excel(excel_file, sheet = 11)
discarded_samples_metadata <- read_excel(excel_file, sheet = 12)
controls_gt <- read_excel(excel_file, sheet = 13)
controls_metadata <- read_excel(excel_file, sheet = 14)

cleaned <- sub("_ampseq_object_f\\.xlsx$", "", basename(excel_file))

gt_long = gt %>% mutate(experiment_sample_name = paste0(Sample_id)) %>%
  pivot_longer(
    cols = -experiment_sample_name,
    names_to = "target_id",
    values_to = "pseudocigar_masked_tmp") %>%
  filter(!is.na(pseudocigar_masked_tmp)) %>%
  mutate(pseudocigar_masked = str_remove(pseudocigar_masked_tmp, ":\\d+$"),
         reads = as.integer(str_extract(pseudocigar_masked_tmp, "(?<=:)\\d+$")),
         fwd_fastq = paste0(experiment_sample_name, "_R1_001.fastq.gz"),
         rev_fastq = paste0(experiment_sample_name, "_R2_001.fastq.gz")) %>% 
  select(-pseudocigar_masked_tmp)

gt_long_hapid <- gt_long %>%
  left_join(
    asv_table %>%
      select(Amplicon, CIGAR, hapid, CIGAR_masked) %>%
      mutate(order_index = row_number()),
    by = c("target_id" = "Amplicon", "pseudocigar_masked" = "CIGAR_masked")
  ) %>%
  arrange(order_index, experiment_sample_name) %>%
  select(-order_index) %>%
  rename(pseudocigar_unmasked = CIGAR)

gt_long_seqs <- gt_long_hapid %>% 
  left_join(
    asv_seqs,
    by = c("hapid" = "asv_id")
  ) %>% 
  rename(asv_raw = asv_seq) %>% 
  left_join(
    asv_seqs_masked,
    by = c("hapid" = "asv_id")
  ) %>% 
  rename(asv_masked = asv_seq) %>%
  mutate(pool = cleaned) %>% 
  select(
    experiment_sample_name,
    target_id,
    asv_raw,
    pseudocigar_unmasked,
    asv_masked,
    pseudocigar_masked,
    reads,
    fwd_fastq,
    rev_fastq,
    pool
  ) 

print(gt_long_seqs)

write.table(gt_long_seqs, file = "microhaplotype_info.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
