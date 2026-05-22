# merge_vf_abundance.R
library(tidyverse)

# 获取所有样本的 vf_abundance.tsv 文件
files <- list.files("D_b2_VF_abundance", pattern = "vf_abundance.tsv", full.names = TRUE)
samples <- gsub(".*/(.*)\\.vf_abundance\\.tsv", "\\1", files)

# 读取并合并
df_list <- map2(files, samples, ~ {
  dat <- read.delim(.x, header = FALSE, stringsAsFactors = FALSE)
  colnames(dat) <- c("gene_id", "coverage", "rpkm", "tpm", "cpm")
  dat$sample <- .y
  # 读取对应注释文件，获取 VF ID
  annot_file <- paste0("D_b1_VFDB_annotation/", .y, ".vfdb.tsv")
  annot <- read.delim(annot_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(annot) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                       "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")
  annot <- annot %>% select(qseqid, sseqid) %>% distinct()
  left_join(dat, annot, by = c("gene_id" = "qseqid"))
})

merged <- bind_rows(df_list) %>%
  select(sample, sseqid, tpm, cpm, rpkm) %>%
  filter(!is.na(sseqid))

# 转为宽格式（以 VF ID 为行，样本为列）
tpm_wide <- merged %>% select(sample, sseqid, tpm) %>%
  pivot_wider(names_from = sample, values_from = tpm, values_fill = 0)

cpm_wide <- merged %>% select(sample, sseqid, cpm) %>%
  pivot_wider(names_from = sample, values_from = cpm, values_fill = 0)

fpkm_wide <- merged %>% select(sample, sseqid, rpkm) %>%
  pivot_wider(names_from = sample, values_from = rpkm, values_fill = 0)

# 保存
write.table(tpm_wide, "VF_TPM_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cpm_wide, "VF_CPM_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fpkm_wide, "VF_FPKM_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)