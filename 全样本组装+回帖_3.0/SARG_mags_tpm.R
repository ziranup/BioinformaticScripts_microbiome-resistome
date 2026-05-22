# build_MAG_ARG_TPM.R
library(dplyr)
library(tidyr)
library(stringr)

# 设置路径
diamond_dir <- "./D_b7_DIAMOND_SARG"
output_tpm <- "./MAG_ARG_TPM.tsv"

# 获取所有 DIAMOND 结果文件
files <- list.files(diamond_dir, pattern = "_vs_SARG\\.tsv$", full.names = TRUE)

# 提前获取每个 MAG 的总 ORF 数（从 *_proteins.faa 行数 / 2）
orf_dir <- "./D_b6_ORFs_MAGs"
orf_files <- list.files(orf_dir, pattern = "_proteins\\.faa$", full.names = TRUE)
orf_counts <- sapply(orf_files, function(f) {
  n_lines <- length(readLines(f, warn = FALSE))
  n_orfs <- n_lines / 2  # 每个蛋白占两行（header + seq）
  mag_id <- str_remove(basename(f), "_proteins.faa")
  return(c(mag_id, n_orfs))
})
orf_df <- data.frame(
  MAG = orf_counts[1, ],
  total_ORFs = as.numeric(orf_counts[2, ]),
  stringsAsFactors = FALSE
)

# 读取并合并所有 DIAMOND 结果
all_hits <- NULL
for (f in files) {
  if (file.size(f) == 0) next
  df <- read.delim(f, header = FALSE, stringsAsFactors = FALSE, quote = "")
  colnames(df) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                    "qlen", "slen", "stitle")
  mag_id <- str_remove(basename(f), "_vs_SARG.tsv")
  df$MAG <- mag_id
  all_hits <- bind_rows(all_hits, df)
}

# 若无任何命中
if (is.null(all_hits)) {
  stop("没有检测到任何 ARG 命中！")
}

# 过滤：identity >= 80%
filtered_hits <- all_hits %>%
  filter(pident >= 80)

# 提取 ARG ID（SARG 数据库中 sseqid 格式通常为 "ARG|Class|Subclass|..."）
# 我们直接使用 sseqid 作为 ARG identifier
arg_table <- filtered_hits %>%
  count(MAG, sseqid, .drop = FALSE) %>%
  rename(ARG = sseqid, hits = n)

# 合并总 ORF 数
arg_table <- left_join(arg_table, orf_df, by = "MAG")

# 计算 TPM: (hits / total_ORFs) * 1e6
arg_table <- arg_table %>%
  mutate(TPM = (hits / total_ORFs) * 1e6)

# 构建宽表（MAG 为行，ARG 为列）
tpm_wide <- arg_table %>%
  select(MAG, ARG, TPM) %>%
  pivot_wider(names_from = ARG, values_from = TPM, values_fill = 0)

# 保存
write.table(tpm_wide, output_tpm, sep = "\t", quote = FALSE, row.names = FALSE)

cat("MAG-ARG TPM 丰度表已保存至:", output_tpm, "\n")