#!/usr/bin/env Rscript
# merge_flower_low.R — 合并 root_vs_flower / leaf_vs_flower 两次 beta-binomial 比较的
# is_dependent 位点（并集 = Flower-low），补 Flower 深度与 Sprinzl 坐标，按代表性核
# tRNA 家族体系拆成入图 / 不入图两张表。输入输出与口径见 README.md。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, "../../../.."))
fisher_dir <- file.path(repo_root, "tRNA/STAT/fisher/ncRNA")
family_top_path <- file.path(repo_root, "tRNA/count/ncRNA/family_top/family_top_tRNA.tsv")

# ---------------------------------------------------------------------------
# 1) 两次比较的 line-level 表 → 各自的 is_dependent 位点，full_join 取并集
# ---------------------------------------------------------------------------
read_dependent <- function(ref_tissue) {
  path <- file.path(script_dir, "betabinomial", paste0(tolower(ref_tissue), "_vs_flower"),
                    paste0("Flower_vs_", ref_tissue, ".dependent-sites.tsv"))
  if (!file.exists(path)) stop("Missing bb line-level output: ", path, "\n先跑 data.sh。")

  readr::read_tsv(path, show_col_types = FALSE) |>
    filter(is_dependent == TRUE) |>
    transmute(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      "{ref_tissue}_Ratio" := as.numeric(WT_Ratio),
      "{ref_tissue}_Delta" := as.numeric(`Delta_Ratio(WT-Mut)`),
      "{ref_tissue}_p" := as.numeric(p.value),
      "{ref_tissue}_q" := as.numeric(q.value),
      Flower_Ratio_from_this_call = as.numeric(Mut_Ratio)
    )
}

root_dep <- read_dependent("Root")
leaf_dep <- read_dependent("Leaf")
if (nrow(root_dep) == 0 && nrow(leaf_dep) == 0) {
  stop("No Flower-low dependent sites in either comparison; nothing to merge.")
}

# 并集：只在一个比较候选池里的位点，另一侧留 NA。两侧都有时 Flower_Ratio 一致
# （同为 Flower 三重复 Treat ratio 均值），coalesce 收成一列。
merged <- full_join(root_dep, leaf_dep, by = c("chr", "pos", "strand")) |>
  mutate(
    Flower_Ratio = coalesce(Flower_Ratio_from_this_call.x, Flower_Ratio_from_this_call.y),
    comparison = case_when(
      !is.na(Root_q) & !is.na(Leaf_q) ~ "Root_vs_Flower,Leaf_vs_Flower",
      !is.na(Root_q) ~ "Root_vs_Flower",
      TRUE ~ "Leaf_vs_Flower"
    )
  ) |>
  select(-Flower_Ratio_from_this_call.x, -Flower_Ratio_from_this_call.y)

# ---------------------------------------------------------------------------
# 2) 补 Flower 深度（加权用）与 Sprinzl 坐标
# ---------------------------------------------------------------------------
flower_fisher <- readr::read_tsv(file.path(fisher_dir, "WT_Flower_sites.Fisher.tsv"),
                                 show_col_types = FALSE) |>
  transmute(chr = as.character(chr), pos = as.integer(pos),
            strand = as.character(strand),
            Flower_Depth = as.numeric(Treat_avg_depth))

# Sprinzl_pos 列在不同组织表里可能被 readr 猜成 character（有空值）或 numeric
# （全数值），统一走字符层再转整数，避免 if_else 类型冲突。
parse_sprinzl <- function(x) {
  chr_x <- as.character(x)
  chr_x[chr_x %in% c("", "NA")] <- NA_character_
  suppressWarnings(as.integer(chr_x))
}

sprinzl_union <- purrr::map_dfr(c("Root", "Leaf", "Flower"), function(t) {
  readr::read_tsv(file.path(fisher_dir, paste0("WT_", t, "_sites.positive.Sprinzl.tsv")),
                  show_col_types = FALSE) |>
    transmute(chr = as.character(chr), pos = as.integer(pos),
              strand = as.character(strand),
              Sprinzl_pos = parse_sprinzl(Sprinzl_pos))
}) |>
  distinct(chr, pos, strand, Sprinzl_pos) |>
  group_by(chr, pos, strand) |>
  summarise(Sprinzl_pos = first(na.omit(Sprinzl_pos)), .groups = "drop")

merged <- merged |>
  left_join(flower_fisher, by = c("chr", "pos", "strand")) |>
  left_join(sprinzl_union, by = c("chr", "pos", "strand")) |>
  mutate(Flower_Depth = tidyr::replace_na(Flower_Depth, 0))

# ---------------------------------------------------------------------------
# 3) 代表性核 tRNA 家族体系：入图 / 不入图
# ---------------------------------------------------------------------------
family_top <- readr::read_tsv(family_top_path, show_col_types = FALSE) |>
  filter(str_starts(feature_id, "^tRNA-Nu-")) |>
  transmute(Rep_ID = family, Representative_ID = feature_id) |>
  distinct(Rep_ID, Representative_ID)
if (nrow(family_top) == 0) stop("family_top 表里没有核 tRNA 家族: ", family_top_path)

classified <- merged |>
  mutate(fam = paste0("Nu-", str_remove(str_remove(chr, "^tRNA-Nu-"), "-[0-9]+$"))) |>
  left_join(family_top, by = c("fam" = "Rep_ID")) |>
  mutate(
    exclusion_reason = case_when(
      !str_starts(chr, "^tRNA-Nu-") ~ "non-nuclear tRNA",
      is.na(Representative_ID) ~ "family not in family_top",
      is.na(Sprinzl_pos) | Sprinzl_pos < 0 | Sprinzl_pos > 76 ~ "missing/invalid Sprinzl_pos",
      TRUE ~ NA_character_
    ),
    in_heatmap = is.na(exclusion_reason)
  )

sites_out <- classified |>
  transmute(
    chr, pos, strand, comparison,
    Root_Ratio, Leaf_Ratio, Flower_Ratio,
    Root_Delta, Leaf_Delta, Root_p, Root_q, Leaf_p, Leaf_q,
    Flower_Depth, Sprinzl_pos,
    Rep_ID = fam, Representative_ID, in_heatmap, exclusion_reason
  ) |>
  arrange(desc(in_heatmap), Rep_ID, pos)

readr::write_tsv(sites_out, file.path(script_dir, "flower_low_sites.tsv"))

mapped_out <- classified |>
  filter(in_heatmap) |>
  transmute(
    comparison, Rep_ID = fam, Representative_ID, chr, pos, strand,
    ratio = Flower_Ratio, depth = Flower_Depth, canonical_pos = Sprinzl_pos
  ) |>
  arrange(Rep_ID, canonical_pos, pos)
readr::write_tsv(mapped_out, file.path(script_dir, "flower_low_mapped_sites.tsv"))

excluded_out <- classified |>
  filter(!in_heatmap) |>
  transmute(chr, pos, strand, comparison, exclusion_reason) |>
  arrange(exclusion_reason, chr, pos)
readr::write_tsv(excluded_out, file.path(script_dir, "flower_low_excluded_sites.tsv"))

message(sprintf(
  "Flower-low sites: %d total | heatmap %d (root-only %d, leaf-only %d, both %d) | excluded %d",
  nrow(sites_out), nrow(mapped_out),
  sum(str_detect(sites_out$comparison, "^Root_vs_Flower$")),
  sum(str_detect(sites_out$comparison, "^Leaf_vs_Flower$")),
  sum(str_detect(sites_out$comparison, ",")),
  nrow(excluded_out)
))
message("Wrote flower_low_sites.tsv / flower_low_mapped_sites.tsv / flower_low_excluded_sites.tsv")
