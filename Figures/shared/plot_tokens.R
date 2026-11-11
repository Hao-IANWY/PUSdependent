# plot_tokens.R — shared plotting constants for Figures (sizes and palettes)
# Pure data file: zero dependencies and zero side effects; just source it.
# Machinery (fonts / theme / export) lives in plot_style.R; the two files are
# independent of each other.
#
# Convention: colors use #RRGGBBCC (alpha = 0.8), except for a few values
# defined without an alpha channel.

# ============================================================================
# 1. Sizes and linewidths (the wildtype_landscape explicit font-size system)
# ============================================================================
SIZE_TOKENS <- list(
  base_linewidth = 0.6,
  base_size = 24,
  axis_text_size = 20,
  axis_title_size = 22,
  annotation_text_size = 16,
  legend_text_size = 18,
  legend_title_size = 18,
  strip_text_size = 18,
  plot_title_size = 18
)

# ============================================================================
# 2. Palettes (only palettes reused by >= 2 scripts; single-use palettes stay
#    local to the plot.R that uses them)
# ============================================================================

# Three-tissue palette, pastel variant (sites_count and ratio_distribution
# panels). Unnamed: callers name the tissues inconsistently (Root/Leaf/Flower
# vs Roots/Leaves/Flowers), so each script applies setNames() itself.
PAL_TISSUE <- c("#FF9FB2CC", "#7ED6B5CC", "#8FB8FFCC")

# Three-tissue palette, opaque NPG variant (kept as an alternative)
PAL_TISSUE_NPG <- c(
  "Root" = "#E64B35FF",
  "Leaf" = "#00A087FF",
  "Flower" = "#3C5488FF"
)

# Two-color site-status palette: non-dependent/WT (gray), dependent (red).
# Unnamed: callers label these differently ("WT sites" / "non-dependent
# sites"), so each script applies setNames() itself.
PAL_STATUS <- c("#B8B8B8CC", "#E64B35CC")

# Three-color mRNA-region palette (region_distribution panels)
PAL_REGION <- c(
  "5'UTR" = "#F6D55CCC",
  "CDS" = "#F4A7B9CC",
  "3'UTR" = "#8FD3F4CC"
)

# Two-color palette for method-comparison Venn figures; values carry no alpha
PAL_METHOD <- c(
  "methylKit" = "#E64B35",
  "dependent" = "#3C5488"
)
