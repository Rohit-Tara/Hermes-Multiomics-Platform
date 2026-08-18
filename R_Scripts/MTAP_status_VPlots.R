library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

# run the Wilcoxon tests directly and draw the significance brackets manually
# using geom_segment and annotate, so no extra significance package is needed

# use separate colours for MTAP wild-type and deleted groups
# so they stay visually distinct from the BAP1 and machine-learning plots
mtap_colours <- c('wt' = 'slateblue', 'del' = 'goldenrod')

df <- read_csv('MEDUSA_Master_ULTIMATE.csv', guess_max = 3000)

df <- df %>%
  mutate(
    MTAP_status = case_when(
      `DriverLoss 9p21 (CDKN2A,MTAP,IFNA)` == 0 ~ 'wt',
      `DriverLoss 9p21 (CDKN2A,MTAP,IFNA)` == 1 ~ 'del'
    ),
    MTAP_status = factor(MTAP_status, levels = c('wt', 'del'))
  )

# convert p-values into the same significance labels used in the Python analysis
p_to_stars <- function(p) {
  if (is.na(p)) return('n.s.')
  if (p < 0.001) return('***')
  if (p < 0.01) return('**')
  if (p < 0.05) return('*')
  return('n.s.')
}

# create a violin plot for a single feature, used here for activated dendritic cells
violin_panel <- function(data, feature_col, y_label = 'Cell Type Abundance') {
  
  plot_df <- data %>%
    select(MTAP_status, value = all_of(feature_col)) %>%
    filter(!is.na(MTAP_status), !is.na(value))
  
  n_labels <- plot_df %>%
    count(MTAP_status) %>%
    mutate(label = paste0('n=', n))
  
  # run the Wilcoxon test and position the significance bracket above the data
  wt_vals <- plot_df$value[plot_df$MTAP_status == 'wt']
  del_vals <- plot_df$value[plot_df$MTAP_status == 'del']
  p_val <- if (length(wt_vals) >= 3 && length(del_vals) >= 3) {
    wilcox.test(wt_vals, del_vals)$p.value
  } else NA
  sig_label <- p_to_stars(p_val)
  
  y_max <- max(plot_df$value, na.rm = TRUE)
  y_range <- diff(range(plot_df$value, na.rm = TRUE))
  bracket_y <- y_max + y_range * 0.08
  tick <- y_range * 0.02
  
  ggplot(plot_df, aes(MTAP_status, value, fill = MTAP_status)) +
    geom_violin(colour = 'black', linewidth = 0.4, alpha = 0.85, trim = TRUE) +
    geom_jitter(width = 0.08, size = 0.6, alpha = 0.35, colour = 'black') +
    stat_summary(fun = median, geom = 'crossbar', width = 0.4, colour = 'black', linewidth = 0.3) +
    # draw the significance bracket using two vertical ticks and a horizontal line
    annotate('segment', x = 1, xend = 1, y = bracket_y - tick, yend = bracket_y) +
    annotate('segment', x = 2, xend = 2, y = bracket_y - tick, yend = bracket_y) +
    annotate('segment', x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
    annotate('text', x = 1.5, y = bracket_y + tick, label = sig_label, size = 3.2) +
    geom_text(data = n_labels, aes(MTAP_status, -Inf, label = label),
              inherit.aes = FALSE, vjust = -0.6, size = 2.6) +
    scale_fill_manual(values = mtap_colours, guide = 'none') +
    labs(title = feature_col, x = 'MTAP Status', y = y_label) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(size = 8, face = 'bold', hjust = 0.5))
}

# apply the same approach across all features from one immune deconvolution method
immune_grid <- function(data, prefix, ncol = 4) {
  
  feat_cols <- names(data)[startsWith(names(data), prefix)]
  
  long_df <- data %>%
    select(MTAP_status, all_of(feat_cols)) %>%
    pivot_longer(-MTAP_status, names_to = 'feature', values_to = 'value') %>%
    filter(!is.na(MTAP_status), !is.na(value)) %>%
    mutate(feature = str_remove(feature, prefix))
  
  n_labels <- long_df %>%
    count(feature, MTAP_status) %>%
    mutate(label = paste0('n=', n))
  
  # run a separate Wilcoxon test for each feature
  # build the significance information for each facet without an extra package
  # and add a bracket and significance label to each feature
  sig_df <- long_df %>%
    group_by(feature) %>%
    summarise(
      p_val = if (sum(MTAP_status == 'wt') >= 3 & sum(MTAP_status == 'del') >= 3) {
        wilcox.test(value[MTAP_status == 'wt'], value[MTAP_status == 'del'])$p.value
      } else NA_real_,
      y_max = max(value, na.rm = TRUE),
      y_range = diff(range(value, na.rm = TRUE)),
      .groups = 'drop'
    ) %>%
    mutate(
      sig_label = sapply(p_val, p_to_stars),
      bracket_y = y_max + y_range * 0.08,
      tick = y_range * 0.02
    )
  
  ggplot(long_df, aes(MTAP_status, value, fill = MTAP_status)) +
    geom_violin(colour = 'black', linewidth = 0.4, alpha = 0.85, trim = TRUE) +
    geom_jitter(width = 0.08, size = 0.5, alpha = 0.3, colour = 'black') +
    stat_summary(fun = median, geom = 'crossbar', width = 0.4, colour = 'black', linewidth = 0.3) +
    geom_segment(data = sig_df, aes(x = 1, xend = 1, y = bracket_y - tick, yend = bracket_y),
                 inherit.aes = FALSE) +
    geom_segment(data = sig_df, aes(x = 2, xend = 2, y = bracket_y - tick, yend = bracket_y),
                 inherit.aes = FALSE) +
    geom_segment(data = sig_df, aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
                 inherit.aes = FALSE) +
    geom_text(data = sig_df, aes(x = 1.5, y = bracket_y + tick, label = sig_label),
              inherit.aes = FALSE, size = 3) +
    facet_wrap(~ feature, scales = 'free_y', ncol = ncol) +  # allow each cell type to use its own y-axis scale
    scale_fill_manual(values = mtap_colours, guide = 'none') +
    labs(x = 'MTAP Status', y = 'Cell Type Abundance') +
    theme_classic(base_size = 9) +
    theme(strip.text = element_text(size = 8, face = 'bold'),
          strip.background = element_rect(fill = 'grey92', colour = NA))
}

# use the same plotting approach for a selected list of features
# the HRD/TF panel contains both ScarHRD scores and decoupleR transcription-factor features
# so the columns are supplied explicitly rather than selected using a shared prefix
custom_grid <- function(data, col_names, display_names = col_names, ncol = 5) {
  
  long_df <- data %>%
    select(MTAP_status, all_of(col_names)) %>%
    rename(setNames(col_names, display_names)) %>%
    pivot_longer(-MTAP_status, names_to = 'feature', values_to = 'value') %>%
    filter(!is.na(MTAP_status), !is.na(value)) %>%
    mutate(feature = factor(feature, levels = display_names))  # keep the panels in the specified order rather than alphabetical order
  
  sig_df <- long_df %>%
    group_by(feature) %>%
    summarise(
      p_val = if (sum(MTAP_status == 'wt') >= 3 & sum(MTAP_status == 'del') >= 3) {
        wilcox.test(value[MTAP_status == 'wt'], value[MTAP_status == 'del'])$p.value
      } else NA_real_,
      y_max = max(value, na.rm = TRUE),
      y_range = diff(range(value, na.rm = TRUE)),
      .groups = 'drop'
    ) %>%
    mutate(
      sig_label = sapply(p_val, p_to_stars),
      bracket_y = y_max + y_range * 0.08,
      tick = y_range * 0.02
    )
  
  ggplot(long_df, aes(MTAP_status, value, fill = MTAP_status)) +
    geom_violin(colour = 'black', linewidth = 0.4, alpha = 0.85, trim = TRUE) +
    geom_jitter(width = 0.08, size = 0.5, alpha = 0.3, colour = 'black') +
    stat_summary(fun = median, geom = 'crossbar', width = 0.4, colour = 'black', linewidth = 0.3) +
    geom_segment(data = sig_df, aes(x = 1, xend = 1, y = bracket_y - tick, yend = bracket_y),
                 inherit.aes = FALSE) +
    geom_segment(data = sig_df, aes(x = 2, xend = 2, y = bracket_y - tick, yend = bracket_y),
                 inherit.aes = FALSE) +
    geom_segment(data = sig_df, aes(x = 1, xend = 2, y = bracket_y, yend = bracket_y),
                 inherit.aes = FALSE) +
    geom_text(data = sig_df, aes(x = 1.5, y = bracket_y + tick, label = sig_label),
              inherit.aes = FALSE, size = 3) +
    facet_wrap(~ feature, scales = 'free_y', ncol = ncol) +
    scale_fill_manual(values = mtap_colours, guide = 'none') +
    labs(x = 'MTAP Status', y = 'Score') +
    theme_classic(base_size = 9) +
    theme(strip.text = element_text(size = 8, face = 'bold'),
          strip.background = element_rect(fill = 'grey92', colour = NA))
}

# generate the ConsensusTME violin plots
p1 <- immune_grid(df, 'ConsensusTME:')
ggsave('figure_immune_violin_consensustme_R.png', p1, width = 16, height = 20, dpi = 300, limitsize = FALSE)

# generate the CIBERSORT violin plots
p2 <- immune_grid(df, 'Cibersort:')
ggsave('figure_immune_violin_cibersort_R.png', p2, width = 16, height = 24, dpi = 300, limitsize = FALSE)

# generate the xCell violin plots
p3 <- immune_grid(df, 'xCell:')
ggsave('figure_immune_violin_xcell_R.png', p3, width = 16, height = 40, dpi = 300, limitsize = FALSE)

# plot the activated dendritic-cell result separately because it remains significant after FDR correction
p4 <- violin_panel(df, 'Cibersort:Dendritic cells activated',
                   y_label = 'Cibersort activated dendritic cell score') +
  labs(title = 'Activated dendritic cell abundance by MTAP status\n(FDR = 0.021)') +
  theme(plot.title = element_text(size = 11, face = 'bold'))
ggsave('figure_dendritic_cells_mtap_violin_R.png', p4, width = 5, height = 5, dpi = 300)

# generate violin plots for the same HRD and transcription-factor features used in the corresponding comparison
hrd_tf_cols <- c('ScarHRD LOH score', 'ScarHRD LST score', 'ScarHRD sum score',
                 'Clonal TMB', 'Loss ratio',
                 'decoupleR:ARID1A', 'decoupleR:HCFC1', 'decoupleR:PAX5',
                 'decoupleR:SOX5', 'decoupleR:STOX1')
hrd_tf_labels <- c('HRD LOH score', 'HRD LST score', 'HRD sum score',
                   'Clonal TMB', 'Loss ratio',
                   'TF: ARID1A', 'TF: HCFC1', 'TF: PAX5', 'TF: SOX5', 'TF: STOX1')

p5 <- custom_grid(df, hrd_tf_cols, hrd_tf_labels, ncol = 5)
ggsave('figure_hrd_tf_violin_grid_R.png', p5, width = 20, height = 8, dpi = 300, limitsize = FALSE)

cat('done - 5 figures saved\n')

