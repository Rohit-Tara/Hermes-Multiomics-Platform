library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(survival)  # used for Cox models, Kaplan-Meier curves and log-rank tests

df <- read_csv('/home/rt334/Downloads/MEDUSA_Master_ULTIMATE.csv', guess_max = 3000)

clocks <- c('ageAcc.Horvath', 'ageAcc.Hannum', 'ageAcc.Levine', 'ageAcc.Hovarth2',
            'ageAcc.PedBE', 'ageAcc.Wu', 'ageAcc.TL', 'ageAcc.BLUP', 'ageAcc.EN')
clock_display <- c(
  'ageAcc.Horvath' = 'Horvath', 'ageAcc.Hannum' = 'Hannum', 'ageAcc.Levine' = 'Levine (PhenoAge)',
  'ageAcc.Hovarth2' = 'Horvath (skin & blood)', 'ageAcc.PedBE' = 'PedBE', 'ageAcc.Wu' = 'Wu',
  'ageAcc.TL' = 'Telomere Length (TL)', 'ageAcc.BLUP' = 'BLUP', 'ageAcc.EN' = 'Elastic Net (EN)'
)

genomic_measures <- c('ScarHRD sum score', 'Clonal TMB', 'Subclonal TMB')

os_col <- 'Overall survival (Surgery-to-Death) in days'
os_status_col <- 'Status (0=Censor 1=Dead)'
pfs_col <- 'Progression-free survival (Surgery-to-Progression) in days'
pfs_status_col <- 'Progression (0=censor 1=progression)'
bap1_col <- 'Overall BAP1 IHC (1=positive/retained, 0=negative/loss)'
mtap_col <- 'DriverLoss 9p21 (CDKN2A,MTAP,IFNA)'

# convert the variables used in the analyses to numeric values
num_cols <- c(clocks, genomic_measures, os_col, os_status_col, pfs_col, pfs_status_col, bap1_col, mtap_col)
df[num_cols] <- lapply(df[num_cols], as.numeric)


# compare methylation age acceleration with genomic instability
# each clock is tested against HRD, clonal TMB and subclonal TMB
results <- list()
for (c in clocks) {
  for (g in genomic_measures) {
    sub <- df %>% select(all_of(c(c, g))) %>% drop_na()
    if (nrow(sub) < 10) next
    test <- cor.test(sub[[c]], sub[[g]], method = 'spearman')
    results[[length(results) + 1]] <- data.frame(
      `Epigenetic clock` = clock_display[[c]], `Genomic measure` = g,
      n = nrow(sub), `Spearman rho` = round(unname(test$estimate), 3),
      `p-value` = round(test$p.value, 4), check.names = FALSE
    )
  }
}

integ_df <- bind_rows(results)
integ_df$FDR <- round(p.adjust(integ_df$`p-value`, method = 'BH'), 4)
integ_df <- integ_df %>% arrange(`p-value`)

write_csv(integ_df, 'table_methylation_vs_genomic_instability.csv')
print(integ_df)


# plot the Horvath skin and blood clock against HRD
# this was the strongest association identified in the comparisons above
plot_col <- 'ageAcc.Hovarth2'
plot_df <- df %>% select(x = `ScarHRD sum score`, y = all_of(plot_col)) %>% drop_na()
r <- cor.test(plot_df$y, plot_df$x, method = 'spearman')

p1 <- ggplot(plot_df, aes(x, y)) +
  geom_point(colour = 'steelblue4', alpha = 0.75, size = 2.5) +
  geom_smooth(method = 'lm', se = FALSE, colour = 'darkorange3', linewidth = 1) +  # the fitted line is for visualisation; the statistical test is Spearman
  labs(x = 'ScarHRD sum score', y = 'Horvath (skin & blood) age acceleration',
       title = sprintf('Methylation age acceleration vs genomic instability\n(Spearman rho=%.2f, p=%.3f, n=%d)',
                       unname(r$estimate), r$p.value, nrow(plot_df))) +
  theme_classic(base_size = 11)

ggsave('figure_methylation_vs_hrd_scatter.png', p1, width = 7, height = 6, dpi = 300)
cat('saved figure_methylation_vs_hrd_scatter.png\n')


# examine survival after stratifying patients by BAP1 and MTAP status
# this tests whether the clock-survival relationship differs within these genomic subgroups
# a separate median split is calculated within each subgroup before the survival tests
stratified_survival <- function(subgroup_col, subgroup_labels, duration_col, event_col, clock_col) {
  
  out <- list()
  for (level in names(subgroup_labels)) {
    label <- subgroup_labels[[level]]
    sub <- df %>%
      filter(.data[[subgroup_col]] == as.numeric(level)) %>%
      select(clock = all_of(clock_col), duration = all_of(duration_col), event = all_of(event_col)) %>%
      drop_na()
    
    if (nrow(sub) < 10) {
      out[[label]] <- data.frame(Subgroup = label, n = nrow(sub), note = 'too few patients to test')
      next
    }
    
    med <- median(sub$clock)
    sub$group <- ifelse(sub$clock >= med, 'High', 'Low')
    
    lr <- survdiff(Surv(duration, event) ~ group, data = sub)
    lr_p <- 1 - pchisq(lr$chisq, df = 1)  # convert the survdiff chi-square statistic into a p-value
    
    cox <- coxph(Surv(duration, event) ~ clock, data = sub)
    cox_summary <- summary(cox)
    
    out[[label]] <- data.frame(
      Subgroup = label, n = nrow(sub),
      `n (High)` = sum(sub$group == 'High'), `n (Low)` = sum(sub$group == 'Low'),
      `Log-rank p` = round(lr_p, 4),
      `Cox HR` = round(cox_summary$coefficients[1, 'exp(coef)'], 3),
      `Cox p` = round(cox_summary$coefficients[1, 'Pr(>|z|)'], 4),
      check.names = FALSE
    )
  }
  
  bind_rows(out)
}

bap1_labels <- c(`1` = 'BAP1 retained', `0` = 'BAP1 loss')
mtap_labels <- c(`0` = 'MTAP intact', `1` = 'MTAP loss')

os_by_bap1 <- stratified_survival(bap1_col, bap1_labels, os_col, os_status_col, plot_col)
os_by_mtap <- stratified_survival(mtap_col, mtap_labels, os_col, os_status_col, plot_col)
pfs_by_bap1 <- stratified_survival(bap1_col, bap1_labels, pfs_col, pfs_status_col, plot_col)
pfs_by_mtap <- stratified_survival(mtap_col, mtap_labels, pfs_col, pfs_status_col, plot_col)

strat_table <- bind_rows(
  os_by_bap1 %>% mutate(Outcome = 'OS', Split = 'BAP1'),
  os_by_mtap %>% mutate(Outcome = 'OS', Split = 'MTAP'),
  pfs_by_bap1 %>% mutate(Outcome = 'PFS', Split = 'BAP1'),
  pfs_by_mtap %>% mutate(Outcome = 'PFS', Split = 'MTAP')
) %>%
  select(Outcome, Split, Subgroup, n, `n (High)`, `n (Low)`, `Log-rank p`, `Cox HR`, `Cox p`)

write_csv(strat_table, 'table_survival_stratified_by_subgroup.csv')
print(strat_table)


# generate Kaplan-Meier curves for the BAP1 and MTAP subgroups
plot_km_by_subgroup <- function(subgroup_col, subgroup_labels, duration_col, event_col,
                                clock_col, title_text, filename) {
  
  # analyse each subgroup separately so that each has its own median split
  # and its own log-rank result displayed on the figure
  all_curves <- list()
  
  for (level in names(subgroup_labels)) {
    label <- subgroup_labels[[level]]
    sub <- df %>%
      filter(.data[[subgroup_col]] == as.numeric(level)) %>%
      select(clock = all_of(clock_col), duration = all_of(duration_col), event = all_of(event_col)) %>%
      drop_na()
    
    if (nrow(sub) < 10) next  # skip subgroups with fewer than 10 complete cases
    
    med <- median(sub$clock)
    sub$group <- ifelse(sub$clock >= med, 'High', 'Low')
    
    lr <- survdiff(Surv(duration, event) ~ group, data = sub)
    lr_p <- 1 - pchisq(lr$chisq, df = 1)
    panel_title <- sprintf('%s\n(log-rank p=%.3f)', label, lr_p)
    
    fit <- survfit(Surv(duration, event) ~ group, data = sub)
    curve_df <- data.frame(
      time = fit$time, surv = fit$surv, lower = fit$lower, upper = fit$upper,
      group = rep(names(fit$strata), fit$strata)
    )
    
    curve_df$group <- sub('group=', '', curve_df$group)  # remove the model prefix before displaying group names in the legend
    curve_df$panel <- panel_title
    all_curves[[label]] <- curve_df
  }
  
  km_df <- bind_rows(all_curves)
  
  p <- ggplot(km_df, aes(time, surv, colour = group, fill = group)) +
    geom_step(linewidth = 0.8) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    facet_wrap(~ panel, scales = 'free_x') +
    scale_colour_manual(values = c('High' = 'mediumpurple', 'Low' = 'seagreen')) +
    scale_fill_manual(values = c('High' = 'mediumpurple', 'Low' = 'seagreen')) +
    labs(x = 'Days', y = 'Survival probability', title = title_text, colour = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(strip.text = element_text(face = 'bold'))
  
  ggsave(filename, p, width = 6 * length(all_curves), height = 5.5, dpi = 300, limitsize = FALSE)
  cat('saved', filename, '\n')
}

plot_km_by_subgroup(bap1_col, bap1_labels, os_col, os_status_col, plot_col,
                    'Overall survival by Horvath (skin & blood) age acceleration,\nstratified by BAP1 status',
                    'figure_km_os_by_bap1.png')
plot_km_by_subgroup(mtap_col, mtap_labels, os_col, os_status_col, plot_col,
                    'Overall survival by Horvath (skin & blood) age acceleration,\nstratified by MTAP status',
                    'figure_km_os_by_mtap.png')

cat('\ndone - tables and figures saved for section 4.7\n')


# render the results tables as PNG files using the same style as the other tables
render_table_png <- function(tbl, title_text, filename) {
  
  n_row <- nrow(tbl)
  n_col <- ncol(tbl)
  
  tbl_chr <- as.data.frame(lapply(tbl, as.character), stringsAsFactors = FALSE)
  
  # set the width of each column using its longest entry so longer labels remain readable
  col_chars <- sapply(seq_len(n_col), function(j) {
    max(nchar(c(names(tbl)[j], tbl_chr[[j]])))
  })
  col_width <- pmax(col_chars * 0.11, 0.9)  # keep a minimum width for shorter columns
  col_left <- c(0, cumsum(col_width)[-n_col])
  col_centre <- col_left + col_width / 2
  total_width <- sum(col_width)
  
  cell_df <- data.frame()
  for (i in seq_len(n_row)) {
    for (j in seq_len(n_col)) {
      cell_df <- rbind(cell_df, data.frame(
        row = i, y = n_row - i + 1, x = col_centre[j], w = col_width[j],
        label = tbl_chr[i, j], is_header = FALSE
      ))
    }
  }
  
  header_df <- data.frame(row = 0, y = n_row + 1, x = col_centre, w = col_width,
                          label = names(tbl), is_header = TRUE)
  all_df <- rbind(cell_df, header_df)
  
  all_df$fill <- ifelse(all_df$is_header, 'steelblue',
                        ifelse(all_df$row %% 2 == 0, 'gray95', 'white'))
  all_df$text_colour <- ifelse(all_df$is_header, 'white', 'black')
  all_df$font_face <- ifelse(all_df$is_header, 'bold', 'plain')
  
  p <- ggplot(all_df, aes(x, y, width = w, fill = fill)) +
    geom_tile(colour = 'white', linewidth = 1.2, height = 1) +
    geom_text(aes(label = label, colour = text_colour, fontface = font_face), size = 3.3) +
    scale_fill_identity() +
    scale_colour_identity() +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    labs(title = title_text) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = 'bold', hjust = 0, size = 11, margin = margin(b = 10)),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  # leave additional space when the title is wider than the table itself
  title_lines <- strsplit(title_text, '\n')[[1]]
  title_width_est <- max(nchar(title_lines)) * 0.095
  canvas_width <- max(total_width * 0.55 + 1.5, title_width_est + 1)
  
  ggsave(filename, p, width = canvas_width, height = 0.35 * (n_row + 1) + 1, dpi = 300, limitsize = FALSE)
  cat('saved', filename, '\n')
}

render_table_png(
  head(integ_df, 10),
  'Table . Top 10 of 27 methylation vs genomic instability correlations\n(Spearman, ranked by p-value; FDR-corrected across all 27 clock x measure pairs)',
  'table_methylation_vs_hrd.png'
)

render_table_png(
  strat_table,
  'Table. Overall/progression-free survival by Horvath (skin & blood) age\nacceleration, stratified by BAP1 and MTAP subgroup (median split within subgroup)',
  'table_survival_stratified.png'
)

