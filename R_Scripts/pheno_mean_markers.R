# load the CONFIRM immune marker data
data <- read.csv(
  "~/Downloads/CONFIRM pheno_percent suppressor.csv",
  header = TRUE
)

# keep only the immune marker columns used for this analysis
markers <- data[, c("CD68", "CD163", "PDL1", "CD86", "VISTA", "CD16", "AF")]

# convert the selected columns to numeric values
markers <- as.data.frame(lapply(markers, as.numeric))

# calculate the mean abundance for each marker
marker_means <- colMeans(markers, na.rm = TRUE)
print(marker_means)

# calculate the summary statistics that sit behind each boxplot
# this includes the median, quartiles, range, sample size and mean
# so the numerical values can be checked alongside the visual output
marker_summary <- data.frame(
  Marker = names(markers),
  n = sapply(markers, function(x) sum(!is.na(x))),
  Mean = sapply(markers, mean, na.rm = TRUE),
  Median = sapply(markers, median, na.rm = TRUE),
  Q1 = sapply(markers, function(x) quantile(x, 0.25, na.rm = TRUE)),
  Q3 = sapply(markers, function(x) quantile(x, 0.75, na.rm = TRUE)),
  IQR = sapply(markers, IQR, na.rm = TRUE),
  Min = sapply(markers, min, na.rm = TRUE),
  Max = sapply(markers, max, na.rm = TRUE)
)

# round a copy of the summary table for display while keeping the original values unchanged
marker_summary_display <- marker_summary
marker_summary_display[, -1] <- round(marker_summary_display[, -1], 4)

print(marker_summary_display)
write.csv(marker_summary, "table_immune_marker_summary_stats.csv", row.names = FALSE)

# create a horizontal boxplot of immune marker abundance
boxplot(
  markers,
  horizontal = TRUE,
  las = 1,
  col = "blue",
  border = "black",
  main = "Mean abundance of immune markers",
  xlab = "Abundance (Proportion/Frequency)",
  ylab = "Immune Markers",
  outline = TRUE
)

# add gridlines to make the marker values easier to compare
grid(
  nx = NULL,
  ny = NA,
  col = "lightgrey",
  lty = "dashed"
)

# add the mean for each marker as a red diamond
points(
  marker_means,
  seq_along(marker_means),
  pch = 18,
  col = "red",
  cex = 1.4
)

