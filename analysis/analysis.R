# Libraries
library(data.table)
library(magrittr)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scales)
library(tools)

# Settings
setwd("/Users/edsouza")
setDTthreads(0)

# Read and process data
dt <- fread("possible_variants_all_chrs.txt")
dt[is.na(methyl_level), methyl_level := 0]
dt[, frequency_status := ifelse(
  observed_in_gnomad == TRUE,
  ifelse(info_ac == 1, 'Singleton', 'AC > 1'),
  'Not observed in gnomAD'
)]

# Create a summary version of the dataset
summary_dt <- dt[,
                 .N,
                 by = .(observed_in_gnomad,
                        methyl_level,
                        consequence)]
summary_dt[, consequence := toTitleCase(gsub("PRIME", "'", gsub("_", " ", consequence)))]
setorder(summary_dt,
         observed_in_gnomad,
         consequence,
         methyl_level)


# Plotting
ggplot(summary_dt,
       aes(x = methyl_level,
           y = N,
           fill = observed_in_gnomad)) +
  geom_col() +
  facet_wrap(~ consequence, scales = 'free_y') +
  theme_minimal() +
  labs(x = "Methylation Level",
       y = "Count of CpG Variants",
       fill = "Variant Observed") +  scale_fill_brewer(palette = "Set1")  +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA)
    
  )

ggsave(
  'plot_consequences_observed_gnomad.png',
  p,
  width = 10,
  height = 10,
  dpi = 300
)

meth_summary_dt <- dt[,
                      .N,
                      by = .(frequency_status,
                             methyl_level,
                             consequence)]

# Calculate the sum of N for each combination of methyl_level and consequence
total_N <-
  meth_summary_dt[, .(total_N = sum(N)), by = .(methyl_level)]

# Join the sum back to the main data.table and calculate the proportion
meth_summary_dt[total_N, proportion := N / total_N, on = .(methyl_level)]

# Create the frequency status as an ordered factor
meth_summary_dt$frequency_status <-
  factor(
    meth_summary_dt$frequency_status,
    levels = c('Not observed in gnomAD', 'Singleton', 'AC > 1')
  )

# Plot
ggplot(meth_summary_dt,
       aes(x = methyl_level, y = proportion, fill = frequency_status)) +
  geom_col() + theme_minimal() +
  labs(x = "Methylation Level",
       y = "Fraction of variant observed",
       fill = "Frequency status") +  scale_fill_brewer(palette = "Paired") +
  ggtitle("Frequency of C->T variants in gnomAD by methylation level") +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA)
  )




# Calculate proportions of observed and not observed in gnomAD


prop_dt <- summary_dt[methyl_level > 4
                      , .(num_observed = sum(.SD[observed_in_gnomad == T]$N),
                          num_not_observed = sum(.SD[observed_in_gnomad == F]$N)),
                      by = .(consequence)]

prop_dt[, prop := num_observed / (num_observed + num_not_observed)]
prop_dt[, lower_ci :=
          prop - 1.96 * sqrt(prop * (1 - prop) / (num_observed + num_not_observed))]
prop_dt[, upper_ci :=
          prop + 1.96 * sqrt(prop * (1 - prop) / (num_observed + num_not_observed))]
prop_dt[, total_num := num_observed + num_not_observed]
setorder(prop_dt, prop)
prop_dt$consequence <- factor(prop_dt$consequence,
                              levels = unique(prop_dt$consequence))

p <- ggplot(prop_dt, aes(x = consequence, y = prop)) +
  geom_col(
    aes(fill = consequence),
    width = 0.6,
    color = "grey30",
    size = 0.2
  ) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.2,
    color = "grey30",
    size = 0.5
  ) +
  geom_text(
    aes(
      y = 0,
      label = paste(comma(num_observed), " (", comma(total_num), ")", sep = "")
    ),
    hjust = "inward",
    color = "black",
    size = 3.5
  ) +
  coord_flip() +
  labs(title = "Fraction of methylated CpG sites with an observed variant in gnomAD",
       x = NULL,
       y = NULL) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_fill_brewer(palette = "Set3")

print(p)

# Save the plot as an image
ggsave(
  'pretty_proportions_plot.png',
  p,
  width = 10,
  height = 10,
  dpi = 150
)




prop_ac_dt <- meth_summary_dt[methyl_level > 4
                              , .(
                                num_observed_singleton = sum(.SD[frequency_status == 'Singleton']$N),
                                num_observed_ac_geq_1 = sum(.SD[frequency_status == 'AC > 1']$N),
                                num_not_observed = sum(.SD[frequency_status == 'Not observed in gnomAD']$N)
                              ),
                              by = .(consequence)]

prop_ac_dt[, total_num := (num_observed_singleton + num_not_observed + num_observed_ac_geq_1)]
prop_ac_dt[,
           `:=` (
             prop_singleton = num_observed_singleton / total_num,
             prop_ac_geq_1 = num_observed_ac_geq_1 / total_num,
             prop_unobserved = num_not_observed / total_num
           )]

# Melt and clean up data.table
prop_ac_dt_long <- melt(
  prop_ac_dt[, .(consequence,
                 prop_singleton,
                 prop_ac_geq_1,
                 prop_unobserved)],
  id.vars = 'consequence',
  measure.vars = c("prop_singleton",
                   "prop_ac_geq_1", "prop_unobserved")
)

names(prop_ac_dt_long) <-
  c("consequence", 'frequency_level', 'proportion')
prop_ac_dt_long[, consequence := toTitleCase(gsub("PRIME", "'", gsub("_", " ", consequence)))]
prop_ac_dt_long$frequency_status = factor(
  prop_ac_dt_long$frequency_level,
  levels = c('prop_unobserved', 'prop_ac_geq_1', 'prop_singleton'),
  labels = c('Unobserved in gnomAD', 'AC > 1', 'Singleton'),
  ordered = T
)

# Create an order for the plot
prop_ac_dt_long_ord <-
  prop_ac_dt_long[frequency_level == 'prop_unobserved', .(consequence, ord = frank(proportion))]
setkey(prop_ac_dt_long_ord, consequence)
setkey(prop_ac_dt_long, consequence)
prop_ac_dt_long <- prop_ac_dt_long_ord[prop_ac_dt_long]

# Data visualization
ggplot(prop_ac_dt_long,
       aes(
         x = reorder(consequence, -ord),
         y = proportion,
         fill = frequency_status
       )) +
  
  # Bar chart
  geom_col(position = "dodge") +
  
  # Text labels
  geom_text(aes(
    y = proportion,
    label = scales::percent(proportion, accuracy = 0.01)
  ),
  position = position_dodge(width = 0.9)) +
  
  # Convert to horizontal bar chart
  coord_flip() +
  
  # Labels and titles
  labs(title = "Fraction of methylated CpG sites with an observed variant in gnomAD",
       x = NULL,
       y = 'Fraction of variants',
       fill = 'Frequency Status') +
  
  # Theming and color palette
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  scale_fill_brewer(palette = "Set1")


#########



# Sample a 1000 positions for each group and write to file
# Separate them by methylation level
sampled_dt <- dt[,
                 .SD[sample(.N, 3000, replace = T)], by = .(consequence, observed_in_gnomad)]

# Plot GERP, PhyloP, and PhastCons scores between observed and not observed
# in gnomAD for each consequence

# Reshape the data for plotting
melted_data <- melt(
  sampled_dt[methyl_level > 4],
  id.vars = c("consequence", "observed_in_gnomad"),
  measure.vars = c("gerp_rs", "mam_phylop", "cons_score"),
  variable.name = "Score_Type",
  value.name = "Score"
)
# Faceted plot
# Faceted plot
melted_data[, consequence := toTitleCase(gsub("PRIME", "'", gsub("_", " ", consequence)))]

# Faceted plot with enhancements
ggplot(melted_data,
       aes(x = consequence, y = Score, fill = observed_in_gnomad)) +
  geom_boxplot() +  # this hides the outliers; remove if you want to show them
  labs(
    title = "Distribution of Evolutionary Scores by Consequence (Methylation Level > 4)",
    subtitle = "Separated by GERP++, PhyloP and
    CADD Scores",
    y = "Score Value",
    x = "Consequence",
    fill = "Observed in gnomAD"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5)
  ) +
  facet_wrap(~ Score_Type, scales = "free_y", ncol = 1)
# Save the plot as an image
ggsave(
  'conservation_plot.png',
  cons_plot,
  width = 10,
  height = 10,
  dpi = 200
)


# Compute
effect_sizes <-
  melted_data[, .(cohens_d = (mean(Score[observed_in_gnomad == TRUE], na.rm = TRUE) -
                                mean(Score[observed_in_gnomad == FALSE], na.rm = TRUE)) /
                    sqrt(((.N - 1) * sd(Score[observed_in_gnomad == TRUE], na.rm = TRUE) ^
                            2 +
                            (.N - 1) * sd(Score[observed_in_gnomad == FALSE], na.rm = TRUE) ^
                            2
                    ) /
                      (2 * (.N - 1)))),
              by = consequence]
