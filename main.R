# A faire
# - correllogrammes par especes et par grands groupes (Nassellaires Spumellaires Acanthaires Gpes env)
# - echelle de temps lineaire sur plots (breaks la ou echantillons pour montrer abondance 0)
# - montrer dcm sur plots d'abondance individuelle
# - faire 1 - Pielou ?
# - evolution d'un indice de similarite entre les 2 replicats d'abondance totale
# - regrouper parametres sur correlogrammes
# - heatmap par periode
# - NMDS


# Importation des librairies ----------

library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyft)
library(ggfortify)
library(gridExtra)
library(vegan)
library(ggcorrplot)
library(tseries)
library(ggh4x)
library(ggpubr)
library(pheatmap)


# Fonctions utiles par la suite ----------

NumberOf = function(class) {
  # Returns the number of different families present in a Radiolarian class
  return(length(unique(TAX$Family[TAX$Class == class])))
}

GetPeriod = function(samples) {
  # Generates a period column assigning a period to each sample of the dataframe
  m <- as.numeric(substr(samples, 5, 6))
  p <- paste(periods$period[match(m, periods$month_id)], "_20", substr(samples, 3, 4), sep = "")
  
  p <- factor(p, levels = c("M_2_2020", "M_1_2021", "B_2021", "S_2021", "M_2_2021", "M_1_2022", "B_2022", "S_2022"))
  
  return(p)
}

GetMonth = function(samples) {
  # Generates a month id column
  m <- substr(samples, 5, 6)
  m[m != "10"] <- gsub("0", "", m[m != "10"])
  
  m <- factor(m, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
  
  return(m)
}

# Importation des jeux de donnees ----------

TAX <- read.csv("Data/TAX_pr2_Radiolaria_GoA.csv")
META <- read.csv("Data/META_Radiolaria_GoA.csv")
OTU <- read.csv("Data/OTU_Radiolaria_GoA.csv")

TAX <- cbind(OTU_name = paste("OTU_", sprintf("%03d", seq(1, 531, 1)), sep = ""), TAX)
META <- subset(META, select = -c(phaeo, Si_TIN, Si_P, TIN_P, air_temp, solar_radiation))

source("palettes.R")

months_equi <- data.frame("month" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                          "month_str" = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
periods <- data.frame("month_id" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                      "period" = c("M_1", "M_1", "B", "B", "S", "S", "S", "S", "S", "M_2", "M_2", "M_2"))

period_labs <- c("M", "M", "B", "S", "M", "M", "B", "S")
names(period_labs) <- c("M_2_2020", "M_1_2021", "B_2021", "S_2021",
                        "M_2_2021", "M_1_2022", "B_2022", "S_2022")


OTU_to_group <- data.frame(OTU = TAX$OTU_name, Class = TAX$Class, Order = TAX$Order, group = "")
OTU_to_group$group <- OTU_to_group$Class
OTU_to_group[OTU_to_group$Class == "Polycystinea", ]$group <- OTU_to_group[OTU_to_group$Class == "Polycystinea", ]$Order


# Description des jeux de données ----------

summ <- data.frame("nb_OTUs" = nrow(OTU),
                  "max_abund" = max(OTU[-1]),
                  "mean_abund" = mean(as.matrix(OTU[-1])),
                  "nb_samples" = ncol(OTU[-1]),
                  "nb_days" = length(unique(substr(colnames(OTU[-1]), 3, 8))),
                  "nb_variables" = ncol(META[8:16]),
                  "nb_families" = length(unique(TAX$Family)),
                  "nb_Acantharea" = NumberOf("Acantharea"),
                  "nb_Polycystinea" = NumberOf("Polycystinea"),
                  "nb_RAD-A" = NumberOf("RAD-A"),
                  "nb_RAD-B" = NumberOf("RAD-B"),
                  "nb_RAD-C" = NumberOf("RAD-C"),
                  "nb_Radiolaria_X" = NumberOf("Radiolaria_X"))


# Normalisation du tableau OTU (abondance) ----------

OTU_norm <- mapply('/', OTU[-1], colSums(OTU[, -1]))
OTU_norm <- cbind(OTU[1], OTU_norm)


# Tableau d'abondance ----------

abund <- cbind("class" = TAX$Class[which(OTU_norm$X == TAX$X)],
               "order" = TAX$Order[which(OTU_norm$X == TAX$X)],
               OTU_norm[-1])

abund <- gather(abund, "sample", "value", -c(1, 2))

abund <- abund %>%
  group_by(sample, class, order) %>%
  summarize(across(value, sum))

abund$depth <- "surface"
abund$depth[substring(abund$sample, 10, 12) == "dcm"] <- "dcm"

abund$month <- GetMonth(abund$sample)
abund$period <- GetPeriod(abund$sample)

abund$year <- paste("20", substr(abund$sample, 3, 4), sep = "")

abund <- abund %>%
  group_by(order, month, year, period, depth) %>%
  summarize(across(value, mean))

for (i in unique(abund$order)) {
  abund <- rbind(abund, data.frame(order = i, month = "7", year = "2021", period = "S_2021", value = 0, depth = "all"))
  abund <- rbind(abund, data.frame(order = i, month = "9", year = "2021", period = "S_2021", value = 0, depth = "all"))
}

abund$period <- factor(abund$period, levels = names(period_labs))

abund$month_str <- months_equi$month_str[match(abund$month, months_equi$month)]
abund$month_str <- factor(abund$month_str, levels = months_equi$month_str)

# Stacked barplot

ggplot(data = subset(abund, depth %in% c("surface")), aes(x = month_str, fill = order, y = value)) +
  geom_bar(stat = "identity") +
  facet_nested(cols = vars(year, period), space = "free", scales = "free", labeller = labeller(period = period_labs)) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 80, hjust = 0, vjust = 0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.1, "lines"),) +
  xlab("month") +
  ylab("abundance") +
  ggtitle("Relative abundance of Radiolarian orders for each period") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = order_palette)


# Divesite de Shannon et equitabilite de Pielou

div <- data.frame("Shannon" = diversity(t(OTU[-1]), MARGIN = 1, index = "shannon"))
div$Pielou <- div$Shannon / log(rowSums(t(OTU[-1]) != 0))

div$month <- GetMonth(rownames(div))
div$period <- GetPeriod(rownames(div))
div$year <- paste("20", substr(rownames(div), 3, 4), sep = "")

div <- drop_na(div)

div <- div %>%
  group_by(month, period, year) %>%
  summarize(across(c(Shannon, Pielou), mean))

div <- rbind(div, data.frame(month = "7", period = "S_2021", year = "2021", Shannon = 0, Pielou = 0))
div <- rbind(div, data.frame(month = "9", period = "S_2021", year = "2021", Shannon = 0, Pielou = 0))

div$month_str <- months_equi$month_str[match(div$month, months_equi$month)]
div$month_str <- factor(div$month_str, levels = months_equi$month_str)

div$period <- factor(div$period, levels = names(period_labs))

div <- gather(div, "index", "value", -c(1, 2, 3, 6))

ggplot(data = div, aes(x = month_str, y = value, fill = index)) +
  geom_col(position = position_dodge()) +
  facet_nested(cols = vars(year, period), scales = "free", space = "free", labeller = labeller(period = period_labs)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("index value") +
  theme(axis.text.x = element_text(size = 10, angle = 80, hjust = 0, vjust = 0),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())


# Courbes d'abondance relative

ab_plots <- lapply(unique(abund$order), function(i) {
  ggplot(subset(abund, order == i), aes(x = month_str, y = value, fill = depth)) +
    geom_bar(stat = "identity") +
    facet_nested(cols = vars(year, period), scales = "free", space = "free", labeller = labeller(period = period_labs)) +
    theme(axis.text.x = element_text(size = 7, angle = 80, hjust = 0, vjust = 0),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1,'lines'),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.text.x = element_text(size = 7)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0, 0)) +
    scale_fill_manual(values = depth_palette) +
    ggtitle(i)
})

grid.arrange(grobs = ab_plots, ncol = 4, nrow = 4)

# Courbes de parametres physicochimiques ----------

metadata = drop_na(META)
metadata <- subset(metadata, sample == "a")
metadata$year <- paste("20", metadata$year, sep = "")
metadata <- cbind(period = GetPeriod(metadata$sample_id), metadata)

var_plots <- lapply(unique(colnames(metadata))[-c(1:8)], function(i) {
  ggplot(data = metadata, aes(x = date, y = .data[[i]], group = depth, color = depth)) +
    geom_line() +
    geom_point() +
    facet_nested(cols = vars(year, period), scales = "free", labeller = labeller(period = period_labs)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1,'lines'),
          panel.grid.major.x = element_blank(),
          strip.text.x = element_text(size = 7)) +
    scale_color_manual(values = depth_palette) +
    ggtitle(i)
})

suppressMessages(grid.arrange(grobs = var_plots, ncol = 3)) # bof bof

# ACP ----------

PC_data <- drop_na(META)
PC_data <- cbind(period = GetPeriod(PC_data$sample_id), PC_data)

PCAForPeriod = function(p, title = paste(p, collapse = " - ")) {
  data = subset(PC_data, period %in% p)
  # Enlever colonnes avec une variance nulle
  pca <- prcomp(data[9:17][, which(apply(data[9:17], 2, var) != 0)], scale = TRUE)
  
  a <- autoplot(pca, data = data, size = 3, colour = "depth", loadings = TRUE, loadings.label = TRUE, loadings.label.colour = "black")
  return(a + ggtitle(title) + scale_color_manual(values = depth_palette))
}

mixing_PCA <- PCAForPeriod(c("M_2_2020", "M_1_2021", "M_2_2021", "M_1_2022"), "mixing")
bloom_PCA <- PCAForPeriod(c("B_2021", "B_2022"), "bloom")
strat_PCA <- PCAForPeriod(c("S_2021", "S_2022"), "stratification")
all_PCA <- PCAForPeriod(unique(PC_data$period), "all periods")


ggarrange(mixing_PCA, bloom_PCA, strat_PCA, all_PCA, ncol = 2, nrow = 2)

# Correlogrammes

correlation_table <- gather(OTU, "sample", "value", -1)
correlation_table$X <- paste("OTU_", sprintf("%03d", seq(1, 531, 1)), sep = "")
correlation_table$order <- OTU_to_group$Order[match(correlation_table$X, OTU_to_group$OTU)]
correlation_table$group <- OTU_to_group$group[match(correlation_table$X, OTU_to_group$OTU)]

correlation_table$period <- GetPeriod(correlation_table$sample)
correlation_table <- cbind(correlation_table, META[match(correlation_table$sample, META$sample_id), ][7:16])

correlation_table <- drop_na(correlation_table)


ComputeCorrelation = function(ind_table, var_table) {
  df <- data.frame(matrix(NA, nrow = ncol(ind_table), ncol = ncol(var_table)))
  colnames(df) <- colnames(var_table)
  rownames(df) <- colnames(ind_table)
  
  for (var in colnames(df)) {
    for (ind in rownames(df)) {
      df[ind, var] <- round(cor(ind_table[ind], var_table[var]), 2)
    }
  }
  
  return(df)
}


# Par OTU
PlotOTUCorr = function(table, periods, depths) {
  table <- subset(table, period %in% periods & depth %in% depths)
  corr <- ComputeCorrelation(as.data.frame(df_mat(table, sample, X, value)),
                             table[match(unique(table$sample), table$sample), ][8:16])
  corr <- cbind(group = table$group[match(rownames(corr), table$X)], corr)
  corr <- corr[order(corr$group), ]
  
  pheatmap(t(corr[-1]),
           annotation_col = corr[1],
           annotation_colors = list(group = group_palette),
           main = paste(paste(periods, collapse = " - "), " / ", paste(depths, collapse = " - ")),
           show_colnames = FALSE,
           cluster_cols = FALSE,
           cluster_rows = FALSE)
}

PlotOTUCorr(correlation_table, c("M_2_2020", "M_1_2021", "M_2_2021", "M_1_2022"), c("surface"))
PlotOTUCorr(correlation_table, c("B_2021", "B_2022"), c("surface"))
PlotOTUCorr(correlation_table, c("S_2021", "S_2022"), c("dcm"))


# Par ordre
PlotOrderCorr = function(table, periods, depths) {
  table <- subset(table, period %in% periods & depth %in% depths) %>%
    group_by(sample, order, period, depth, temp, salinity, chl, PO4, Silica, TIN, wind_speed, PAR, UV) %>%
    summarize(across(value, sum))
  corr <- ComputeCorrelation(as.data.frame(df_mat(table, sample, order, value)),
                             as.data.frame(df_mat(table, sample, order, value)))
  pheatmap(corr,
           main = paste(paste(periods, collapse = " - "), " / ", paste(depths, collapse = " - ")),
           cluster_cols = FALSE,
           cluster_rows = FALSE)
}


PlotOrderCorr(correlation_table, c("M_2_2020", "M_1_2021", "M_2_2021", "M_1_2022"), c("surface"))
PlotOrderCorr(correlation_table, c("B_2021", "B_2022"), c("surface"))
PlotOrderCorr(correlation_table, c("S_2021", "S_2022"), c("dcm"))


# Par variable
PlotVariableCorr = function(table, periods, depths) {
  table <- subset(table, period %in% periods & depth %in% depths)
  corr <- ComputeCorrelation(table[match(unique(table$sample), table$sample), ][8:16],
                             table[match(unique(table$sample), table$sample), ][8:16])
  pheatmap(corr,
           main = paste(paste(periods, collapse = " - "), " / ", paste(depths, collapse = " - ")),
           cluster_cols = FALSE,
           cluster_rows = FALSE)
}

PlotVariableCorr(correlation_table, c("M_2_2020", "M_1_2021", "M_2_2021", "M_1_2022"), c("surface"))
PlotVariableCorr(correlation_table, c("B_2021", "B_2022"), c("surface"))
PlotVariableCorr(correlation_table, c("S_2021", "S_2022"), c("dcm"))
