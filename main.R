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

TAX <- cbind(OTU_name = paste("OTU_", seq(1, 531, 1), sep = ""), TAX)
META <- select(META, -c("phaeo", "Si_TIN", "Si_P", "TIN_P", "air_temp", "solar_radiation"))
META_srf <- META[!grepl("dcm", fixed = TRUE, META$sample_id),]

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


# Creation d'un tableau avec les ordres en colonnes et les echantillons en lignes

orders_table <- cbind("order" = TAX$Order[which(OTU$X == TAX$X)], OTU[-1])
orders_table <- select(orders_table, -contains("dcm"))

orders_table <- orders_table %>%
  group_by(order) %>%
  summarize_all(sum)

orders_table <- as.data.frame(orders_table)

rownames(orders_table) <- orders_table$order
orders_table <- orders_table[-1]

orders_table <- data.frame(t(orders_table))

# Normalisation du tableau OTU (abondance) ----------

OTU_norm <- mapply('/', OTU[-1], colSums(OTU[, -1]))
OTU_norm <- cbind(OTU[1], OTU_norm)


# Tableau d'abondance ----------

abund <- cbind("class" = TAX$Class[which(OTU_norm$X == TAX$X)],
               "order" = TAX$Order[which(OTU_norm$X == TAX$X)],
               OTU_norm[-1])
abund <- select(abund, -contains("dcm"))

abund <- gather(abund, "sample", "value", -c(1, 2))

abund <- abund %>%
  group_by(sample, class, order) %>%
  summarize(across(value, sum))

abund$month <- GetMonth(abund$sample)
abund$period <- GetPeriod(abund$sample)

abund$year <- paste("20", substr(abund$sample, 3, 4), sep = "")

abund <- abund %>%
  group_by(order, month, year, period) %>%
  summarize(across(value, mean))

for (i in unique(abund$order)) {
  abund <- rbind(abund, data.frame(order = i, month = "7", year = "2021", period = "S_2021", value = 0))
  abund <- rbind(abund, data.frame(order = i, month = "9", year = "2021", period = "S_2021", value = 0))
}

abund$period <- factor(abund$period, levels = names(period_labs))

abund$month_str <- months_equi$month_str[match(abund$month, months_equi$month)]
abund$month_str <- factor(abund$month_str, levels = months_equi$month_str)

# Stacked barplot

ggplot(data = abund, aes(x = month_str, fill = order, y = value)) +
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

div <- data.frame("Shannon" = diversity(orders_table, MARGIN = 1, index = "shannon"))
div$Pielou <- div$Shannon / log(rowSums(orders_table != 0))

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
  ggplot(subset(abund, order == i), aes(x = month_str, y = value)) +
    geom_col() +
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
  pca <- prcomp(data[9:23][, which(apply(data[9:23], 2, var) != 0)], scale = TRUE)
  
  a <- autoplot(pca, data = data, size = 3, colour = "depth", loadings = TRUE, loadings.label = TRUE, loadings.label.colour = "black")
  return(a + ggtitle(title) + scale_color_manual(values = depth_palette))
}

mixing_PCA <- PCAForPeriod(c("M_2_2020", "M_1_2021", "M_2_2021", "M_1_2022"), "mixing")
bloom_PCA <- PCAForPeriod(c("B_2021", "B_2022"), "bloom")
strat_PCA <- PCAForPeriod(c("S_2021", "S_2022"), "stratification")
all_PCA <- PCAForPeriod(unique(PC_data$period), "all periods")


ggarrange(mixing_PCA, bloom_PCA, strat_PCA, all_PCA, ncol = 2, nrow = 2)

# Correlogrammes

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

PlotCorrelogram = function(matrix, title = "") {
  p <- ggcorrplot(matrix, title = title) +
    theme(axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 6),
          legend.position = "none")
  
  return(p)
}

variables_corr_table <- META_srf[8:16]
rownames(variables_corr_table) <- META_srf$sample_id
variables_corr_table <- drop_na(variables_corr_table)
variables_corr_table <- cbind("period" = GetPeriod(rownames(variables_corr_table)), variables_corr_table)

orders_corr_table <- orders_table
orders_corr_table$NA_DROP <- META_srf$TIN
orders_corr_table <- drop_na(orders_corr_table)
orders_corr_table <- subset(orders_corr_table, select = -NA_DROP)
orders_corr_table <- cbind("period" = GetPeriod(rownames(orders_corr_table)), orders_corr_table)

orders_corr_plot <- ggcorrplot(cor(orders_corr_table[-1]), hc.order = TRUE, title = "orders")
variables_corr_plot <- ggcorrplot(cor(variables_corr_table[-1]), hc.order = TRUE, title = "variables")

ggarrange(orders_corr_plot, variables_corr_plot, ncol = 2, nrow = 1)


# Pour chaque periode
corr_mixing <- ComputeCorrelation(subset(orders_corr_table, period %in% c("M_2_2020", "M_1_2021", "M_2_2021", "M_1_2022"))[-1],
                                  subset(variables_corr_table, period %in% c("M_2_2020", "M_1_2021", "M_2_2021", "M_1_2022"))[-1])
corr_m_plot <- PlotCorrelogram(corr_mixing, "mixing")

corr_bloom <- ComputeCorrelation(subset(orders_corr_table, period %in% c("B_2021", "B_2022"))[-1],
                                 subset(variables_corr_table, period %in% c("B_2021", "B_2022"))[-1])
corr_b_plot <- PlotCorrelogram(corr_bloom, title = "bloom")

corr_strat <- ComputeCorrelation(subset(orders_corr_table, period %in% c("S_2021", "S_2022"))[-1],
                                 subset(variables_corr_table, period %in% c("S_2021", "S_2022"))[-1])
corr_s_plot <- PlotCorrelogram(corr_strat, title = "stratification")

corr_all <- ComputeCorrelation(orders_corr_table[-1],
                               variables_corr_table[-1])
corr_all_plot <- PlotCorrelogram(corr_all, title = "all periods")

ggarrange(corr_m_plot, corr_b_plot, corr_s_plot, corr_all_plot, ncol = 2, nrow = 2)

# Par OTU
OTU_corr_table <- data.frame(select(OTU, -contains("dcm"))[-1])
rownames(OTU_corr_table) <- paste("OTU_", seq(1, 531, 1), sep = "")
OTU_corr_table <- rbind(TIN_NA_DROP = META_srf$TIN, OTU_corr_table)
OTU_corr_table <- OTU_corr_table[, colSums(is.na(OTU_corr_table)) == 0]
OTU_corr_table <- OTU_corr_table[-1, ]

OTU_corr_table <- cbind(OTU_name = rownames(OTU_corr_table), OTU_corr_table)

OTU_corr_table <- cbind(group = OTU_to_group$group[match(OTU_to_group$OTU, rownames(OTU_corr_table))], OTU_corr_table)
OTU_corr_table$group <- as.factor(OTU_corr_table$group)


OTU_correlation <- ComputeCorrelation(as.data.frame(t(OTU_corr_table[, -c(1, 2)])), as.data.frame(drop_na(META_srf[8:16])))
OTU_correlation <- cbind(group = OTU_corr_table$group[match(OTU_corr_table$OTU_name, rownames(OTU_correlation))], OTU_correlation)

OTU_correlation <- drop_na(OTU_correlation)

OTU_correlation <- OTU_correlation[order(OTU_correlation$group), ]

pheatmap(t(OTU_correlation[-1]),
         annotation_col = OTU_correlation[1], 
         annotation_colors = list(group = group_palette), 
         show_colnames = FALSE,
         cluster_cols = FALSE)

OTU_norm[OTU_norm$X == "df70605d2e764f7e03a6814509d358d9", ]
