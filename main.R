# Importation des librairies ----------

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggfortify)
library(gridExtra)


# Importation des jeux de donnees ----------

TAX <- read.csv("Data/TAX_pr2_Radiolaria_GoA.csv")
META <- read.csv("Data/META_Radiolaria_GoA.csv")
OTU <- read.csv("Data/OTU_Radiolaria_GoA.csv")

source("palettes.R")


# Verification de la correspondnce des tableaux entre eux ----------

# Correspondance des OTUs
TAX_to_OTU <- data.frame(OTU[1], TAX[1])
TAX_to_OTU$match <- TAX_to_OTU[1] == TAX_to_OTU[2]

length(TAX_to_OTU$match)
which(TAX_to_OTU$match == FALSE)

# Correspondance des samples
META_to_OTU <- data.frame(META[1], colnames(OTU[-1]))
META_to_OTU$match <- META_to_OTU[1] == META_to_OTU[2]

length(META_to_OTU$match)
which(META_to_OTU$match == FALSE)

# Correspondance des replicats entre eux (META)
META_match = as.data.frame(META[META$sample == "a", ] == META[META$sample == "b", ])
which(META_match[c(-1, -6)] == FALSE) # Les 2 réplicats (a/b) sont exactement identiques partout


# Description des jeux de données ----------

NumberOf = function(class) {
  return(length(unique(TAX$Family[TAX$Class == class])))
}

summ <- data.frame("nb_OTUs" = nrow(OTU),
                  "max_abund" = max(OTU[-1]),
                  "mean_abund" = mean(as.matrix(OTU[-1])),
                  "nb_samples" = ncol(OTU[-1]),
                  "nb_days" = length(unique(substr(colnames(OTU[-1]), 3, 8))),
                  "nb_variables" = ncol(META[7:22]),
                  "nb_families" = length(unique(TAX$Family)),
                  "nb_Acantharea" = NumberOf("Acantharea"),
                  "nb_Polycystinea" = NumberOf("Polycystinea"),
                  "nb_RAD-A" = NumberOf("RAD-A"),
                  "nb_RAD-B" = NumberOf("RAD-B"),
                  "nb_RAD-C" = NumberOf("RAD-C"),
                  "nb_Radiolaria_X" = NumberOf("Radiolaria_X"))

unique(TAX$Class)


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

abund$month <- substr(abund$sample, 5, 6)
abund$month[abund$month != 10] <- gsub("0", "", abund$month[abund$month != 10])
abund$month <- factor(abund$month, levels = c("1", "2", "3", "4", "5", "6", "8", "10", "11", "12"))


# Stacked barplot

ggplot(data = abund, aes(paste0(sample, "&", month), fill = order, y = value)) +
  geom_bar(stat = "identity") +
  guides(x = ggh4x::guide_axis_nested(delim = "&"), vjust = 1) +
  guides(fill = guide_legend(ncol = 2)) +
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 80, hjust = 1)) +
  xlab("sample") +
  ggtitle("Relative abundance of Radiolarian orders for each sampling date") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = order_palette)

# Courbes d'abondance relative

ab_plots <- lapply(unique(abund$order), function(i) {
  ggplot(subset(abund, order == i), aes(x = sample, y = value, fill = month)) +
    geom_col() +
    guides(fill = guide_legend(ncol = 2)) +
    theme(axis.text.x = element_blank()) +
    ylab("relative abundance") +
    ggtitle(i) +
    scale_fill_manual(values = months_palette)
})

grid.arrange(grobs = ab_plots, ncol=4)

PlotAbundance = function(order) {
  ggplot(data = abund[abund$order == order, ], aes(x = sample, y = value, fill = month)) +
    geom_col() +
    theme(axis.text.x = element_text(size = 8, angle = 80, hjust = 1)) +
    ylab("relative abundance") +
    scale_fill_manual(values = months_palette)
}

PlotAbundance("Acantharea_B")


# Courbes de parametres physicochimiques ----------

var_plots <- lapply(unique(colnames(META))[-c(1:7)], function(i) {
  ggplot(data = META, aes(x = sample_id, y = .data[[i]], group = 1)) +
    geom_line(aes(color = factor(month))) +
    geom_point(aes(color = factor(month))) +
    theme(axis.text.x = element_blank()) +
    ylab("value") +
    labs(color = "month") +
    ggtitle(i) +
    scale_color_manual(values = months_palette)
})

grid.arrange(grobs = var_plots, ncol=4)


# ACP ----------

PCA_for_month = function(month) {
  PC_data <- drop_na(META)
  PC_data <- PC_data[PC_data$month %in% month, ]
  
  # Enlever colonnes avec une variance nulle
  pca <- prcomp(PC_data[8:22][, which(apply(PC_data[8:22], 2, var) != 0)], scale = TRUE)
  PC_data$month = as.factor(PC_data$month)
  
  a <- autoplot(pca, data = PC_data, colour = "month", size = 3, loadings = TRUE, loadings.label = TRUE, loadings.label.colour = "black")
  a + scale_colour_manual(values = months_palette)
}

PCA_for_month(c(5, 6, 8, 10)) # Numero de mois ou plusieurs mois avec c()
