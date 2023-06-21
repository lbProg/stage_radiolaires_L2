# Importation des librairies ----------

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggfortify)


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

abund <- cbind("order" = TAX$Order[which(OTU_norm$X == TAX$X)], OTU_norm[-1])
abund <- gather(abund, "sample", "value", -1)

abund <- abund %>%
  group_by(sample, order) %>%
  summarize(across(value, sum))

abund$month = substr(abund$sample, 5, 6)

ggplot(data = abund, aes(x = sample, fill = order, y = value)) +
  facet_wrap(~month, scales="free") +
  geom_bar(stat = "identity") +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(legend.text = element_text(size = 10)) +
  guides(fill = guide_legend(ncol = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(guide = guide_axis(angle = 80)) +
  scale_fill_manual(values = order_palette)


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

PCA_for_month(unique(META$month)) # Numero de mois ou plusieurs mois avec c()
