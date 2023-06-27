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


# Importation des jeux de donnees ----------

TAX <- read.csv("Data/TAX_pr2_Radiolaria_GoA.csv")
META <- read.csv("Data/META_Radiolaria_GoA.csv")
OTU <- read.csv("Data/OTU_Radiolaria_GoA.csv")

source("palettes.R")

months_equi <- data.frame("month_id" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                          "month_str" = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
periods <- data.frame("month_id" = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                      "period" = c("M_1", "M_1", "B", "B", "S", "S", "S", "S", "S", "M_2", "M_2", "M_2"))
META_srf <- META[!grepl("dcm", fixed = TRUE, META$sample_id),]


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
abund$date <- META_srf$date[match(abund$sample, META_srf$sample_id)]
abund$date <- factor(abund$date, levels = unique(abund$date))

abund$period <- GetPeriod(abund$sample)

abund$year <- paste("20", substr(abund$sample, 3, 4), sep = "")

# Stacked barplot

ggplot(data = abund, aes(x = sample, fill = order, y = value)) +
  geom_bar(stat = "identity") +
  facet_nested(cols = vars(year, period), scales = "free") +
  #guides(x = ggh4x::guide_axis_nested(delim = "&"), vjust = 1) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 10, angle = 80, hjust = 1)) +
  xlab("sample") +
  ylab("abundance") +
  ggtitle("Relative abundance of Radiolarian orders for each sampling date") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = order_palette)

# Divesite de Shannon et equitabilite de Pielou

div <- cbind("order" = TAX$Order[which(OTU$X == TAX$X)], OTU[-1])
div <- select(div, -contains("dcm"))

div <- div %>%
  group_by(order) %>%
  summarize_all(sum)

div <- as.data.frame(div)

rownames(div) <- div$order
div <- div[-1]

div <- data.frame(t(div))

div <- cbind("Shannon" = diversity(div, MARGIN = 1, index = "shannon"), div)
div <- cbind("Pielou" = div$Shannon / log(rowSums(div != 0)), div)

div <- cbind("date" = META[!grepl("dcm", fixed = TRUE, META$sample_id),]$date, div)
div <- cbind("sample" = rownames(div), div)

div$period <- GetPeriod(div$sample)
div$year <- paste("20", substr(div$sample, 3, 4), sep = "")

div$date <- factor(div$date, levels = unique(div$date)) # Sinon ordre alphabetique nuul


ggplot(data = div, aes(x = date)) +
  geom_line(aes(y = Shannon, group = 1), colour = "red") +
  geom_line(aes(y = 2*Pielou, group = 1), colour = "green") +
  facet_nested(cols = vars(year, period), scales = "free") +
  scale_y_continuous(name = "Shannon diversity", 
                     sec.axis = sec_axis(trans = ~./2, name = "Pielou evenness", breaks = seq(0, 1, 0.25))) +
  theme(axis.text.x = element_text(size = 10, angle = 80, hjust = 1),
        axis.title.y = element_text(color = "red"),
        axis.title.y.right = element_text(color = "green"))

# Courbes d'abondance relative

ab_plots <- lapply(unique(abund$order), function(i) {
  ggplot(subset(abund, order == i), aes(x = sample, y = value)) +
    geom_col() +
    facet_nested(cols = vars(year, period), scales = "free") +
    guides(fill = guide_legend(ncol = 2)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1,'lines'),
          panel.grid.major.x = element_blank(),
          strip.text.x = element_text(size = 7)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggtitle(i)
})

grid.arrange(grobs = ab_plots, ncol = 4)


# Courbes de parametres physicochimiques ----------

metadata = drop_na(META)
metadata <- subset(metadata, sample == "a")
metadata$year <- paste("20", metadata$year, sep = "")
metadata <- cbind(period = GetPeriod(metadata$sample_id), metadata)

var_plots <- lapply(unique(colnames(metadata))[-c(1:8)], function(i) {
  ggplot(data = metadata, aes(x = date, y = .data[[i]], group = depth, color = depth)) +
    geom_line() +
    geom_point() +
    facet_nested(cols = vars(year, period), scales = "free") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1,'lines'),
          panel.grid.major.x = element_blank(),
          strip.text.x = element_text(size = 7)) +
    scale_color_manual(values = c("#27a123", "black")) +
    ggtitle(i)
})

suppressMessages(grid.arrange(grobs = var_plots, ncol = 4)) # bof bof

# ACP ----------

PC_data <- drop_na(META)
PC_data <- cbind(period = GetPeriod(PC_data$sample_id), PC_data)

PCA_for_period = function(p) {
  data = subset(PC_data, period %in% p)
  # Enlever colonnes avec une variance nulle
  pca <- prcomp(data[9:23][, which(apply(data[9:23], 2, var) != 0)], scale = TRUE)
  
  a <- autoplot(pca, data = data, size = 3, loadings = TRUE, loadings.label = TRUE, loadings.label.colour = "black")
  a + ggtitle(paste(p, collapse = " - "))
}

PCA_for_period(c("B_2021", "B_2022")) # periode


# Correlogrammes

cor_table <- cbind(div[5:20], META_srf[8:22])
cor_table <- drop_na(cor_table)

cor_mat <- data.frame(matrix(NA, nrow = 16, ncol = 15))
rownames(cor_mat) <- colnames(cor_table[1:16])
colnames(cor_mat) <- colnames(cor_table[17:31])


for (var in colnames(cor_mat)) {
  for (order in rownames(cor_mat)) {
    cor_mat[order, var] <- round(cor(cor_table[order], cor_table[var]), 2)
  }
}

ggcorrplot(cor_mat) # En fonction des periodes
ggcorrplot(cor(cor_table[1:16]), hc.order = TRUE)
ggcorrplot(cor(cor_table[17:31]), hc.order = TRUE)



autocorrplots <- lapply(colnames(cor_table[1:16]), function(i) {
  a <- acf(cor_table[[i]])
  b <- with(a, data.frame(lag, acf))
  
  ggplot(data = b, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    ggtitle(i)
})

grid.arrange(grobs = autocorrplots, ncol = 4)
