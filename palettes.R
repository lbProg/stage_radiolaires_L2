library(RColorBrewer)
library(scales)

order_palette <- c(colorRampPalette(brewer.pal(9, "Blues"))(10),
                   "red", 
                   colorRampPalette(brewer.pal(3, "Oranges"))(4),
                   "black")

names(order_palette) <- c("Acantharea_1",
                          "Acantharea_2",
                          "Acantharea_3",
                          "Acantharea_4",
                          "Acantharea_A",
                          "Acantharea_B",
                          "Acantharea_C",
                          "Acantharea_D",
                          "Acantharea_E",
                          "Acantharea_F",
                          "Nassellaria",
                          "RAD-A_X",
                          "RAD-B_X",
                          "RAD-C_X",
                          "Radiolaria_XX",
                          "Spumellaria")
show_col(order_palette)


months_palette <- rainbow(10)
names(months_palette) <- c("1", "2", "3", "4", "5", "6", "8", "10", "11", "12")
show_col(months_palette)
