setwd("Z://Users//RMunn//R_analysis//Speed")

library('tidyverse')
library('data.table')
library('ggplot2')

Open_Field <- read.csv(file = "Z://Users//RMunn//R_analysis//Speed//Speed_Cells_OF.csv", header = TRUE, sep=",")
Squish <- read.csv(file = "Z://Users//RMunn//R_analysis//Speed//Speed_Cells_Squish.csv", header = TRUE, sep=",")

Score_Merge <- data.frame("Speed_Score_OF" = Open_Field$Speed_Score_2, 
                  "Speed_Score_SQ" = Squish$Speed_Score_2, 
                  "Stability_OF" = Open_Field$Speed_Stability, 
                  "Stability_SQ" = Squish$Speed_Stability,
                  "Speed_Score_OF_X" = Open_Field$Speed_Score_X,
                  "Speed_Score_SQ_X" = Squish$Speed_Score_X,
                  "Speed_Score_OF_Y" = Open_Field$Speed_Score_Y,
                  "Speed_Score_SQ_Y" = Squish$Speed_Score_Y)
Slope_Merge <- data.frame("Slope_OF" = Open_Field$Slope,
                          "Slope_SQ" = Squish$Slope,
                          "Slope_OF_X" = Open_Field$Slope_X,
                          "Slope_SQ_X" = Squish$Slope_X,
                          "Slope_OF_Y" = Open_Field$Slope_Y,
                          "Slope_SQ_Y" = Squish$Slope_Y)

dev.off()

Score = ggplot(Score_Merge, aes(Speed_Score_OF,Speed_Score_SQ, color = Stability_OF)) +
  geom_point(shape = 2) +
  geom_smooth(method = lm)+
  xlab("Speed Score in the Open Field") +
  ylab("Speed Score in Environmental Compression")

Score_X = ggplot(Score_Merge, aes(Speed_Score_OF_X,Speed_Score_SQ_X, color = Stability_OF)) +
  geom_point(shape = 2) +
  geom_smooth(method = lm) +
  xlab("Speed Score in the X axis of the Open Field") +
  ylab("Speed Score in X axis of the Environmental Compression")

Score_Y = ggplot(Score_Merge, aes(Speed_Score_OF_Y,Speed_Score_SQ_Y, color = Stability_OF)) +
  geom_point(shape = 2) +
  geom_smooth(method = lm) +
  xlab("Speed Score in the Y axis of the Open Field") +
  ylab("Speed Score in Y axis of the Environmental Compression")


   