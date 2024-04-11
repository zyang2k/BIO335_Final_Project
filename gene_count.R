# Lab 9 generating the gene counts matrices for our individuals

library(ggpubr)

# Make the dataframe with your sample names and the calculated stats
## MAKE SURE YOU REPLACE THESE BEFORE RUNNING ##

stats<-as.data.frame(cbind(rbind("C1","C2","C3","HS1","HS2","HS3"),
                           rbind("Control","Control","Control","Heat_Shock","Heat_Shock","Heat_Shock"),
                           rbind(77.04,74.70,76.63,75.81,77.01,74.78),
                           rbind(73.86,79.52,88.62,92.94,71.79,83.39)))

# Add in the column names
colnames(stats)<-c("Sample_ID","Treatment","Trim_Perc","Align_Perc")

# Convert the numbers you entered to numerical values
stats$Trim_Perc<-as.numeric(stats$Trim_Perc)
stats$Align_Perc<-as.numeric(stats$Align_Perc)

# Plot the trimming percentages and view
plot1<-ggboxplot(stats, x = "Treatment", y = "Trim_Perc", add = "jitter",
                 fill = "Treatment", legend = "none",
                 ylab = "Reads retained after trimming (%)")
plot1

# Plot the alignment percentages and view
plot2<-ggboxplot(stats, x = "Treatment", y = "Align_Perc", add = "jitter",
                 fill = "Treatment", legend = "none",
                 ylab = "Reads aligned to reference (%)")
plot2

# Combine the plots and save the output (CHANGE THE OUTPUT DIRECTORY)
plot3<-ggarrange(plot1, plot2, ncol = 1, nrow = 2)
plot3
# ggsave("trim_align_stats.png", plot3,
#        height = 20, width = 10, units = "cm") # try changind the height and width to suit


#####


library(ggpubr)

library(pheatmap)
# Assuming you've set the correct path to your files
flist <- list.files(pattern = "_counts.txt$")

elist <- list()
glist <- list()

for(file in 1:length(flist)){
  full_path <- paste0("~/Desktop/BIO335/Final_Project/", flist[file])
  df <- read.table(file = full_path, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("Gene", "count")
  
  # Extract sample ID by removing "_counts.txt" from the filename
  sample_id <- gsub("_counts.txt", "", flist[file])
  df$Sample_ID <- sample_id
  
  # Add the modified dataframe to the list
  elist[[sample_id]] <- df
  
  # Extract the top gene by count and display its name and count
  df_sorted <- df[order(-df$count), ]
  second_highest_gene <- df_sorted[2, ] # Get the second row
  cat("highest gene for", sample_id, "is", second_highest_gene$Gene, "with count", second_highest_gene$count, "\n")
  # Prepare a simplified version of the data frame for gene comparison
  df2 <- df[, 1:2] # Assuming the first two columns are Gene and count
  colnames(df2) <- c("Gene", sample_id)
  glist[[sample_id]] <- df2
}


all_files<-do.call(rbind,elist)
all_files2<-do.call(cbind, glist)

all_files2<-all_files2[c(1:27091),c(1,2,4,6,8,10,12)]

biplot(princomp(all_files2[2:7]))


write.csv(all_files2, "lab9_all_counts.csv", row.names=TRUE)
