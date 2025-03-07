#ADDING THE METADATA
#loading metadata
additional_meta <- read.csv("/stor/work/FRI-BigDataBio/alcohol_addiction_transcriptomics/scratch/rnn384/additional_meta.csv", sep = ",")

#changing the sample names:
alcohol_samples <- grep("Alcohol", additional_meta$Group)

# Step 2: Rename these samples to a1, a2, a3, etc.
additional_meta$Group[alcohol_samples] <- paste0("a", seq_along(alcohol_samples))

meta_data_reference <- human.combined@meta.data

colnames(additional_meta)[colnames(additional_meta) == "Group"] <- "sample_id"

merged_new_meta_data <- meta_data_reference |>
  left_join(additional_meta, by = "sample_id")

merged_new_meta_data$DSMIV
#turn NAs into controls
merged_new_meta_data$DSMIV[is.na(merged_new_meta_data$DSMIV)] <- "Control"

human.combined@meta.data <- merged_new_meta_data
save(human.combined, file="seurat.merged_v4.Rdata")
