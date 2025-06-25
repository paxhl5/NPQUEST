```{r}
# Load necessary library
library(dplyr)
library(shiny)

# Define the NPQUEST function
NPQUEST <- function(node_table, edge_table) {
  # Step 1: Split strings into individual objects
  node_table$ATTRIBUTE_SampleStrain <- strsplit(node_table$ATTRIBUTE_SampleStrain, split = ",")
  node_table$ATTRIBUTE_MediaType <- strsplit(node_table$ATTRIBUTE_MediaType, split = ",")
  node_table$ATTRIBUTE_SampleExtractionMethod <- strsplit(node_table$ATTRIBUTE_SampleExtractionMethod, split = ",")
  node_table$ATTRIBUTE_GrowthTimeInDays <- strsplit(node_table$ATTRIBUTE_GrowthTimeInDays, split = ",")
  node_table$UniqueFileSources <- strsplit(node_table$UniqueFileSources, split = ",")
  
  # Step 1: Filter node_table
  node_table_filtered <- node_table %>%
    filter(RTMean < 12.0, lengths(ATTRIBUTE_SampleStrain) == 1 |
    !grepl("CPCC95", ATTRIBUTE_SampleStrain)  , Compound_Name == "", parent.mass < 1500) %>%
    select(!starts_with("GNPSGROUP"))
 
  # Initialize lists to store the scores for each node
  node_names <- character()
  node_media_scores <- numeric()
  node_solvent_scores <- numeric()
  good_connection_scores <- numeric()
  node_intensity_scores <- numeric()
  node_files_scores <- numeric()
  node_spectra_scores <- numeric()
  final_scores <- numeric()
  
  # Step 2: Iterate through each node in the filtered node table
  for (i in 1:nrow(node_table_filtered)) {
    # Calculate the scores
    node_media_score <- lengths(node_table_filtered$ATTRIBUTE_MediaType[i]) * 0.1
    node_solvent_score <- lengths(node_table_filtered$ATTRIBUTE_SampleExtractionMethod[i]) * 0.2
    node_time_score <- lengths(node_table_filtered$ATTRIBUTE_GrowthTimeInDays[i]) * 0.2
    node_files_score <- lengths(node_table_filtered$UniqueFileSources[i]) * 0.05
    
    # Calculate the intensity score
    node_table_intensities <- node_table_filtered[, grep("Peak.area", colnames(node_table_filtered))]
    node_table_intensities$MaxIntensity <- apply(node_table_intensities, 1, max, na.rm = TRUE)
    node_intensity_score <- node_table_intensities$MaxIntensity[i]/max(node_table_intensities$MaxIntensity)
    #node_spectra_score <- node_table_filtered$number.of.spectra[i] * 0.02
    # Initialize the good connections counter
    good_connections <- 0
    
    # Step 3: Iterate through the edge table to find good connections
    for (k in 1:nrow(edge_table)) {
      # Check if node2 in the edge table matches the current node
      if (edge_table$node2[k] == node_table_filtered$name[i]) {
        # Find the index of node1 in the filtered node table
        node1_name <- edge_table$node1[k]
        y <- which(node_table_filtered$name == node1_name)
        
        # Check if the ATTRIBUTE_SampleStrain matches
        if (length(y) > 0 && identical(node_table_filtered$ATTRIBUTE_SampleStrain[[y]], node_table_filtered$ATTRIBUTE_SampleStrain[[i]])) {
          good_connections <- good_connections + 1
        } else {
          good_connections <- good_connections - 2
        }
      }
      
      # Check if node1 in the edge table matches the current node
      if (edge_table$node1[k] == node_table_filtered$name[i]) {
        # Find the index of node2 in the filtered node table
        node2_name <- edge_table$node2[k]
        y <- which(node_table_filtered$name == node2_name)
        
        # Check if the ATTRIBUTE_SampleStrain matches
        if (length(y) > 0 && identical(node_table_filtered$ATTRIBUTE_SampleStrain[[y]], node_table_filtered$ATTRIBUTE_SampleStrain[[i]])) {
          good_connections <- good_connections + 1
        } else {
          good_connections <- good_connections - 2
        }
      }
    }
    
    # Calculate good_connections_score
    if (good_connections < 0) {
      good_connections_score <- -0.3
    } else if (good_connections > 5) {
      good_connections_score <- 0.80
    } else {
      good_connections_score <- good_connections * 0.16
    }
    
    # Calculate the final score
    final_score <- node_media_score + node_time_score + node_files_score + node_intensity_score + node_solvent_score
    
    # Append the results for this node
    node_names <- c(node_names, node_table_filtered$name[i])
    node_media_scores <- c(node_media_scores, node_media_score)
    node_solvent_scores <- c(node_solvent_scores, node_solvent_score)
    good_connection_scores <- c(good_connection_scores, good_connections_score)
    node_intensity_scores <- c(node_intensity_scores, node_intensity_score)
    node_files_scores <- c(node_files_scores, node_files_score)
    #node_spectra_scores <- c(node_spectra_scores, node_spectra_score)
    final_scores <- c(final_scores, final_score)
  }
  
  # Add the results to the filtered node table
  node_table_filtered$node_media_score <- node_media_scores
  node_table_filtered$node_solvent_score <- node_solvent_scores
  node_table_filtered$good_connection_score <- good_connection_scores
  node_table_filtered$intensityscore <- node_intensity_scores
  node_table_filtered$filescore <- node_files_scores
  #node_table_filtered$spectrascore <- node_spectra_scores
  node_table_filtered$final_score <- final_scores/max(final_scores)
  
  return(node_table_filtered)
}

# Apply the function to the provided data
result_table <- NPQUEST(node_table, edge_table)

# Display the result
head(result_table)
 
```