# Source: https://doi.org/10.1016/j.rse.2024.114507
# FUNCTION: FUNCTIONAL BETA DIVERSITY - SPATIAL ----

  # Built on the work of Ricotta and Pavoine (2022) (https://www.sciencedirect.com/science/article/pii/S0304380022000084#sec0008)

calc.beta.functional.spatial <- function(composition.years){
  
  
  # Import functional dissimilarity functions
  source("scripts/EX1_adiv_functions.txt")
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  functional.output.df <- data.frame()
  
  # Create empty list for assimiliating dissimilarity matrix as matrix
  functional.output.list <- list()
  
  
  # Run loop to calculate for each year
  for (i in composition.years){
    
    # Filter saddle traits to correct year
    composition.traits.year <- filter(composition.traits, YEAR == i)
    
    # Modify species name so unique for each plot (allows for intraspecific trait variation)
    composition.traits.sp <- composition.traits.year %>% 
      mutate(SPECIES_unique = paste0(SPECIES, "_", PlotYear)) %>% 
      relocate(SPECIES_unique, .after = SPECIES)
    
    # Transform abundance data so plot is the row name and columns are species
    saddle.abundance <- composition.traits.sp %>% 
      dplyr::select(PLOT, SPECIES_unique, RelativeCover) %>% 
      pivot_wider(names_from = "SPECIES_unique", values_from = "RelativeCover") %>% 
      column_to_rownames(var = "PLOT") %>% # Make plot name the row names
      dplyr::select(order(colnames(.))) %>% # Order species columns alphabetically
      mutate_all(~replace(., is.na(.), 0)) # Replace all NAs with 0
    
    # Transform trait data so have a trait per species
    saddle.traits <- composition.traits.sp %>% 
      dplyr::select(SPECIES_unique, Height, SLA, Chlorophyll, D13C, D15N, LDMC, LeafN, LeafC) %>% # Don't iunclude leaf area as SLA derived from it
      distinct() %>% 
      arrange(SPECIES_unique) %>% 
      column_to_rownames(var = "SPECIES_unique")
    
    # Calculate functional distance between species using Euclidean distance between each species' traits
    saddle.fdis <- dist(scale(saddle.traits)) # Standardized to zero mean and unit standard deviation
    
    # Scale to the unit range by dividing distances by the maximum value in the distance matrix
    saddle.fdis.scaled <- saddle.fdis / max(saddle.fdis)
    
    # Run the functional dissimilarity function from Ricotta and Pavoine (2022) (https://www.sciencedirect.com/science/article/pii/S0304380022000084#sec0008)
    functional.ds.1 <- fundisparam(comm = saddle.abundance,
                                   dis = saddle.fdis.scaled,
                                   method = "D",
                                   abundance = "relative",
                                   alpha = 2,
                                   tol = 1e-8)
    
    # Append matrix to list for outputting
    functional.output.list[[paste0(i)]] <- functional.ds.1
    
    # Convert to dataframe
    functional.ds.2 <- melt(as.matrix(functional.ds.1), varnames = c("row", "col")) %>% 
      rename(Dissimilarity = value)
    
    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    functional.ds.3 <- functional.ds.2 %>% 
      mutate(Temp_ID = row.names(.)) %>% 
      relocate(Temp_ID, .before = ) %>% 
      pivot_longer(names_to = "PlotIDType", values_to = "PlotID", cols = c("row", "col")) %>% 
      group_by(Temp_ID) %>%
      mutate(minPlotID = min(PlotID), maxPlotID = max(PlotID)) %>% 
      ungroup() %>% 
      distinct(minPlotID, maxPlotID, .keep_all = TRUE) %>% 
      filter(minPlotID != maxPlotID) %>% 
      dplyr::select(minPlotID, maxPlotID, Dissimilarity)
    
    # Prepare output for combination with other dataframes
    functional.ds.4 <- functional.ds.3 %>% 
      rename(PLOT_1 = minPlotID, PLOT_2 = maxPlotID) %>% 
      mutate(Year = i,
             Type = "Functional",
             Method = "FDis") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(PLOT_1, PLOT_2)
    
    # Join dataframe output to overall output dataframe
    functional.output.df <- rbind(functional.output.df, functional.ds.4)
    
    # Remove intermediate objects
    rm(composition.traits.year, saddle.abundance, saddle.traits, saddle.fdis, saddle.fdis.scaled,
       functional.ds.1, functional.ds.2, functional.ds.3, functional.ds.4)
    
    
  } # End loop
  
  
  # Export dataframe of functional dissimilarity across all years
  write.csv(functional.output.df, file = paste0("outputs/output_beta_functional_spatial", filepath.brightness,
                                                filepath.smoothing, filepath.37, filepath.top.hits, ".csv"), row.names = FALSE)
  
  # Export list of functional dissimilarity across all years
  save(functional.output.list, file = paste0("outputs/output_beta_functional_spatial", filepath.brightness,
                                             filepath.smoothing, filepath.37, filepath.top.hits, ".RData"))
  
  # Remove intermediate objects
  rm(functional.output.list, functional.output.df)
  
  
} # End of function
