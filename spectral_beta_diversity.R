# Source: https://doi.org/10.1016/j.rse.2024.114507
# FUNCTION: SPECTRAL BETA DIVERSITY - SPATIAL ----

calc.beta.spectral.spatial <- function(composition.years){
  
  
  # Create empty dataframe for assimilating dissimilarity matrix outputs as dataframes
  spectra.output.df <- data.frame()
  
  # Create empty lists for assimiliating dissimilarity matrices as matrix
  spectra.output.list.euc <- list()
  spectra.output.list.man <- list()
  spectra.output.list.sam <- list()
  
  
  # Run spectral euclidean distance calculations as a loop
  for (i in composition.years){
    
    
    # Filter to get correct year
    spectra.input.1 <- filter(spectra, Year == i) %>% 
      mutate(Wavelength = round(Wavelength, digits = 1)) # Round to 1.d.p. as some off by 0.0001 in 2020
    
    # Create an NIR key to know what dissimilarity values to put as NA (NIR < 0.2 likely shaded)
    spectra.key.NIR <- spectra.input.1 %>% 
      dplyr::select(PLOT, mean_NIR) %>% 
      mutate(mean_NIR = round(mean_NIR, digits = 5)) %>% # Strange error in plot 201
      arrange(PLOT) %>%
      distinct() %>% 
      mutate(mean_NIR = ifelse(mean_NIR >= 0.2, "Not_Shaded", "Shaded"))
    
    # Create an NDVI key to know what dissimilarity values to put as NA (NDVI >= 0.2 retain)
    spectra.key.NDVI <- spectra.input.1 %>% 
      dplyr::select(PLOT, NDVI_broad) %>% 
      mutate(NDVI_broad = round(NDVI_broad, digits = 5)) %>% # Strange error in plot 201
      arrange(PLOT) %>%
      distinct() %>% 
      mutate(NDVI_broad = ifelse(NDVI_broad >= 0.2, "NDVI_keep", "NDVI_remove"))
    
    # Create an input dataframe for running the euclidean distance calculations
    spectra.input.2 <- spectra.input.1 %>% 
      dplyr::select(PLOT, Band, smooth_Reflectance) %>% 
      pivot_wider(names_from = "Band", values_from = "smooth_Reflectance")
    
    # Create a plot key for rejoining plot information below
    spectra.key.plot <- spectra.input.2 %>% 
      dplyr::select(PLOT) %>% 
      mutate(ID = as.numeric(row.names(.)))
    
    # Remove plot column for calculations
    spectra.input.3 <- dplyr::select(spectra.input.2, -PLOT)
    
    # Create vector of wavelengths
    spectra.wavelengths <- sort(unique(spectra.input.1$Wavelength))
    
    # Turn into class object 'speclib' (DON'T LOAD "hsdar" PACKAGE (confuses tidyverse)
    spectra.speclib <- hsdar::speclib(as.matrix(spectra.input.3), spectra.wavelengths)
    
    # Calculate three types of distance matrix using different methods
    spectra.distance.euc.1 <- hsdar::dist.speclib(spectra.speclib, method = "euclidean") # Euclidean distance
    spectra.distance.man.1 <- hsdar::dist.speclib(spectra.speclib, method = "manhattan") # Manhattan distance
    spectra.distance.sam.1 <- hsdar::dist.speclib(spectra.speclib, method = "sam") # Spectral Angle Metic (SAM)
    
    # Append matrices to lists for outputting
    spectra.output.list.euc[[paste0(i)]] <- spectra.distance.euc.1
    spectra.output.list.man[[paste0(i)]] <- spectra.distance.man.1
    spectra.output.list.sam[[paste0(i)]] <- spectra.distance.sam.1

    # Convert matrices to dataframes
    spectra.distance.euc.2 <- melt(as.matrix(spectra.distance.euc.1), varnames = c("row", "col")) %>% rename(Euclidean = value)
    spectra.distance.man.2 <- melt(as.matrix(spectra.distance.man.1), varnames = c("row", "col")) %>% rename(Manhattan = value)
    spectra.distance.sam.2 <- melt(as.matrix(spectra.distance.sam.1), varnames = c("row", "col")) %>% rename(SAM = value)
    
    # Join dataframes together
    spectra.distance.3 <- spectra.distance.euc.2 %>% 
      left_join(., spectra.distance.man.2, by = c("row" = "row", "col" = "col")) %>% 
      left_join(., spectra.distance.sam.2, by = c("row" = "row", "col" = "col"))

    # Remove duplicate rows in dataframe (from conversion to dataframe from distance matrix)
    spectra.distance.4 <- spectra.distance.3 %>% 
      mutate(Temp_ID = row.names(.)) %>% 
      relocate(Temp_ID, .before = ) %>% 
      pivot_longer(names_to = "PlotIDType", values_to = "PlotID", cols = c("row", "col")) %>% 
      group_by(Temp_ID) %>%
      mutate(minPlotID = min(PlotID), maxPlotID = max(PlotID)) %>% 
      ungroup() %>% 
      distinct(minPlotID, maxPlotID, .keep_all = TRUE) %>% 
      filter(minPlotID != maxPlotID) %>% 
      dplyr::select(minPlotID, maxPlotID, Euclidean, Manhattan, SAM)
    
    # Rename to actual plot numbers and prepare output for combination with other dataframes
    spectra.distance.5 <- left_join(spectra.distance.4, spectra.key.plot, by = c("minPlotID" = "ID")) %>% 
      rename(PLOT_1 = PLOT) %>% 
      left_join(., spectra.key.plot, by = c("maxPlotID" = "ID")) %>%
      rename(PLOT_2 = PLOT) %>% 
      dplyr::select(-c(minPlotID, maxPlotID)) %>% 
      pivot_longer(names_to = "Method", values_to = "Dissimilarity", cols = c("Euclidean", "Manhattan", "SAM")) %>% 
      mutate(Year = i,
             Type = "Spectral") %>% 
      relocate(Year, Type, Method, PLOT_1, PLOT_2, .before = ) %>% 
      arrange(Method, PLOT_1, PLOT_2)
    
    # Replace shaded NIR values (< 0.2) with NAs
    spectra.distance.6 <- left_join(spectra.distance.5, spectra.key.NIR, by = c("PLOT_1" = "PLOT")) %>% 
      rename(Shaded_1 = mean_NIR) %>% 
      left_join(., spectra.key.NIR, by = c("PLOT_2" = "PLOT")) %>% 
      rename(Shaded_2 = mean_NIR) %>%
      mutate(Dissimilarity = ifelse(Shaded_1 == "Shaded" | Shaded_2 == "Shaded", NA, Dissimilarity)) %>% 
      dplyr::select(-c(Shaded_1, Shaded_2))
    
    # If using NDVI threshold (NDVI >= 0.2), replace plots that don't meet it with NAs
    spectra.distance.7 <- left_join(spectra.distance.6, spectra.key.NDVI, by = c("PLOT_1" = "PLOT")) %>% 
      rename(NDVI_1 = NDVI_broad) %>% 
      left_join(., spectra.key.NDVI, by = c("PLOT_2" = "PLOT")) %>% 
      rename(NDVI_2 = NDVI_broad) %>%
      mutate(Dissimilarity = ifelse(NDVI_1 == "NDVI_remove" | NDVI_2 == "NDVI_remove", NA, Dissimilarity)) %>% 
      dplyr::select(-c(NDVI_1, NDVI_2))

    # Join dataframe output to overall output dataframe
    spectra.output.df <- rbind(spectra.output.df, spectra.distance.7)
    
    # Remove intermediate objects
    rm(spectra.input.1, spectra.input.2, spectra.input.3, spectra.key.NIR, spectra.key.plot,
       spectra.speclib, spectra.distance.euc.1, spectra.distance.euc.2, spectra.distance.man.1,
       spectra.distance.man.2, spectra.distance.sam.1, spectra.distance.sam.2, spectra.distance.3,
       spectra.distance.4, spectra.distance.5, spectra.distance.6, spectra.distance.7)


  } # End of years loop
