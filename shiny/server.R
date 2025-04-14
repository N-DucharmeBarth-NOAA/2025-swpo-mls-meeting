
# Function to create comparative composition plots (length or weight) across multiple models
plot_composition_comparison = function(model_ids, model_stem, 
                                    comp_type = "length", # Can be "length" or "weight"
                                    show_fit = TRUE, 
                                    n_col = 2, 
                                    free_y_scale = TRUE,  # Option to have free y scales
                                    model_labels = NULL, 
                                    custom_colors = NULL) {
  
  # Required packages
  required_packages <- c("data.table", "ggplot2", "viridis")
  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is needed for this function to work. Please install it."))
    }
  }
  
  # Set comp_type to lowercase for consistent handling
  comp_type = tolower(comp_type)
  
  # Validate comp_type
  if(!comp_type %in% c("length", "weight")) {
    stop("comp_type must be either 'length' or 'weight'")
  }
  
  # Set file names and plot titles based on comp_type
  if(comp_type == "length") {
    file_name = "comp_len.csv"
    kind_filter = "LEN"
    x_label = "Length bin"
    y_label = "Proportion"
    plot_title = "Length Composition Comparison"
  } else { # weight
    file_name = "comp_size.csv"
    kind_filter = "SIZE"
    x_label = "Weight bin"
    y_label = "Proportion"
    plot_title = "Weight Composition Comparison"
  }
  
  # Validate inputs
  if(length(model_ids) < 1) {
    stop("Please provide at least one model ID")
  }
  
  # Handle model_stem as vector or single path
  if(length(model_stem) == 1) {
    model_stem = rep(model_stem, length(model_ids))
  } else if(length(model_stem) != length(model_ids)) {
    stop("If model_stem is a vector, it must have the same length as model_ids")
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  model_labels = sort(model_labels)
  
  # Create vector of model specifications with labels
  selected_models = data.frame(
    id = model_ids,
    path = model_stem,
    label = model_labels,
    stringsAsFactors = FALSE
  )
  
  # Check if composition file exists for all models
  file_exists = sapply(1:nrow(selected_models), function(i) {
    file.exists(file.path(selected_models$path[i], selected_models$id[i], file_name))
  })
  
  if(sum(file_exists) == 0) {
    stop(paste0("The file ", file_name, " does not exist for any of the models you selected."))
  } else if(sum(file_exists) < length(file_exists)) {
    warning(paste0("The file ", file_name, " does not exist for some models. These models will be excluded."))
    selected_models = selected_models[file_exists, ]
  }
  
  # Read and compile composition data
  comp_data_list = lapply(1:nrow(selected_models), function(i) {
    # Read composition data
    comp_file = file.path(selected_models$path[i], selected_models$id[i], file_name)
    comp_dt = fread(comp_file)
    
    # Filter by Kind
    comp_dt = comp_dt[Kind == kind_filter]
    
    # Add model identifiers
    comp_dt[, model_id := selected_models$id[i]]
    comp_dt[, model_label := selected_models$label[i]]
    
    # Round observed values for consistent ID generation
    comp_dt[, Obs := round(Obs, digits=6)]
    
    return(comp_dt)
  })
  
  # Combine all composition data
  plot_dt = rbindlist(comp_data_list, fill = TRUE)
  
  # Check if we have any data after filtering
  if(nrow(plot_dt) == 0) {
    stop(paste0("No ", comp_type, " composition data found with Kind = ", kind_filter))
  }
  
  # Calculate bin widths for each fleet
  bin_info = plot_dt[, .(
    Bin = sort(unique(Bin))
  ), by = .(Fleet_name)]
  
  # Calculate bin widths as differences between consecutive bins
  # For the last bin, use the width of the previous bin
  bin_info[, bin_width := c(diff(Bin), tail(diff(Bin), 1)), by = .(Fleet_name)]
  
  # Calculate bin centers and edges for plotting
  bin_info[, bin_center := Bin + bin_width/2]
  
  # Merge bin widths back to the main data table
  plot_dt = merge(plot_dt, bin_info, by = c("Fleet_name", "Bin"))
  
  # Check consistency of observed values across models for the same fleet/bin
  consistency_check = plot_dt[, .(
    obs_values = list(unique(Obs)),
    models = list(unique(model_label))
  ), by = .(Fleet_name, Bin)]
  
  # Check if observed values are consistent
  obs_inconsistent = consistency_check[sapply(obs_values, length) > 1]
  if(nrow(obs_inconsistent) > 0) {
    error_msg = paste0(
      "Observed ", comp_type, " composition values differ across models for the same fleet and bin:\n",
      paste(
        apply(obs_inconsistent[, .(Fleet_name, Bin)], 1, function(row) {
          return(paste0("  - Fleet: ", row[1], ", Bin: ", row[2]))
        }),
        collapse = "\n"
      )
    )
    stop(error_msg)
  }
  
  # Create a version of the data with just the first model's observations
  # This ensures we only plot the observations once
  obs_dt = plot_dt[, .SD[1], by = .(Fleet_name, Bin)]
  
  # Check if we have data to plot
  if(nrow(plot_dt) == 0) {
    stop("No data available after applying filters")
  }
  
  # Set model_label as factor with correct order
  plot_dt[, model_label := factor(model_label, levels = model_labels)]
  
  # Use Fleet_name as facet label
  plot_dt[, facet_label := Fleet_name]
  obs_dt[, facet_label := Fleet_name]
  
  # Combine data for plotting
  combined_obs = obs_dt[, .(Fleet_name, Bin, bin_width, bin_center, Obs, facet_label)]
  combined_exp = plot_dt[, .(Fleet_name, Bin, bin_width, bin_center, Exp, model_label, facet_label)]
  
  # Create a combined plot with facets
  combined_plot = ggplot() +
    # Add bars for observed data using geom_rect
    geom_rect(data = combined_obs, 
              aes(xmin = Bin, xmax = Bin + bin_width, 
                  ymin = 0, ymax = Obs),
              fill = "gray80", color = "gray40", alpha = 0.7) +
    # Add fitted lines if requested
    {if(show_fit) geom_line(data = combined_exp, 
                           aes(x = bin_center, y = Exp, color = model_label, group = model_label),
                           linewidth = 1)} +
    # Facet by fleet, with option for free y scales
    facet_wrap(~ facet_label, ncol = n_col, scales = if(free_y_scale) "free_y" else "fixed") +
    # Labels
    labs(x = x_label, y = y_label, title = plot_title)
  
  # Apply color scales - either custom colors or viridis palette
  color_count = length(unique(plot_dt$model_label))
  color_legend = "Model"
  
  if(!is.null(custom_colors)) {
    # Validate that enough colors are provided
    if(length(custom_colors) < color_count) {
      warning("Not enough custom colors provided. Falling back to viridis palette.")
      combined_plot = combined_plot + 
        viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                   direction = 1, option = "H", discrete = TRUE)
    } else {
      # Use custom colors
      combined_plot = combined_plot + scale_color_manual(color_legend, values = custom_colors)
    }
  } else {
    # Use default viridis palette
    combined_plot = combined_plot + 
      viridis::scale_color_viridis(color_legend, begin = 0.1, end = 0.8, 
                                 direction = 1, option = "H", discrete = TRUE)
  }
  
  # Apply consistent theme elements
  combined_plot = combined_plot + theme(
    panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
    panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
    panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
    strip.background = element_rect(fill = "white"),
    legend.key = element_rect(fill = "white")
  )
  
    return(combined_plot)

}


# Function to create comparative length-based selectivity plots
plot_model_comparison_selex = function(model_ids, model_stem, 
                                      n_col = 2, model_labels = NULL, custom_colors = NULL) {
  
  # Validate inputs
  if(length(model_ids) < 1) {
    stop("Please provide at least one model ID")
  }

  # Handle model_stem as vector or single path
  if(length(model_stem) == 1) {
    model_stem = rep(model_stem, length(model_ids))
  } else if(length(model_stem) != length(model_ids)) {
    stop("If model_stem is a vector, it must have the same length as model_ids")
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  model_labels = sort(model_labels)
  
  # Create vector of model specifications with labels
  selected_models = data.frame(
    id = model_ids,
    path = model_stem,
    label = model_labels,
    stringsAsFactors = FALSE
  )
  
  # Read and process selectivity data from each model
  selex_list = lapply(1:nrow(selected_models), function(i) {
    # Read selectivity data
    selex_file = file.path(selected_models$path[i], selected_models$id[i], "selex_l.csv")
    if (!file.exists(selex_file)) {
      warning(paste("File not found:", selex_file))
      return(NULL)
    }
    
    dt = fread(selex_file)
    dt[, model_id := selected_models$id[i]]
    dt[, model_label := selected_models$label[i]]
    return(dt)
  })
  
  # Combine all selectivity data
  selex_data = rbindlist(selex_list, use.names = TRUE, fill = TRUE)
  
  # Check if we have data
  if(nrow(selex_data) == 0) {
    stop("No selectivity data found for the specified models")
  }
  
  # Read fleet summary data for each model to get fleet names
  fleet_list = lapply(1:nrow(selected_models), function(i) {
    # Read fleet summary data
    fleet_file = file.path(selected_models$path[i], selected_models$id[i], "fleet_summary.csv")
    if (!file.exists(fleet_file)) {
      warning(paste("Fleet summary file not found:", fleet_file))
      return(NULL)
    }
    
    dt = fread(fleet_file)
    dt[, model_id := selected_models$id[i]]
    return(dt)
  })
  
  # Combine all fleet data
  fleet_data = rbindlist(fleet_list, use.names = TRUE, fill = TRUE)
  
  # Check if we have fleet data
  if(nrow(fleet_data) == 0) {
    warning("No fleet summary data found, using original Fleet_name")
  } else {
    # Map fleet names from fleet_summary to selex_data
    selex_data = merge(
      selex_data,
      fleet_data[, .(fleet, fleetname, model_id)],
      by.x = c("Fleet", "model_id"),
      by.y = c("fleet", "model_id"),
      all.x = TRUE
    )
    
    # Use fleetname from fleet_summary if available, otherwise keep original Fleet_name
    selex_data[!is.na(fleetname), Fleet_name := fleetname]
  }
  
  # Check if fleet names are consistent across models
  fleet_consistency = selex_data[, .(
    unique_names = length(unique(Fleet_name))
  ), by = .(Fleet)]
  
  inconsistent_fleets = fleet_consistency[unique_names > 1, Fleet]
  if(length(inconsistent_fleets) > 0) {
    warning(paste("Inconsistent fleet names found for fleet(s):", 
                 paste(inconsistent_fleets, collapse = ", ")))
  }
  
  # Add sex label for better readability
  selex_data[Sex == 1, sex_label := "Female"]
  selex_data[Sex == 2, sex_label := "Male"]
  selex_data[is.na(sex_label), sex_label := "NA"]
  
  # For each Fleet, if there's only one sex, change the label to "Aggregated"
  fleet_sex_counts <- selex_data[, .(sex_count = uniqueN(Sex)), by = .(Fleet)]
  single_sex_fleets <- fleet_sex_counts[sex_count == 1, Fleet]
  
  if(length(single_sex_fleets) > 0) {
    selex_data[Fleet %in% single_sex_fleets, sex_label := "Aggregated"]
  }
  
  # Convert labels to factors to ensure consistent ordering
  selex_data[, model_label := factor(model_label, levels = model_labels)]
  # Sort unique sex labels alphabetically for factor levels
  sex_levels <- sort(unique(selex_data$sex_label))
  selex_data[, sex_label := factor(sex_label, levels = sex_levels)]
  
  # Create plot
  p = ggplot(selex_data, aes(x = variable, y = value, color = model_label, linetype = sex_label)) +
    geom_line(linewidth = 1.2) +
    facet_wrap(~ Fleet_name, scales = "free_x", ncol = n_col) +
    ylim(0, 1) +
    xlab("Length") +
    ylab("Selectivity") +
    labs(color = "Model", linetype = "Sex") +
    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
          panel.grid.major = element_line(color = 'gray70', linetype = "dotted"),
          panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
          strip.background = element_rect(fill = "white"),
          legend.key = element_rect(fill = "white"))
  
  # Apply color scales - either custom colors or viridis palette
  if(!is.null(custom_colors)) {
    # Validate that enough colors are provided
    if(length(custom_colors) < length(unique(selex_data$model_label))) {
      warning("Not enough custom colors provided. Falling back to viridis palette.")
      p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                         direction = 1, option = "H", discrete = TRUE)
    } else {
      # Use custom colors
      p = p + scale_color_manual("Model", values = custom_colors)
    }
  } else {
    # Use default viridis palette
    p = p + viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                                       direction = 1, option = "H", discrete = TRUE)
  }

    return(p)
}

plot_catch_comparison <- function(model_ids, model_stem = "data", 
                                fleets = NULL,        
                                show_fit = TRUE,
                                n_col = 2,
                                free_y_scale = TRUE,  
                                model_labels = NULL) {
  
  # Validate inputs
  if(length(model_ids) < 1) {
    return(NULL)
  }
  
  # Get paths for selected models
  selected_models <- sapply(model_ids, function(x) all_dirs[grep(x, all_dirs, fixed=TRUE)])
  
  # Check if catch.csv exists for all models
  keep_selected <- which(sapply(file.path(model_stem, selected_models, "catch.csv"), file.exists) == TRUE)
  if(length(keep_selected) == 0) {
    warning("The file catch.csv does not exist for any of the files you selected.")
    return(NULL)
  } else {
    selected_models <- selected_models[keep_selected]
    model_ids <- model_ids[keep_selected]
  }
  
  # If model_labels is not provided, use model_ids as labels
  if(is.null(model_labels)) {
    model_labels = model_ids
  } else if(length(model_labels) != length(model_ids)) {
    stop("If model_labels is provided, it must have the same length as model_ids")
  }
  model_labels = sort(model_labels)
  
  # Read and combine catch data
  plot_dt <- rbindlist(lapply(seq_along(selected_models), function(i) {
    model_id <- model_ids[i]
    model_path <- selected_models[i]
    catch_data <- fread(file.path(model_stem, model_path, "catch.csv"))
    
    # Add model identifiers
    catch_data[, id := model_id]
    catch_data[, model_label := model_labels[i]]
    
    # Round values for consistent comparison
    catch_data[, Obs := round(Obs, digits=6)]
    catch_data[, Exp := round(Exp, digits=6)]
    
    return(catch_data)
  }))
  
  # If no fleet names, create them
  if(!"Fleet_name" %in% names(plot_dt)) {
    # Try to read fleet names from fleet_summary.csv if available
    fleet_names <- lapply(seq_along(selected_models), function(i) {
      model_id <- model_ids[i]
      model_path <- selected_models[i]
      fleet_path <- file.path(model_stem, model_path, "fleet_summary.csv")
      
      if(file.exists(fleet_path)) {
        fleet_data <- fread(fleet_path) %>%
                      .[,id:=model_ids[i]] %>%
                      .[, .(id,fleet, fleetname)]
        return(fleet_data)
      } else {
        return(NULL)
      }
    })
    
    # Combine all fleet names or create default names
    if(all(sapply(fleet_names, is.null))) {
      plot_dt[, Fleet_name := paste("Fleet", Fleet)]
    } else {
      fleet_dt <- rbindlist(fleet_names[!sapply(fleet_names, is.null)])
      setnames(fleet_dt, c("fleet", "fleetname"), c("Fleet", "Fleet_name"))
      
      # Merge fleet names with catch data
      plot_dt <- merge(plot_dt, fleet_dt, by = c("id","Fleet"), all.x = TRUE)
      
      # For any missing fleet names, create default
      plot_dt[is.na(Fleet_name), Fleet_name := paste("Fleet", Fleet)]
    }
  }
  
  # Filter by fleet if specified
  if(!is.null(fleets) && length(fleets) > 0) {
    plot_dt <- plot_dt[Fleet %in% fleets]
    
    # Check if we still have data after filtering
    if(nrow(plot_dt) == 0) {
      warning("No data available after filtering for the specified fleets")
      return(NULL)
    }
  }
  
  # Check for consistency in observed values across models
  consistency_check <- plot_dt[, .(
    obs_values = list(unique(Obs)),
    models = list(unique(model_label))
  ), by = .(Fleet, Time)]
  
  # Warn about inconsistencies if found
  obs_inconsistent <- consistency_check[sapply(obs_values, length) > 1]
  if(nrow(obs_inconsistent) > 0) {
    warning("Observed catch values differ across models for some fleet/time combinations.")
  }
  
  # Create a version of the data with just the first model's observations
  # This ensures we only plot the observations once
  obs_dt <- plot_dt[, .SD[1], by = .(Fleet, Fleet_name, Time)]
  
  # Calculate bins for bar width
  unique_times <- sort(unique(plot_dt$Time))
  time_diffs <- diff(unique_times)
  non_zero_diffs <- time_diffs[time_diffs > 0]
  bin_width <- if(length(non_zero_diffs) > 0) min(non_zero_diffs) else 1
  
  # Add bin information for plotting
  plot_dt[, bin_start := Time]
  plot_dt[, bin_end := Time + bin_width]
  plot_dt[, bin_center := Time + bin_width/2]
  
  obs_dt[, bin_start := Time]
  obs_dt[, bin_end := Time + bin_width]
  obs_dt[, bin_center := Time + bin_width/2]
  
  # Set model_label as factor with correct order
  plot_dt[, model_label := factor(model_label, levels = model_labels)]
  
  # Create the plot
  p <- ggplot() +
    # Add bars for observed data
    geom_rect(data = obs_dt, 
              aes(xmin = bin_start, xmax = bin_end, 
                  ymin = 0, ymax = Obs),
              fill = "gray80", color = "gray40", alpha = 0.7) +
    # Add fitted lines if requested
    {if(show_fit) geom_line(data = plot_dt, 
                           aes(x = bin_center, y = Exp, color = model_label, group = model_label),
                           linewidth = 1)} +
    # Facet by fleet
    facet_wrap(~ Fleet_name, ncol = n_col, scales = if(free_y_scale) "free_y" else "fixed") +
    # Labels
    labs(x = "Year", y = "Catch", title = "Observed vs. Fitted Catch")
  
  # Apply consistent theme and color scheme
  p <- p + 
    viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                               direction = 1, option = "H", discrete = TRUE) +
    theme(
      panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
      panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
      panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
      strip.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white")
    )
  
  return(p)
}

create_cpue_comparison_plot <- function(model_ids, model_stem = "data", 
                                      show_se = TRUE, 
                                      show_fit = TRUE, 
                                      use_log_scale = FALSE,
                                      apply_varadj = TRUE,
                                      filter_lambda = TRUE) {
  
  # Validate inputs
  if(length(model_ids) < 1) {
    return(NULL)
  }
  
  # Get paths for selected models
  selected_models <- sapply(model_ids, function(x) all_dirs[grep(x, all_dirs, fixed=TRUE)])
  
  # Check if cpue.csv exists for all models
  keep_selected <- which(sapply(file.path(model_stem, selected_models, "cpue.csv"), file.exists) == TRUE)
  if(length(keep_selected) == 0) {
    warning("The file cpue.csv does not exist for any of the files you selected.")
    return(NULL)
  } else {
    selected_models <- selected_models[keep_selected]
    model_ids <- model_ids[keep_selected]
  }
  
  # Read and combine CPUE data
  plot_dt <- rbindlist(lapply(seq_along(selected_models), function(i) {
    model_id <- model_ids[i]
    model_path <- selected_models[i]
    cpue_data <- fread(file.path(model_stem, model_path, "cpue.csv"))
    
    # Rename columns to be consistent
    setnames(cpue_data, 
             c("Fleet", "Fleet_name", "Time", "Obs", "Exp", "SE"), 
             c("fleet", "fleetname", "time", "obs", "exp", "se"), 
             skip_absent = TRUE)
    
    # Add model identifiers
    cpue_data[, id := model_id]
    cpue_data[, model_label := model_id] # Can be customized if needed
    
    # Round observed values for consistent comparison
    cpue_data[, obs := round(obs, digits=3)]
    cpue_data[, id3 := as.numeric(as.factor(paste0(id, "_", fleetname)))]
    
    return(cpue_data)
  }))
  
  # Read variance adjustment and lambda values
  varlambda_dt <- rbindlist(lapply(seq_along(selected_models), function(i) {
    model_id <- model_ids[i]
    model_path <- selected_models[i]
    
    # Read variance adjustments
    tmp_var <- fread(file.path(model_stem, model_path, "varadj.csv"))
    if(nrow(tmp_var) == 0) {
      tmp_var <- data.table(
        fleet = unique(plot_dt[id == model_id]$fleet),
        id = model_id,
        factor = 1,
        value = 0
      )[, .(id, factor, fleet, value)]
    } else {
      tmp_var <- tmp_var[factor == 1]
      tmp_var2 <- data.table(
        fleet = unique(plot_dt[id == model_id]$fleet),
        id = model_id,
        factor = 1,
        value = 0
      )[, .(id, factor, fleet, value)]
      
      tmp_var2 <- tmp_var2[!(fleet %in% tmp_var$fleet)]
      tmp_var <- rbind(tmp_var, tmp_var2)[order(fleet)]
    }
    
    tmp_var <- tmp_var[, .(id, fleet, value)]
    setnames(tmp_var, "value", "varadj")
    
    # Read lambdas
    tmp_lambda <- fread(file.path(model_stem, model_path, "lambdas.csv"))
    if(nrow(tmp_lambda) == 0) {
      tmp_lambda <- data.table(
        fleet = unique(plot_dt[id == model_id]$fleet),
        id = model_id,
        value = 0
      )[, .(id, fleet, value)]
    } else {
      tmp_lambda <- tmp_lambda[like_comp == 1][, .(id, fleet, value)]
      tmp_lambda2 <- data.table(
        fleet = unique(plot_dt[id == model_id]$fleet),
        id = model_id,
        value = 0
      )[, .(id, fleet, value)]
      
      tmp_lambda2 <- tmp_lambda2[!(fleet %in% tmp_lambda$fleet)]
      tmp_lambda <- rbind(tmp_lambda, tmp_lambda2)[order(fleet)]
    }
    
    setnames(tmp_lambda, "value", "lambda")
    
    return(merge(tmp_var, tmp_lambda, by = c("id", "fleet")))
  }))
  
  # Merge CPUE data with variance adjustment and lambda data
  plot_dt <- merge(plot_dt, varlambda_dt, by = c("id", "fleet"))
  
  # Apply variance adjustment if specified
  if(!apply_varadj) {
    plot_dt[, se := se - varadj]
  }
  
  # Consistency checks for observations
  obs_check <- plot_dt[, .(
    obs_values = list(unique(obs)),
    models = list(unique(id))
  ), by = .(fleetname, time)]
  
  obs_inconsistent <- obs_check[sapply(obs_values, length) > 1]
  if(nrow(obs_inconsistent) > 0) {
    warning(paste0(
      "Observed CPUE values differ across models for the same fleet and time point. ",
      "This may indicate model configuration differences."
    ))
  }
  
  # Consistency checks for standard errors
  se_check <- plot_dt[, .(
    se_values = list(unique(se)),
    models = list(unique(id))
  ), by = .(fleetname, time)]
  
  se_inconsistent <- se_check[sapply(se_values, length) > 1]
  if(nrow(se_inconsistent) > 0 && show_se) {
    warning(paste0(
      "Standard errors differ across models for the same fleet and time point. ",
      "Error bars may not be consistent."
    ))
  }
  
  # Create a version of the data with just the first model's observations
  # This ensures we only plot the observations once
  obs_dt <- plot_dt[, .SD[1], by = .(fleetname, time)]
  
  # Apply log transformation if requested
  if(use_log_scale) {
    plot_dt[, obs := log(obs)]
    plot_dt[, exp := log(exp)]
    plot_dt[, lse := obs - se]
    plot_dt[, use := obs + se]
    
    obs_dt[, obs := log(obs)]
    obs_dt[, lse := obs - se]
    obs_dt[, use := obs + se]
    
    ylab_txt <- "Index (log-scale)"
    yint <- 0
  } else {
    plot_dt[, lse := exp(log(obs) - se)]
    plot_dt[, use := exp(log(obs) + se)]
    
    obs_dt[, lse := exp(log(obs) - se)]
    obs_dt[, use := exp(log(obs) + se)]
    
    ylab_txt <- "Index"
    yint <- 1
  }
  
  # Filter by lambda if requested
  if(filter_lambda) {
    plot_dt <- plot_dt[lambda == 1]
  }
  
  # Check if we still have data to plot after filtering
  if(nrow(plot_dt) == 0) {
    warning("No data available after applying filters")
    return(NULL)
  }
  
    # Set label as factor with correct order
  plot_dt[, label := factor(model_label, levels = sort(model_ids))]
  
  # Create base plot
  p = ggplot(plot_dt[order(id3, time)])  

  # Set y-axis limits for non-log scale
  if(!use_log_scale) {
    p = p + ylim(0, NA)
  }
  
  # Add common elements
  if(uniqueN(plot_dt$fleetname)>=3){
    n_col = 3
  } else {
    n_col = uniqueN(plot_dt$fleetname)
  }
  
  p = p + 
    xlab("Time") +
    ylab(ylab_txt) +
    geom_hline(yintercept = yint, linetype = "dashed", linewidth = 1) +
    facet_wrap(~fleetname, ncol = n_col)
  
  # Add error bars if requested
  if(show_se) {
    p = p + geom_errorbar(data = obs_dt, aes(x = time, ymin = lse, ymax = use), 
                       width = 0.2, linewidth = 0.7, color = "gray60")
  }
  
  # Add observed points
  p = p + geom_point(data = obs_dt, aes(x = time, y = obs), 
                   shape = 21, size = 3, stroke = 0.5, fill = "gray90")
  
  # Add fitted lines if requested
  if(show_fit) {
    p = p + geom_line(aes(x = time, y = exp, color = label), linewidth = 1)
  }
  
  # Apply consistent theme and color scheme
  p <- p + 
    viridis::scale_color_viridis("Model", begin = 0.1, end = 0.8, 
                               direction = 1, option = "H", discrete = TRUE) +
    viridis::scale_fill_viridis("Model", begin = 0.1, end = 0.8, 
                              direction = 1, option = "H", discrete = TRUE) +
    theme(
      panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
      panel.grid.major = element_line(color = 'gray70', linetype = "dotted"), 
      panel.grid.minor = element_line(color = 'gray70', linetype = "dotted"),
      strip.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white")
    )
  
  return(p)
}

server = function(input, output){
  # pixel height for each panel. i.e row height when plotting by species
  height_per_panel = 350

  ref_table_reduced = summary_dt %>%
                as.data.frame(.)

  output$summarytable = DT::renderDataTable({
    summary_df = summary_dt %>%
                 as.data.frame(.,stringsAsFactors=FALSE)
    summary_DT = DT::datatable(summary_df, filter = 'top',rownames=FALSE,
    options = list(scrollX = TRUE, search = list(regex = TRUE, caseInsensitive = FALSE),pageLength = 25))
    return(summary_DT)
  })
  outputOptions(output, "summarytable", suspendWhenHidden = FALSE)

#   filtered_table = reactive({
#     req(input$summarytable_rows_selected)
#     keep_models = c(ref_table_reduced[input$summarytable_rows_selected, ]$id)
#     return(as.data.frame(summary_dt[id%in%keep_models],stringsAsFactors=FALSE))  
#   })

  filtered_id = reactive({
    req(input$summarytable_rows_selected)
    keep_models = c(ref_table_reduced[input$summarytable_rows_selected, ]$id)
    return(keep_models)  
  })

  # define plots
  output$ts_plots = renderPlot({
    input_models = unique(filtered_id())
    input_category = input$category
    if(length(input_models) < 1 | length(input_category) < 1){
      return()
    }
    
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    plot_dt.list = list("SSB"=NA, "F"=NA, "Recruitment"=NA,"Depletion (D)"=NA,"SSB (Unfished)"=NA,"SSB/SSBmsy"=NA,"F/Fmsy"=NA)
    
    if("SSB" %in% input_category){
        plot_dt.list[["SSB"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"ssb.csv")))) %>%
             setnames(.,c("Value","StdDev"),c("value","sd")) %>% 
             .[,cv:=sd/value] %>%
			 .[,l95:=exp(log(value)-2*sqrt(log(cv^2+1)))] %>%
			 .[,u95:=exp(log(value)+2*sqrt(log(cv^2+1)))] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:="SSB"] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("F" %in% input_category){
        plot_dt.list[["F"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"f.csv")))) %>%
             setnames(.,c("Value","StdDev"),c("value","sd")) %>% 
             .[,cv:=sd/value] %>%
			 .[,l95:=exp(log(value)-2*sqrt(log(cv^2+1)))] %>%
			 .[,u95:=exp(log(value)+2*sqrt(log(cv^2+1)))] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:="F"] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("Recruitment" %in% input_category){
        plot_dt.list[["Recruitment"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"rec.csv")))) %>%
             setnames(.,c("Value","StdDev"),c("value","sd")) %>% 
             .[,cv:=sd/value] %>%
			 .[,l95:=exp(log(value)-2*sqrt(log(cv^2+1)))] %>%
			 .[,u95:=exp(log(value)+2*sqrt(log(cv^2+1)))] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:="Recruitment"] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("Depletion (D)" %in% input_category){
        plot_dt.list[["Depletion (D)"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"dyn_dep.csv")))) %>%
             .[Era %in% c("TIME")] %>%
             setnames(.,c("Yr"),c("yr")) %>% 
             .[,value:=SSB/SSB_nofishing] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="Depletion (D)"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("SSB (Unfished)" %in% input_category){
        plot_dt.list[["SSB (Unfished)"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"dyn_dep.csv")))) %>%
             .[Era %in% c("TIME")] %>%
             setnames(.,c("Yr"),c("yr")) %>% 
             .[,value:=SSB_nofishing] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="SSB (Unfished)"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }
    if("SSB/SSBmsy" %in% input_category){
        plot_dt.list[["SSB/SSBmsy"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"kobe.csv")))) %>%
             setnames(.,c("Year"),c("yr")) %>% 
             .[,value:=B.Bmsy] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="SSB/SSBmsy"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }
    if("F/Fmsy" %in% input_category){
        plot_dt.list[["F/Fmsy"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"kobe.csv")))) %>%
             setnames(.,c("Year"),c("yr")) %>% 
             .[,value:=F.Fmsy] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="F/Fmsy"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    plot_dt = rbindlist(plot_dt.list[which(sapply(plot_dt.list,length)>1)]) %>% .[order(id,type,yr)]

    if(nrow(plot_dt) == 0){
      return()
    }
    if(length(input_category)>2)
    {
      n_col = 3
    } else {
      n_col = length(input_category)
    }
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Year") +
			facet_wrap(~type,scales="free_y",ncol=n_col)

      p = p + ylab("Metric") +
              geom_hline(yintercept=0)
      if(input$se)
      {
        p = p + geom_ribbon(aes(x=yr,ymin=l95,ymax=u95,group=id,fill=id),alpha=0.25)
      }
	  p = p + geom_path(aes(x=yr,y=value,group=id,color=id),linewidth=1.5)
    
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
    if(length(input$category)>2)
    {
      n_col = 3
    } else {
      n_col = length(input$category)
    }
    return(max(height_per_panel*1.5, (height_per_panel * ceiling(length(input$category) / n_col))))
  })
  # growth_plot
  output$growth_plot = renderPlot({
    input_models = unique(filtered_id())
    input_growth_var = input$growth_var
    input_bio_sex = input$bio_sex
    if(length(input_models) < 1|length(input_bio_sex)<1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"bio.csv")))) %>%
             .[,.(id,Sex,A1,A2,L_a_A1,L_a_A2,K,CVmin,CVmax)]
    max_age_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"selex_a.csv"))[,.(max_age=max(as.numeric(variable))),by=id]))
    ages_dt.list = as.list(rep(NA,length(selected_models)))
    for(i in 1:length(selected_models)){
        ages_dt.list[[i]] = as.data.table(expand.grid(id=input_models[i],age=0:max_age_dt[id==input_models[i]]$max_age))
    }
    age_dt = rbindlist(ages_dt.list)

    plot_dt = merge(plot_dt,age_dt,by="id",allow.cartesian=TRUE) %>%
              .[,pred_len := L_a_A1 + (L_a_A2-L_a_A1)*((1-exp(-K*(age-A1)))/(1-exp(-K*(A2-A1))))] %>%
              .[age==A1,sd:=pred_len*CVmin] %>%
              .[age<A2,sd:=pred_len*(CVmin+((pred_len-L_a_A1)/(L_a_A2-L_a_A1))*(CVmax-CVmin))] %>%
              .[age==A2,sd:=pred_len*CVmax] %>%
              .[,u95:=pred_len+1.96*sd] %>%
              .[,l95:=pred_len-1.96*sd] %>%
              .[,.(id,Sex,age,pred_len,u95,l95)] %>%
              .[,Sex:=factor(Sex,levels=as.character(2:1),labels=c("Male","Female"))] %>%
              .[Sex %in% input_bio_sex]

    if(nrow(plot_dt) == 0){
      return()
    }
    
    n_col = length(input_bio_sex)
    
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Age (years)") +
			facet_wrap(~Sex,ncol=n_col)

      p = p + ylab("Length (cm)") +
              geom_hline(yintercept=0)
      if(input$growth_var)
      {
        p = p + geom_ribbon(aes(x=age,ymin=l95,ymax=u95,group=id,fill=id),alpha=0.25)
      }
	  p = p + geom_path(aes(x=age,y=pred_len,group=id,color=id),linewidth=1.5)
    
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
        return(height_per_panel*1.5)
  })
# m_plot
  output$m_plot = renderPlot({
    input_models = unique(filtered_id())
    input_bio_sex = input$bio_sex
    if(length(input_models) < 1|length(input_bio_sex)<1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"bio.csv")))) %>%
             .[,.(id,Sex,M_age0,M_nages)] %>%
             melt(.,id.vars=c("id","Sex")) %>%
             setnames(.,c("variable","value"),c("age_cat","m"))

    max_age_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"selex_a.csv"))[,.(M_age0=min(as.numeric(variable)),M_nages=max(as.numeric(variable))),by=id])) %>%
    melt(.,id.vars=c("id")) %>%
    setnames(.,c("variable","value"),c("age_cat","age"))


    plot_dt = merge(plot_dt,max_age_dt,by=c("id","age_cat"),allow.cartesian=TRUE) %>%
              .[,Sex:=factor(Sex,levels=as.character(2:1),labels=c("Male","Female"))] %>%
              .[Sex %in% input_bio_sex]

    if(nrow(plot_dt) == 0){
      return()
    }
    
    n_col = length(input_bio_sex)
    
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Age (years)") +
			facet_wrap(~Sex,ncol=n_col)

      p = p + ylab("Natural mortality (M)") +
              geom_hline(yintercept=0)
	  p = p + geom_path(aes(x=age,y=m,group=id,color=id),linewidth=1.5) +
            geom_point(aes(x=age,y=m,group=id,fill=id),shape=21,size=5)
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
        return(height_per_panel*1.5)
  })
  # lw_plot
  output$lw_plot = renderPlot({
    input_models = unique(filtered_id())
    input_bio_sex = input$bio_sex
    if(length(input_models) < 1|length(input_bio_sex)<1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"bio.csv")))) %>%
             .[,.(id,Sex,WtLen1,WtLen2)] 

    max_len_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"selex_l.csv"))[,.(max_len=max(as.numeric(variable))),by=id])) 
    len_dt.list = as.list(rep(NA,length(selected_models)))
    for(i in 1:length(selected_models)){
        len_dt.list[[i]] = as.data.table(expand.grid(id=input_models[i],len=0:ceiling(max_len_dt[id==input_models[i]]$max_len)))
    }
    len_dt = rbindlist(len_dt.list)

    plot_dt = merge(plot_dt,len_dt,by="id",allow.cartesian=TRUE) %>%
              .[,pred_weight := WtLen1*len^WtLen2] %>%
              .[,.(id,Sex,len,pred_weight)] %>%
              .[,Sex:=factor(Sex,levels=as.character(2:1),labels=c("Male","Female"))] %>%
              .[Sex %in% input_bio_sex]

    if(nrow(plot_dt) == 0){
      return()
    }
    
    n_col = length(input_bio_sex)
    
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Length (cm)") +
			facet_wrap(~Sex,ncol=n_col)

      p = p + ylab("Weight (kg)") +
              geom_hline(yintercept=0)
	  p = p + geom_path(aes(x=len,y=pred_weight,group=id,color=id),linewidth=1.5)
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
        return(height_per_panel*1.5)
  })

  # catch_total_plot
  output$catch_total_plot = renderPlot({
    input_models = unique(filtered_id())
    input_units = input$catch_units
    input_aggregate = input$catch_by_fleet
    input_fill = input$catch_fill
    input_rescale = input$catch_rescale
    input_chart = input$catch_chart
    if(length(input_models) < 1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    # check to make sure catch.csv exists for all the selected models and update if needed
    keep_selected = which(sapply(file.path(model_stem,selected_models,"catch.csv"),file.exists)==TRUE)
    if(length(keep_selected)==0){
      return(warning("The file catch.csv does not exist for any of the files you selected."))
    } else {
      selected_models = selected_models[keep_selected]
    }
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"catch.csv")))) %>%
             .[,.(id,Fleet,Time,Obs,Exp,sel_bio,sel_num)] %>%
             setnames(.,c("Fleet","Time","Obs","Exp"),c("fleet","time","obs","exp"))
    u_id = unique(plot_dt$id)
    # read in fleet_summary to get names and units
    fleet_dt.list = as.list(rep(NA,length(selected_models)))
    for(i in 1:length(selected_models)){
      fleet_dt.list[[i]] = fread(file.path(model_stem,selected_models[i],"fleet_summary.csv")) %>%
                           .[,id:=u_id[i]] %>%
                           .[,.(id,fleet,fleetname,type,units)]

    }
    fleet_dt = rbindlist(fleet_dt.list)
    plot_dt = merge(plot_dt,fleet_dt,by=c("id","fleet")) %>%
              .[,id2:=paste0(fleetname,"_",units)]

    if(nrow(plot_dt) == 0){
      return()
    }

    # plot if models == 1
    if(length(selected_models)==1){
      # define fill
      p = plot_dt[,.(id,fleet,fleetname,units,time,sel_bio,sel_num)] 
      if(input_fill == "None"){
        p = p %>% .[,fill:=id]
      } else if(input_fill == "Fleet"){
         p = p %>% .[,fill:=fleetname] %>% 
            .[,fill:=factor(fill,levels=unique(p$fleetname))]
      } else {
         p = p %>% .[,fill:=c("mt","n (1000s)")[units]]
      }
      # define catch
      if(input_aggregate == "TRUE"){
        if(input_units == "Numbers"){
          p = p[,.(catch=sum(sel_num)),by=.(time,fill)]
          ylab_txt = "Catch (1000s N)"
        } else {
          p = p[,.(catch=sum(sel_bio)),by=.(time,fill)]
          ylab_txt = "Catch (mt)"
        }
        if(input_rescale == "TRUE"){
          p = p %>% 
            .[,catch:=catch/mean(catch)]
          ylab_txt = paste0(ylab_txt,": re-scaled")
        }
      } else {
        if(input_units == "Numbers"){
          p = p %>% 
            .[,facet:=factor(fleetname,levels=unique(p$fleetname))] %>%
            .[,.(catch=sum(sel_num)),by=.(time,facet,fill)]
            ylab_txt = "Catch (1000s N)"
        } else {
          p = p %>% 
            .[,facet:=factor(fleetname,levels=unique(p$fleetname))] %>%
            .[,.(catch=sum(sel_bio)),by=.(time,facet,fill)]
          ylab_txt = "Catch (mt)"
        }
        if(input_rescale == "TRUE"){
          p = p %>% 
            .[,catch:=catch/mean(catch),by=facet]
          ylab_txt = paste0(ylab_txt,": re-scaled")
        }
      }
      p = p %>%
        ggplot() +
			  ylim(0,NA) +
			  xlab("Time")
        if(input_aggregate=="FALSE"){
          n_row = 5
          p = p + facet_wrap(~facet,nrow=n_row)
        }
      p = p +
          ylab(ylab_txt)
        if(input_chart == "Bar"){
          p = p + geom_bar(aes(x=time,y=catch,fill=fill),position="stack",stat = "identity",color="white")
        } else {
          p = p + geom_area(aes(x=time,y=catch,fill=fill))
        }
      p = p + geom_hline(yintercept=0) + 
          viridis::scale_color_viridis("Category",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    viridis::scale_fill_viridis("Category",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
      return(p)
    } else {
        if(length(selected_models)>1&input_aggregate=="FALSE"&uniqueN(plot_dt$id2) != max(plot_dt$fleet)){
            return()
        }
        p = plot_dt[,.(id,fleet,fleetname,units,time,sel_bio,sel_num)] 
        p = p %>%
            .[,id:=factor(id,levels=unique(p$id))]
        if(input_units == "Numbers"){
          p = p[,.(catch=sum(sel_num)),by=.(id,time)]
          ylab_txt = "Catch (1000s N)"
        } else {
          p = p[,.(catch=sum(sel_bio)),by=.(id,time)]
          ylab_txt = "Catch (mt)"
        }
        if(input_rescale == "TRUE"){
          p = p %>% 
            .[,catch:=catch/mean(catch),by=.(id)]
          ylab_txt = paste0(ylab_txt,": re-scaled")
        }

         p = p %>% .[order(id,time)] %>%
        ggplot() +
			  ylim(0,NA) +
			  xlab("Time")
      p = p +
          ylab(ylab_txt)
          p = p + geom_path(aes(x=time,y=catch,color=id),linewidth=2)
        
      p = p + geom_hline(yintercept=0) + 
          viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
      return(p)
    }
  },
  height=function(){
    if(input$catch_by_fleet)
    {
      return((height_per_panel*1.5))
    } else {
      return((height_per_panel*1.5))
    }
  })


  # cpue_obs
  output$cpue_obs = renderPlot({
    input_models <- unique(filtered_id())
    
    if(length(input_models) < 1){
      return()
    }
    
    p <- create_cpue_comparison_plot(
      model_ids = input_models,
      model_stem = model_stem,
      show_se = input$cpue_se,
      show_fit = input$cpue_fit,
      use_log_scale = input$cpue_log,
      apply_varadj = input$cpue_varadj,
      filter_lambda = input$cpue_lambda
    )
    
    if(is.null(p)) {
      return()
    }
    
    return(p)
  },
  height = function() {
    return(height_per_panel * 1.5)
  })

  # Dynamic UI for fleet selection using shinyWidgets
  output$catch_fit_fleet_selector <- renderUI({
    input_models <- unique(filtered_id())
    
    if(length(input_models) < 1) {
      return(NULL)
    }
    
    # Get paths for selected models
    selected_models <- sapply(input_models, function(x) all_dirs[grep(x, all_dirs, fixed=TRUE)])
    
    # Check if catch.csv exists for models
    keep_selected <- which(sapply(file.path(model_stem, selected_models, "catch.csv"), file.exists) == TRUE)
    if(length(keep_selected) == 0) {
      return(NULL)
    } else {
      selected_models <- selected_models[keep_selected]
      input_models <- input_models[keep_selected]
    }
    
    # Read fleet information from all selected models
    all_fleets <- lapply(seq_along(selected_models), function(i) {
      model_path <- selected_models[i]
      catch_file <- file.path(model_stem, model_path, "catch.csv")
      
      if(file.exists(catch_file)) {
        catch_data <- fread(catch_file)
        
        # Try to get fleet names if available
        fleet_summary_file <- file.path(model_stem, model_path, "fleet_summary.csv")
        if(file.exists(fleet_summary_file)) {
          fleet_data <- fread(fleet_summary_file)
          return(data.table(Fleet = fleet_data$fleet, Fleet_name = fleet_data$fleetname))
        } else {
          return(data.table(Fleet = unique(catch_data$Fleet), 
                            Fleet_name = paste("Fleet", unique(catch_data$Fleet))))
        }
      } else {
        return(NULL)
      }
    })
    
    # Combine all fleet information
    all_fleets_dt <- rbindlist(all_fleets[!sapply(all_fleets, is.null)])
    
    # Get unique fleets with names
    unique_fleets <- unique(all_fleets_dt[, .(Fleet, Fleet_name)])
    
    # Store the fleet choices in a named vector
    fleet_choices <- setNames(unique_fleets$Fleet, unique_fleets$Fleet_name)
    
    # Create picker input with buttons for select/deselect all
    shinyWidgets::pickerInput(
      inputId = "catch_fit_fleets",
      label = "Select fleets to display",
      choices = fleet_choices,
      selected = unique_fleets$Fleet,
      multiple = TRUE,
      options = list(
        `actions-box` = TRUE,
        size = 10,
        `selected-text-format` = "count > 3"
      )
    )
  })

  # Render the catch fit plot
  output$catch_fit_plot <- renderPlot({
    input_models <- unique(filtered_id())
    
    if(length(input_models) < 1) {
      return()
    }
    
    # Get selected fleets (if any)
    selected_fleets <- input$catch_fit_fleets
    
    # Create the plot using our function
    p <- plot_catch_comparison(
      model_ids = input_models,
      model_stem = model_stem,
      fleets = selected_fleets,
      show_fit = input$catch_fit_show_lines,
      n_col = input$catch_fit_n_col,
      free_y_scale = input$catch_fit_free_y
    )
    
    return(p)
  }, height = function() {
    # Calculate appropriate height based on number of fleets and columns
    if(is.null(input$catch_fit_fleets)) {
      return(height_per_panel * 1.5)
    }
    
    n_fleets <- length(input$catch_fit_fleets)
    n_rows <- ceiling(n_fleets / input$catch_fit_n_col)
    
    return(max(height_per_panel * 1.5, height_per_panel * n_rows))
  })

  # Render the selectivity plot
  output$selex_plot <- renderPlot({
    input_models <- unique(filtered_id())
    
    if(length(input_models) < 1) {
      return()
    }
    
    # Create the plot using our function
    p <- plot_model_comparison_selex(
      model_ids = input_models,
      model_stem = model_stem,
      n_col = input$selex_n_col
    )
    
    return(p)
  }, height = function() {
    # Calculate appropriate height based on number of fleets and columns
    input_models <- unique(filtered_id())
    
    if(length(input_models) < 1) {
      return(height_per_panel * 1.5)
    }
    
    # Try to estimate the number of fleets
    selected_models <- sapply(input_models, function(x) all_dirs[grep(x, all_dirs, fixed=TRUE)])
    fleet_count <- 0
    
    for(model_path in selected_models) {
      selex_file <- file.path(model_stem, model_path, "selex_l.csv")
      if(file.exists(selex_file)) {
        selex_data <- fread(selex_file)
        fleet_count <- max(fleet_count, length(unique(selex_data$Fleet)))
      }
    }
    
    if(fleet_count == 0) {
      return(height_per_panel * 1.5)
    }
    
    n_rows <- ceiling(fleet_count / input$selex_n_col)
    return(max(height_per_panel * 1.5, height_per_panel * n_rows))
  })

  output$comp_plot <- renderPlot({
    input_models <- unique(filtered_id())
    
    if(length(input_models) < 1) {
      return()
    }
    
    # Convert input$comp_type to lowercase for function
    comp_type <- tolower(input$comp_type)
    
    # Create the plot using the provided function
    tryCatch({
      p <- plot_composition_comparison(
        model_ids = input_models,
        model_stem = model_stem,
        comp_type = comp_type,
        show_fit = input$comp_show_fit,
        n_col = input$comp_n_col,
        free_y_scale = input$comp_free_y
      )
      
      return(p)
    }, error = function(e) {
      # Create empty plot with error message
      p <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                label = paste("No composition data available:", e$message), 
                hjust = 0.5, vjust = 0.5) +
        theme_void()
      return(p)
    })
  }, height = function() {
    input_models <- unique(filtered_id())
    
    if(length(input_models) < 1) {
      return(height_per_panel * 1.5)
    }
    
    # Try to estimate the number of fleets with composition data
    selected_models <- sapply(input_models, function(x) all_dirs[grep(x, all_dirs, fixed=TRUE)])
    fleet_count <- 0
    
    comp_type <- tolower(input$comp_type)
    file_name <- if(comp_type == "length") "comp_len.csv" else "comp_size.csv"
    kind_filter <- if(comp_type == "length") "LEN" else "SIZE"
    
    for(model_path in selected_models) {
      comp_file <- file.path(model_stem, model_path, file_name)
      if(file.exists(comp_file)) {
        comp_data <- fread(comp_file)
        comp_data <- comp_data[Kind == kind_filter]
        fleet_count <- max(fleet_count, length(unique(comp_data$Fleet_name)))
      }
    }
    
    if(fleet_count == 0) {
      return(height_per_panel * 1.5)
    }
    
    n_rows <- ceiling(fleet_count / input$comp_n_col)
    return(max(height_per_panel * 1.5, height_per_panel * n_rows))
  })

} # End of server
