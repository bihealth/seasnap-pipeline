############################################
## define functions for generic calculations
############################################

# get the number suffix for entry ids based on the contrast list (in config file)
entry_num_id <- function(entry_name, config) {
	titles <- sapply(config$contrasts$contrast_list, function(x) x$title)
	return(paste0("_ID", which(titles==entry_name)[1] - 1))
}

