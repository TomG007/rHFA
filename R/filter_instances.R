#' @title Count Instances
#'
#' @description Count the number of times each category of `group_col` appears
#' per `rep_col`.
#'
#' @param data A *data.frame* or *data.table*.
#' @param group_col *character* or *numeric*: Column in `data` for primary grouping.
#' @param rep_col *character* or *numeric*: Column in `data` indicating replicates of
#' `group_col`.
#' @param na.rm *logical* Remove rows containing missing values? Default is `FALSE`.
#'
#' @return A *named vector*. Names correspond to levels in `group_col` and values represent unique counts.
#' @seealso `filter_instances()`, `autofilter_instances()`
#'
#' @import data.table
#' @export
count_instances <- function(data, group_col, rep_col, na.rm = FALSE) {
  # Remove missing values if na.rm is TRUE
  if (na.rm) data <- na.omit(data)

  # Convert to data.table if not already
  setDT(data)

  # Select unique combinations of group_col and rep_col
  unique_data <- unique(data[, .(get(group_col), get(rep_col))])

  # Count instances of rep_col within each group_col
  instance_counts <- unique_data[, .N, by = get(group_col)]$N
  names(instance_counts) <- unique_data[, get(group_col)]

  return(instance_counts)
}

#' @title Filter Instances
#'
#' @description Keep categories of `group_col` that occur in at least `min_times` of `rep_col` categories.
#'
#' @param data A *data.frame* or *data.table*.
#' @param group_col *character* or *numeric*: Column in `data` for primary grouping.
#' @param rep_col *character* or *numeric*: Column in `data` indicating replicates of `group_col`.
#' @param min_times *numeric*: Minimum number of instances (e.g., minimum site-years). Default is 2.
#' @param na.rm *logical*: Remove rows containing missing values? Default is `FALSE`.
#'
#' @return A *data.frame* or *data.table*, a filtered version of `data`.
#' @seealso `count_instances()`, `autofilter_instances()`
#'
#' @import data.table
#' @export
filter_instances <- function(data, group_col, rep_col, min_times = 2, na.rm = FALSE) {
  # Count unique instances of grouping in rep_col
  instance_counts <- count_instances(data, group_col, rep_col, na.rm)

  # Identify those with at least min_times instances
  keep_groups <- names(instance_counts[instance_counts >= min_times])

  # Filter the data to keep only rows from those groups
  setDT(data)
  filtered_data <- data[get(group_col) %in% keep_groups]

  return(droplevels(filtered_data))
}

#' @title Autofilter Instances
#'
#' @description Filters geno-years, site-years, and genos per site-year to `min_times` instances.
#'
#' @param data A *data.frame* or *data.table*.
#' @param site *character* or *numeric*: Column in `data` for site grouping.
#' @param year *character* or *numeric*: Column in `data` for year grouping.
#' @param geno *character* or *numeric*: Column in `data` for genotype grouping.
#' @param min_times *numeric*: Minimum number of instances (e.g., minimum site-years). Default is 2.
#' @param na.rm *logical*: Remove rows containing missing values? Default is `TRUE`.
#' @param max_cycles *numeric*: Maximum number of times to iteratively filter. Default is 999.
#'
#' @details Calls `filter_instances()` iteratively for the combinations of {geno, year}, {site, year}, and {siteyear, geno}.
#' Repeats until no more rows are removed or `max_cycles` is reached.
#' `siteyear` is created using `paste(data[, site], data[, year], sep = "___")`.
#'
#' @return A *data.frame* or *data.table*, a filtered version of `data`.
#' @seealso `filter_instances()`
#'
#' @import data.table
#' @export
autofilter_instances <- function(data, site, year, geno, min_times = 2, na.rm = TRUE, max_cycles = 999) {
  # Set rownames if not present
  if (is.null(rownames(data))) {
    rownames(data) <- seq_len(nrow(data))
  }

  # Create siteyear column for site-year grouping
  data[, siteyear := paste(get(site), get(year), sep = "___")]

  # Initialize cycle control
  nr <- nrow(data)
  new_nr <- nr - 1
  cycles <- 0

  # Repeat filtering until no more rows are removed or max_cycles is reached
  while (nr != new_nr && cycles < max_cycles) {
    nr <- nrow(data)

    # Filter geno-year
    data <- filter_instances(data, geno, year, min_times, na.rm)

    # Filter site-year
    data <- filter_instances(data, site, year, min_times, na.rm)

    # Filter genos within site-year
    data <- filter_instances(data, "siteyear", geno, min_times, na.rm)

    new_nr <- nrow(data)
    cycles <- cycles + 1
  }

  # If cycles were reached, give a message
  if (cycles < max_cycles) {
    message(sprintf("Filtered geno-years, site-years, and genos per site-year to at least %d instances over %d cycle%s.",
                    min_times, cycles, ifelse(cycles > 1, "s", "")))
  } else {
    message(sprintf("Stopped filtering after %d cycles. Consider increasing `max_cycles`.", cycles))
  }

  return(data)
}
