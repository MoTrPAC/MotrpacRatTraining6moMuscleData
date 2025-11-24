#' @title Overrepresentation analysis for specific gene sets
#'
#' @param input the set of features you would like to investigate
#' @param background the set of all features included/studied
#' @param database the gene sets you would like to include. See `rat_gene_sets` for information about the specific format that is expected.
#' @param min_size Minimum number of features in a gene set included in the `background` object
#'
#' @return a data.frame object that details enrichment results for each of the sets describes in \code{rat_gene_sets}
#' or whatever is used for the `database` input.
#' @importFrom dplyr %>% mutate across everything filter arrange
#' @importFrom data.table :=
#' @importFrom stats phyper
#'

muscle_ora <- function(input,
                       background,
                       database = rat_gene_sets,
                       min_size = 5L) {
  on.exit(gc())

  if (!is.vector(input, mode = "character") || !length(input)) {
    stop("`input` must be a character vector of features.")
  }

  if (!is.vector(background, mode = "character") || !length(background)) {
    stop("`background` must be a character vector of features.")
  }

  input <- unique(input[!is.na(input)])
  background <- unique(background[!is.na(background)])

  if (any(!input %in% background)) {
    stop("`input` must be a subset of `background`.")
  }

  ## Prepare molecular signatures ----
  index <- database

  set_size_DB <- lengths(index)

  if (all(grepl("^PTMSIGDB", names(index)))) {
    # This code is from TMSig::filterSets. It was repurposed to work with
    # PTMsigDB, where the sites in each set end with ";u" or ";d", but the
    # background vector of IDs do not.
    set_dt <- data.table(
      sets = rep(names(index), lengths(index)),
      elements = unlist(index),
      stringsAsFactors = FALSE
    )

    set_dt[, elements2 := sub(";.*$", "", elements)]
    set_dt <- subset(set_dt, subset = elements2 %in%  background)

    index <- split(x = set_dt[["elements2"]], f = set_dt[["sets"]])
    set_sizes <- lengths(index)
    keep_sizes <- (set_sizes >= min_size) & (set_sizes < length(background))
    index <- index[keep_sizes]
  } else {
    # Restrict sets to background, filter by size
    index <- TMSig::filterSets(
      x = index,
      background = background,
      min_size = min_size,
      max_size = length(background) - 1L
    )
  }

  # overlap_cutoff does not affect RefMet, PhosphoSitePlus, or PTMSigDB
  # # signatures
  # if (!any(grepl("^REFMET|^PSP|^PTMSIGDB", names(index)))) {
  #   keep <- lengths(index) / set_size_DB[names(index)]
  #
  #   index <- index[keep]
  # }

  # Restrict sets to input vector of features.
  index_filt <- .fast_list_intersect(x = index, y = input)
  # New chnge as of Oct 8th 2025: empty sets are NOT kept.
  #so we dont test any pathways without a matching value.
  index_filt = index_filt[lengths(index_filt) > 0]

  out <- data.frame(set = names(index),
                    set_size = lengths(index)) %>%
    mutate(set_size_DB = set_size_DB[set],
           size_ratio = round(set_size / set_size_DB, digits = 3L),
           set_size_in_input = lengths(index_filt)[set],
           input_size = length(input),
           background_size = length(background),
           p_value = phyper(q = set_size_in_input - 1L,
                            m = set_size,
                            n = background_size - set_size,
                            k = input_size,
                            lower.tail = FALSE),
           across(.cols = everything(),
                  .fns = ~ structure(.x, names = NULL))) %>%
    arrange(p_value) %>%
    filter(!is.na(p_value))

  return(out)
}
# Note: does not check that x is a valid named list of character vectors or that
# y is a character vector.
.fast_list_intersect <- function(x, y) {
  # Convert list to data.table for fast filtering
  dt <- data.table(sets = rep(names(x), lengths(x)),
                   elements = unlist(x),
                   stringsAsFactors = FALSE)

  # Convert to factor to preserve order and keep empty sets when splitting
  dt[, sets := factor(sets, levels = unique(names(x)))]

  setorderv(dt, cols = c("sets", "elements"), order = rep(1L, 2L))

  y <- unique(y)
  y <- y[!is.na(y)]
  dt <- subset(dt, subset = elements %in% y)
  dt <- unique(dt)

  out <- split(x = dt[["elements"]], f = dt[["sets"]])

  return(out)
}

