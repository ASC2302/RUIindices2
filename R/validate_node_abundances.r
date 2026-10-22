
#' Validate tip abundances and prune zero-abundance tips
#'
#' @param tree An object of class "phylo".
#' @param node_abundances A two-column data frame containing
#'   tip labels and non-negative numeric abundances.
#'
#' @return A list containing the validated/pruned tree,
#'   corresponding abundances and removed tip labels.
#'
#' @keywords internal
validate_node_abundances <- function(tree, node_abundances) {

  if (!inherits(tree, "phylo")) {
    stop("`tree` must be an object of class 'phylo'.",
         call. = FALSE)
  }

  if (anyNA(tree$tip.label) ||
      anyDuplicated(tree$tip.label)) {
    stop("Tree tip labels must be unique and non-missing.",
         call. = FALSE)
  }

  if (!is.data.frame(node_abundances)) {
    stop("`node_abundances` must be a data frame.",
         call. = FALSE)
  }

  if (ncol(node_abundances) != 2L) {
    stop(
      "`node_abundances` must have exactly two columns: ",
      "tip label and abundance.",
      call. = FALSE
    )
  }

  if (nrow(node_abundances) != ape::Ntip(tree)) {
    stop(
      "`node_abundances` must have exactly ",
      ape::Ntip(tree),
      " rows, one for each tree tip.",
      call. = FALSE
    )
  }

  labels <- node_abundances[[1]]
  values <- node_abundances[[2]]

  if (!is.character(labels) ||
      anyNA(labels) ||
      any(!nzchar(labels))) {
    stop(
      "The first column must contain non-missing, ",
      "non-empty character tip labels.",
      call. = FALSE
    )
  }

  if (anyDuplicated(labels)) {
    stop(
      "Duplicate tip labels found in `node_abundances`: ",
      paste(unique(labels[duplicated(labels)]),
            collapse = ", "),
      call. = FALSE
    )
  }

  missing_tips <- setdiff(tree$tip.label, labels)
  extra_tips <- setdiff(labels, tree$tip.label)

  if (length(missing_tips) || length(extra_tips)) {
    stop(
      "Tip labels in `node_abundances` do not match ",
      "the tree.\n",
      "Missing labels: ",
      paste(missing_tips, collapse = ", "),
      "\nUnknown labels: ",
      paste(extra_tips, collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.numeric(values) ||
      is.complex(values)) {
    stop(
      "The second column must contain numeric abundances.",
      call. = FALSE
    )
  }

  if (any(!is.finite(values))) {
    stop(
      "Abundances must not contain NA, NaN or Inf.",
      call. = FALSE
    )
  }

  if (any(values < 0)) {
    stop(
      "Abundances must be non-negative.",
      call. = FALSE
    )
  }

  # Match abundances by label, not by row position.
  ordered_values <- values[
    match(tree$tip.label, labels)
  ]

  # Identify and remove zero-abundance tips.
  zero_tips <- tree$tip.label[ordered_values == 0]

  if (length(zero_tips) > 0L) {

    remaining <- sum(ordered_values > 0)

    if (remaining < 2L) {
      stop(
        "Removing zero-abundance tips would leave ",
        "fewer than two tips. Cannot calculate indices.",
        call. = FALSE
      )
    }

    tree <- ape::drop.tip(tree, zero_tips)
  }

  # Reconstruct the abundance table in the new tip order.
  index <- match(tree$tip.label, labels)

  clean_abundances <- data.frame(
    tip_label = tree$tip.label,
    abundance = values[index],
    stringsAsFactors = FALSE
  )

  return(list(
    tree = tree,
    node_abundances = clean_abundances,
    removed_tips = zero_tips
  ))
}