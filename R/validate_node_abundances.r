
#' Validate tip abundances and remove zero-abundance tips
#'
#' @param tree A phylo object.
#' @param node_abundances A two-column data frame containing
#'   tip labels and their non-negative abundances.
#'
#' @return A list containing the validated/pruned tree
#'   and its corresponding node_abundances data frame.
#'
#' @keywords internal
validate_node_abundances <- function(tree, node_abundances) {

  # 1. Validate tree
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be a phylo object.", call. = FALSE)
  }

  tips <- tree$tip.label

  if (anyNA(tips) ||
      any(!nzchar(tips)) ||
      anyDuplicated(tips)) {
    stop(
      "Tree tip labels must be non-missing, non-empty and unique.",
      call. = FALSE
    )
  }

  # 2. Validate input type
  if (!is.data.frame(node_abundances)) {
    stop(
      "`node_abundances` must be a data frame with exactly ",
      "two columns: tip labels and abundances.",
      call. = FALSE
    )
  }

  # 3. Validate dimensions
  if (ncol(node_abundances) != 2L) {
    stop(
      "`node_abundances` must have exactly two columns.",
      call. = FALSE
    )
  }

  if (nrow(node_abundances) != length(tips)) {
    stop(
      "`node_abundances` must contain exactly one row per ",
      "tree tip (expected ", length(tips), ", received ",
      nrow(node_abundances), ").",
      call. = FALSE
    )
  }

  labels <- node_abundances[[1L]]
  values <- node_abundances[[2L]]

  # 4. Validate names
  if (!is.character(labels) ||
      anyNA(labels) ||
      any(!nzchar(labels))) {
    stop(
      "The first column of `node_abundances` must contain ",
      "non-missing, non-empty character tip labels.",
      call. = FALSE
    )
  }

  if (anyDuplicated(labels)) {
    stop(
      "Duplicate tip labels found in `node_abundances`: ",
      paste(
        unique(labels[duplicated(labels)]),
        collapse = ", "
      ),
      call. = FALSE
    )
  }

  # 5. Check exact correspondence
  missing_tips <- setdiff(tips, labels)
  unknown_tips <- setdiff(labels, tips)

  if (length(missing_tips) > 0L ||
      length(unknown_tips) > 0L) {

    message_parts <- c(
      "Tip labels in `node_abundances` do not match tree$tip.label."
    )

    if (length(missing_tips) > 0L) {
      message_parts <- c(
        message_parts,
        paste0(
          "Missing tree tips: ",
          paste(missing_tips, collapse = ", ")
        )
      )
    }

    if (length(unknown_tips) > 0L) {
      message_parts <- c(
        message_parts,
        paste0(
          "Unknown labels: ",
          paste(unknown_tips, collapse = ", ")
        )
      )
    }

    stop(
      paste(message_parts, collapse = "\n"),
      call. = FALSE
    )
  }

  # 6. Validate numeric abundances
  if (!is.numeric(values) ||
      any(!is.finite(values))) {
    stop(
      "The second column of `node_abundances` must ",
      "contain finite numeric values (no NA, NaN or Inf).",
      call. = FALSE
    )
  }

  if (any(values < 0)) {
    stop(
      "Abundance values cannot be negative.",
      call. = FALSE
    )
  }

  # 7. Match the order of tree tip labels
  node_abundances <- node_abundances[
    match(tips, labels),
    ,
    drop = FALSE
  ]

  rownames(node_abundances) <- NULL
  values <- node_abundances[[2L]]

  # 8. Identify zero-abundance tips
  zero_tips <- tips[values == 0]

  if (length(zero_tips) > 0L) {

    remaining <- sum(values > 0)

    if (remaining < 2L) {
      stop(
        "Removing zero-abundance tips would leave fewer ",
        "than two tips. Indices cannot be calculated.",
        call. = FALSE
      )
    }

    # Prune zero-abundance tips
    tree <- ape::drop.tip(tree, zero_tips)

    # Keep only corresponding abundance rows
    node_abundances <- node_abundances[
      match(tree$tip.label, node_abundances[[1L]]),
      ,
      drop = FALSE
    ]

    rownames(node_abundances) <- NULL
  }

  # 9. Return aligned objects
  list(
    tree = tree,
    node_abundances = node_abundances
  )
}
