#' Compute Proximity Matrix
#'
#' @description
#' This function takes a numeric matrix and computes a square proximity matrix (similarity or distance) based on a specified method.
#'
#' @param data A numeric matrix with n rows and p columns. Each row typically represents an observation.
#' @param proxType An integer specifying the type of proximity measure to use.
#' @param side An integer indicating the direction for computing proximity.
#' @param isContainMissingValue An integer indicating whether the input data contains missing values.
#'
#' @return A square matrix representing the proximity between rows or columns, depending on the selected side.
#'
#' @details
#'
#' **proxType**
#'
#' Available proxType options include:
#'
#' \itemize{
#'   \item \code{0}: Euclidean
#'   \item \code{1}: Pearson correlation
#'   \item \code{2}: Kendall correlation
#'   \item \code{3}: Spearman correlation
#'   \item \code{4}: Adjusted tangent correlation (atancorr)
#'   \item \code{5}: City-block (Manhattan) distance
#'   \item \code{6}: Absolute Pearson correlation
#'   \item \code{7}: Uncentered correlation
#'   \item \code{8}: Absolute uncentered correlation
#'   \item \code{20}: Hamman similarity (binary)
#'   \item \code{21}: Jaccard index (binary)
#'   \item \code{22}: Phi coefficient (binary)
#'   \item \code{23}: Rao coefficient (binary)
#'   \item \code{24}: Rogers-Tanimoto similarity (binary)
#'   \item \code{25}: Simple matching coefficient (binary)
#'   \item \code{26}: Sneath coefficient (binary)
#'   \item \code{27}: Yule's Q (binary)
#' }
#'
#' Ensure the data type matches the selected method. For example, binary methods should only be used on binary (0/1) data.
#'
#' **side**
#'
#' Use \code{0} for row-wise proximity and \code{1} for column-wise proximity.
#'
#' **isContainMissingValue**
#'
#' Set to \code{1} if the input data includes missing values; otherwise, use \code{0}.
#'
#' @export
computeProximity <- function(data, proxType, side, isContainMissingValue) {
  .Call('_GAPR_computeProximity_R', as.matrix(data), as.integer(proxType),
        as.integer(side), as.integer(isContainMissingValue))
}

#' Ellipse Sort
#'
#' @description
#' This function applies an Rank-2-Ellipse seriation method to reorder a proximity matrix.
#'
#' @param data A square numeric proximity matrix (either n × n or p × p), representing pairwise distances or similarities between items.
#'
#' @return An integer vector representing the reordered indices of the matrix rows.
#'
#' @export
ellipse_sort <- function(data) {
  .Call('_GAPR_ellipse_sort_R', as.matrix(data))
}

#' HCTree Sort
#'
#' @description
#' This function applies a hierarchical clustering tree (HCT) sorting algorithm
#' to reorder rows of a proximity matrix. It supports external ordering constraints,
#' different linkage-based order types, and optional flipping for optimal layout.
#'
#' @param distance_matrix A square numeric proximity matrix (either n × n or p × p) representing pairwise distances between items.
#' @param externalOrder An integer vector specifying an initial or external ordering of the items (can be empty or NULL if not used).
#' @param orderType An integer indicating the type of hierarchical clustering order to apply.
#' @param flipType An integer indicating the flipping methods.
#'
#' @return A list representing a dendrogram tree structure, containing:
#' \code{left}, \code{right}, and \code{height} for tree construction,
#' and \code{order} for the optimal leaf order.
#'
#' @details
#'
#' **distance_matrix**
#'
#' The input matrix must represent pairwise distances between items.
#' If you start with a similarity matrix (e.g., a correlation matrix), you must convert it to a dissimilarity matrix before use.
#' For example, for correlation-based similarities,
#' use \code{as.matrix(as.dist(1 - cor_matrix))} or other appropriate transformations to convert it to a proper distance matrix.
#' The matrix should also be symmetric and non-negative.
#'
#' **orderType**
#'
#' Specifies the linkage method used for hierarchical clustering:
#' \itemize{
#'   \item \code{0}: Single-linkage
#'   \item \code{1}: Complete-linkage
#'   \item \code{2}: Average-linkage (UPGMA)
#' }
#'
#' **flipType**
#'
#' Controls how the branches of the clustering tree are flipped:
#' \itemize{
#'   \item \code{1}: Flip based on \code{externalOrder} This option should be used only when \code{externalOrder} is provided.
#'   \item \code{2}: Uncle-flipping
#'   \item \code{3}: Grandpa-flipping
#' }
#'
#' \strong{Important:}
#' Do not specify both \code{flipType = 1} and a \code{NULL} or missing \code{externalOrder}.
#' When using \code{flipType = 1}, \code{externalOrder} must be a valid integer vector.
#'
#' @export
hctree_sort <- function(distance_matrix, externalOrder = NULL, orderType, flipType) {
  .Call('_GAPR_hctree_sort_R', as.matrix(distance_matrix), as.integer(externalOrder),
        as.integer(orderType), as.integer(flipType))
}
