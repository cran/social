#' All paths between two nodes
#'
#' Estimate all the possible paths between two nodes in a 
#'     simple graph using the stochastic method described 
#'     by Roberts, B. & Kroese, D.P. (2007) Estimating the 
#'     number of s-t paths in a graph. Journal of Graph 
#'     Algorithms and Applications 11(1), 195-214.
#'
#' @param A a (possibly weighted) adjacency matrix.
#' @param start.node the index of the vertex from which the paths will be calculated.
#' @param end.node the index of the vertex to which the paths will be calculated.
#' @param max.depth the maximum length of the paths to the returned.
#' @param n.pilot the number of naive paths to generate (see Roberts & Kroese, 2007).
#' @param n.estimate the number of paths to generate (see Roberts & Kroese, 2007).
#'
#' @return An estimate of all the unique paths between \code{start.node} and \code{end.node} as an nrow(A)xN matrix padded with zeros.
#'
#' @examples
#' # Using the data from Figure 1 in Roberts & Kroese (2007)
#' A = matrix(c(0,1,0,1,0,
#'              1,0,0,1,1,
#'              0,0,0,1,1,
#'              1,1,1,0,0,
#'              0,1,1,0,0), nrow=5)
#' paths = social.all.paths(A, 1, 5)
#'
#' @export
social.all.paths = function(A, start.node, end.node, max.depth=nrow(A), n.pilot=5000, n.estimate=10000) {
	paths = allpaths(A, start.node, end.node, max.depth, n.pilot, n.estimate)
	if(is.null(paths)) {
		paths = matrix(ncol=ncol(A), nrow=0)
	}
	return(paths)
}
