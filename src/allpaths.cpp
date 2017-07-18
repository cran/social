#include <Rcpp.h>
using namespace Rcpp;

//// [[Rcpp::export]]
int weightedsample(NumericVector nodes, NumericVector weights) {
	// Normalise weights
	int n = weights.length();
	int sumOfWeights = 0;
	for(int i=0; i<n; i++) {
		sumOfWeights += weights(i);
	}
	for(int i=0; i<n; i++) {
		weights(i) = weights(i)/sumOfWeights;
	}
	// Choose random weighted sample
	NumericVector randomNumber = runif(1);
	for(int i=0; i<n; i++) {
		if(randomNumber(0) < weights(i)) {
			return nodes(i);
		}
		else {
			randomNumber(0) -= weights(i);
		}
	}
	return nodes(n-1); // Shouldn't need this
}

//// [[Rcpp::export]]
NumericVector findnodes(NumericMatrix adjMatrix, int currentNode) {
	int numNodes = adjMatrix.nrow();
	NumericVector tmp(numNodes);
	int numNeighbours = 0;
	for(int j=0; j<numNodes; j++) {
		if(adjMatrix(currentNode, j) > 0) {
			tmp(numNeighbours) = j;
			numNeighbours++;
		}
	}
	NumericVector nextNodes(numNeighbours);
	for(int j=0; j<numNeighbours; j++) {
		nextNodes(j) = tmp(j);
	}
	return nextNodes;
}

//// [[Rcpp::export]]
bool isuniquepath(NumericMatrix paths, int i, LogicalVector keep) {
	int numNodes = paths.ncol();
	for(int j=0; j<i; j++) {
		if(keep(j)) {
			int sum = 0;
			for(int k=0; k<numNodes; k++) {
				if(paths(i, k) == paths(j, k)) {
					sum++;
				}
			}
			if(sum == numNodes) {
				return false;
			}
		}
	}
	return true;
}

//// [[Rcpp::export]]
bool containsendnode(NumericMatrix paths, int i, int endNode) {
	int numNodes = paths.ncol();
	for(int j=0; j<numNodes; j++) {
		if(paths(i, j) == endNode) {
			return true;
		}
	}
	return false;
}

//// [[Rcpp::export]]
NumericMatrix clean(NumericMatrix paths, int endNode) {
	int numPaths = paths.nrow();
	int numNodes = paths.ncol();
	// Keep path?
	LogicalVector keep(numPaths, false);
	int numUniquePaths = 0;
	for(int i=0; i<numPaths; i++) {
		if(containsendnode(paths, i, endNode) && isuniquepath(paths, i, keep)) {
			keep(i) = true;
			numUniquePaths++;
		}
	}
	// Remove unwanted paths
	NumericMatrix goodPaths(numUniquePaths, numNodes);
	int n = 0;
	for(int i=0; i<numPaths; i++) {
		if(keep(i)) {
			for(int k=0; k<numNodes; k++) {
				goodPaths(n, k) = paths(i, k);
			}
			n++;
		}
	}
	return goodPaths;
}

//// [[Rcpp::export]]
NumericMatrix naivepaths(NumericMatrix adjMatrix, int startNode, int endNode, int nPilot) {
	int numNodes = adjMatrix.nrow();
	NumericMatrix naivePaths(nPilot, numNodes);
	for(int i=0; i<nPilot; i++) {
		for(int j=0; j<numNodes; j++) {
			naivePaths(i, j) = -1;
		}
	}
	for(int i=0; i<nPilot; i++) {
		NumericMatrix adjMatrixCopy = clone(adjMatrix);
		naivePaths(i, 0) = startNode;
		int pathLength = 1;
		int currentNode = startNode;
		while(currentNode != endNode) {
			NumericVector nextNodes = findnodes(adjMatrixCopy, currentNode);
			if(nextNodes.length() == 0) {
				break;
			}
			else {
				for(int j=0; j<numNodes; j++) {
					adjMatrixCopy(currentNode, j) = 0;
					adjMatrixCopy(j, currentNode) = 0;
				}
				NumericVector weights(nextNodes.length(), 1.0);
				currentNode = weightedsample(nextNodes, weights);
				naivePaths(i, pathLength) = currentNode;
				pathLength++;
			}
		}
	}
	NumericMatrix goodPaths = clean(naivePaths, endNode);
	return goodPaths;
}

//// [[Rcpp::export]]
NumericVector lengthdistribution(NumericMatrix adjMatrix, NumericMatrix goodPaths) {
	int numNodes = adjMatrix.nrow();
	int numUniquePaths = goodPaths.nrow();
	// Get the lengths of all paths
	NumericVector pathLengths(numUniquePaths);
	for(int i=0; i<numUniquePaths; i++) {
		int sum = 0;
		for(int j=0; j<numNodes; j++) {
			if(goodPaths(i, j) > -1) {
				sum++;
			}
		}
		pathLengths(i) = sum;
	}
	// Calculate the length-distribution vector, p
	NumericVector p(numNodes);
	for(int k=0; k<numNodes; k++) {
		double n1 = 0.0;
		double n2 = 0.0;
		for(int i=0; i<numUniquePaths; i++) {
			if(pathLengths(i) == k) {
				n1 = n1 + 1.0;
			}
			if(pathLengths(i) >= k) {
				n2 = n2 + 1.0;
			}
		}
		p(k) = n1/n2;
		if(p(k) == 0) {
			p(k) = 1.0/(double)numUniquePaths;
		}
		if(p(k) == 1) {
			p(k) = 1.0 - 1.0/(double)numUniquePaths;
		}
	}
	return p;
}

//// [[Rcpp::export]]
NumericMatrix estimatedpaths(NumericMatrix adjMatrix, int startNode, int endNode, int maxDepth, int nEstimation, NumericVector p) {
	int numNodes = adjMatrix.nrow();
	NumericMatrix paths(nEstimation, numNodes);
	for(int i=0; i<nEstimation; i++) {
		for(int j=0; j<numNodes; j++) {
			paths(i, j) = -1;
		}
	}
	for(int i=0; i<nEstimation; i++) {
		NumericMatrix adjMatrixCopy = clone(adjMatrix);
		paths(i, 0) = startNode;
		int pathLength = 1;
		int currentNode = startNode;
		while(currentNode != endNode) {
			if(adjMatrixCopy(currentNode, endNode) == 0) {
				NumericVector nextNodes = findnodes(adjMatrixCopy, currentNode);
				if(nextNodes.length() == 0) {
					break;
				}
				else {
					for(int j=0; j<numNodes; j++) {
						adjMatrixCopy(currentNode, j) = 0;
						adjMatrixCopy(j, currentNode) = 0;
					}
					NumericVector weights(nextNodes.length(), 1.0);
					currentNode = weightedsample(nextNodes, weights);
					paths(i, pathLength) = currentNode;
					pathLength++;
					if(pathLength > maxDepth) {
						break;
					}
				}
			}
			else {
				NumericVector nextNodes = findnodes(adjMatrixCopy, currentNode);
				if(nextNodes.length() == 1) {
					paths(i, pathLength) = endNode;
					break;
				}
				else {
					NumericVector randomNumber = runif(1);
					if(randomNumber(0) < p(pathLength)) {
						paths(i, pathLength) = endNode;
						break;
					}
					else {
						adjMatrixCopy(currentNode, endNode) = 0;
						adjMatrixCopy(endNode, currentNode) = 0;
					}
				}
			}
		}
	}
	NumericMatrix allPaths = clean(paths, endNode);
	return allPaths;
}

// [[Rcpp::export]]
NumericMatrix allpaths(NumericMatrix adjMatrix, int startNode, int endNode, int maxDepth, int nPilot, int nEstimation) {

	// One to zero indexing (R to Cpp)
	startNode--;
	endNode--;
	
	// Naive path generation
	NumericMatrix goodPaths = naivepaths(adjMatrix, startNode, endNode, nPilot);	
	
	// Calculate the length-distribution vector
	NumericVector p = lengthdistribution(adjMatrix, goodPaths);
	
	// Length-distribution method
	NumericMatrix allPaths = estimatedpaths(adjMatrix, startNode, endNode, maxDepth, nEstimation, p);
	
 	// Zero to one indexing (Cpp to R)
 	for(int i=0; i<allPaths.nrow(); i++) {
 		for(int j=0; j<adjMatrix.nrow(); j++) {
 			allPaths(i, j) += 1;
 		}
 	}
    return allPaths;
}
