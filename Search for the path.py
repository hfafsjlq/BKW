### A valid rounded reduction path with optimal max-complexity

import math
from typing import List, Tuple, Optional

def search_optimal_path(n: int, sigma: float, S: List[int], eta: float, q: int, b: int, epsilon: float, d: int) -> Tuple[Optional[List[Tuple]], bool]:
   
    # Step 1: Construct the graph G using the provided parameters
    V, E = Graph_const(n, S, sigma, eta, b, q, d)

    # Initialize data structures
    Sigma = {}  # Stores noise standard deviation for each vertex
    Edge = {}   # Stores the edge reaching each vertex

    # Initialize for all η1 ∈ S
    for eta1 in S:
        Sigma[(n, eta1, 0)] = sigma**2
        Edge[(n, eta1, 0)] = None
        for st in range(1, 5):
            Sigma[(n, eta1, st)] = -float('inf')
            Edge[(n, eta1, st)] = None

    # Iterate from j = n down to 1
    for j in range(n, 0, -1):
        for eta2 in sorted(S, reverse=True):
            # Initialize for current (j, eta2)
            for st in range(1, 5):
                Sigma[(j, eta2, st)] = 0
                Edge[(j, eta2, st)] = None

            # Process each edge to (j, eta2, st')
            for st_prime in range(1, 5):
                for e in E:
                    if e[1] == (j, eta2, st_prime):  # Edge e reaches (j, eta2, st')
                        i, eta1, st = e[0]  # Origin of edge e
                        c1, c2, _ = e[2]    # Coefficients from edge e

                        # Update Sigma and Edge if a better path is found
                        new_sigma = c1 * Sigma[(i, eta1, st)] + c2
                        if new_sigma <= Sigma[(j, eta2, st_prime)]:
                            Sigma[(j, eta2, st_prime)] = new_sigma
                            Edge[(j, eta2, st_prime)] = e

            # Check if the current path satisfies the conditions
            if eta2 > 4 * math.log((2 * d + 1)**j / epsilon) * (1 - (2 * math.pi**2 * Sigma[(j, eta2, st_prime)]) / q**2) ** (-2 ** (n / b)) and j + math.log2(j) <= eta:
                # Reconstruct the path
                path = []
                current_edge = Edge[(j, eta2, st_prime)]
                while current_edge is not None:
                    path.append(current_edge)
                    current_edge = Edge[current_edge[0]]  # Move to the previous edge
                path.reverse()
                return path, True

    # If no valid path is found
    return None, False



def find_valid_rounded_path(n: int, sigma: float, eta: float, q: int, b: int, epsilon: float, d: int, prec: float = 0.1) -> Optional[List[Tuple]]:
    
    found = None
    rise = n  

    while rise > prec:
        rise /= 2  # Halve the step size
        S = [x * prec for x in range(int((eta - rise) / prec) + 1)]  
        S = [x for x in S if x <= eta - rise]  

        # Call (search_optimal_path) 
        output, success = search_optimal_path(n, sigma, S, eta - rise, q, b, epsilon, d)

        if success:
            found = output  
            eta -= rise  

    return found
