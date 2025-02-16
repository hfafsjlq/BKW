# Construction of the graph for the LWE-solving BKW algorithm

import math

# Define the approximation function for S
def Approx_S(x, S):
    
    # Ensure S is sorted and contains no duplicates
    S = sorted(set(S))
    # Compute ceiling and floor approximation
    ceil_approx = min([s for s in S if s >= x], default=None)
    floor_approx = max([s for s in S if s <= x], default=None)
    if ceil_approx is None and floor_approx is None:
        raise ValueError(f"No valid approximation found for {x} in {S}")
    return ceil_approx if ceil_approx is not None else floor_approx


# Construction of the graph G
def construct_graph(n, S, sigma, eta, b, q, d):
    
    V = [(i, eta1, st) for i in range(1, n+1) for eta1 in S for st in range(1, 5)]
    E = []

    # LF1-reduce
    for i in range(1, n+1):
        for eta1 in S:
            if i - b > 0:
                j = i - b
                try:
                    term = math.log2(2**eta1 - (q**b - 1) / 2)
                    eta2 = Approx_S(term, S)
                except ValueError:
                    continue # Skip invalid eta2
                rop = math.log2(i) + eta1
                if rop > eta:
                    continue  # Skip if complexity exceeds eta
                for st1 in range(1, 5):
                    for st2 in range(1, 5):
                        E.append(((i, eta1, st1), (j, eta2, st2), (2, 0, 'LF1-reduce')))

    # LF2-reduce
    for i in range(1, n+1):
        for eta1 in S:
            if i - b > 0:
                j = i - b
                try:
                    term = math.log2((2**(2*eta1) / (q**b - 1) -2**eta1/ 2))
                    eta2 = Approx_S(term, S)
                except ValueError:
                    continue  # Skip invalid eta2
                rop = math.log2(i) + max(eta1, eta2)
                if rop > eta:
                    continue  # Skip if complexity exceeds eta
                for st1 in range(1, 5):
                    for st2 in range(1, 5):
                        E.append(((i, eta1, st1), (j, eta2, st2), (2, 0, 'LF2-reduce')))


    # Discard-reduce
    for i in range(1, n+1):
        for eta1 in S:
            if i - b > 0:
                j = i - b
                try:
                    eta2 = Approx_S(eta1 - math.log2(q**b), S)
                except ValueError:
                    continue
                rop = math.log2(i) + eta1 + math.log2((q - 1 / q**b) / (q - 1))
                if rop > eta:
                    continue  # Skip if complexity exceeds eta
                for st1 in range(1, 5):
                    for st2 in range(1, 5):
                        E.append(((i, eta1, st1), (j, eta2, st2), (1, 0, 'Discard-reduce')))

    # Transform-secret
    for eta1 in S:
        i = n
        if 2**eta1 > i:
            term = math.log2(2**eta1 - i)
            eta2 = Approx_S(term, S)
        else:
            continue
        rop = math.log2(i**2 + (2**eta1 - i) * i**2 / (math.log2(i) - math.log2(math.log2(i))))
        if rop > eta:
            continue  # Skip if complexity exceeds eta
        for st1 in range(1, 5):
            for st2 in range(1, 5):
                E.append(((i, eta1, st1), (j, eta2, st2),(1, 0, 'Transform-secret')))

    # Guess-reduce
    for i in range(1, n+1):
        for eta1 in S:
            if i - b > 0:
                j = i - b
                eta2 = Approx_S(eta1, S)
                rop = eta1 + math.log2(b) + b * math.log2(2 * d + 1)
                if rop > eta:
                    continue  # Skip if complexity exceeds eta
                for st1 in range(1, 5):
                    for st2 in range(1, 5):
                        E.append(((i, eta1, st1), (j, eta2, st2), (1, 0, 'Guess-reduce')))

    # Code-reduce
    for i in range(1, n+1):
        for eta1 in S:
            for j in range(1, i):
                eta2 = eta1
                c1 = 1
                c2 = sigma**2  # Assume sigma_cod^2 = sigma^2 here for simplicity
                rop = math.log2(i) + eta1
                if rop > eta:
                    continue  # Skip if complexity exceeds eta
                for st1 in range(1, 5):
                    for st2 in range(1, 5):
                        E.append(((i, eta1, st1), (j, eta2, st2), (c1, c2, 'Code-reduce')))

    return V, E









