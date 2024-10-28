"""
Python script for finding a non-SIBO matroid
"""

import itertools
import warnings
from pysat.solvers import Cadical195

warnings.simplefilter("ignore", SyntaxWarning)

solver = Cadical195()

r = 5  # matroid rank
n = 2 * r
excludeR10 = True

E = list(range(n))
id_to_subset = {
    i + 1: frozenset(subset) for i, subset in enumerate(itertools.combinations(E, r))
}
subset_to_id = {subset: id for id, subset in id_to_subset.items()}


# R_{10} from K5 representation: a 5-cycle and its complement are labeled by A and B, respectively
def R10(A, B):
    is_base = {id: True for id in id_to_subset.keys()}
    edge_name = [[0] * 5 for _ in range(5)]
    for i in range(5):
        edge_name[i][(i + 1) % 5] = edge_name[(i + 1) % 5][i] = A[i]
        edge_name[i][(i + 2) % 5] = edge_name[(i + 2) % 5][i] = B[i]

    for v in range(5):
        V = [i for i in range(5) if i != v]
        cycles4 = [
            [V[0], V[1], V[2], V[3]],
            [V[0], V[1], V[3], V[2]],
            [V[0], V[2], V[1], V[3]],
        ]
        for cycle in cycles4:
            C_set = set()
            for i in range(4):
                C_set.add(edge_name[cycle[i]][cycle[(i + 1) % 4]])
            for e in range(n):
                if e not in C_set:
                    S = C_set.copy()
                    S.add(e)
                    is_base[subset_to_id[frozenset(S)]] = False
    return is_base


# basis exchange property
for A, id_A in subset_to_id.items():
    for B, id_B in subset_to_id.items():
        for e in A.difference(B):
            clause = [-id_A, -id_B]
            for f in B.difference(A):
                Aef = set(A)
                Aef.remove(e)
                Aef.add(f)
                clause.append(subset_to_id[frozenset(Aef)])
            solver.add_clause(clause)

# A and B are bases such that no ordering is a SIBO
A = list(range(r))
B = list(range(r, n))
solver.add_clause([subset_to_id[frozenset(A)]])
solver.add_clause([subset_to_id[frozenset(B)]])
for pA in itertools.permutations(range(r)):
    for pB in itertools.permutations(range(r)):
        clause = []
        for i in range(-1, r):
            for j in range(i, r):
                act = set()
                for k in range(r):
                    if k <= i or k >= j + 1:
                        act.add(A[pA[k]])
                    else:
                        act.add(B[pB[k]])
                clause.append(-subset_to_id[frozenset(act)])
        solver.add_clause(clause)

if r == 5 and excludeR10:
    print("expect to wait a few minutes")
    for pA in itertools.permutations(A):
        for pB in itertools.permutations(B):
            M = R10(pA, pB)
            clause = [-id if is_base else id for id, is_base in M.items()]
            solver.add_clause(clause)


b = solver.solve()
if b:
    ans = solver.get_model()
    print(
        f"No SIBO for the basis pair {A}, {B} in the matroid with the following bases:",
    )
    print(" ".join(str(set(id_to_subset[x])) for x in ans if x > 0))
else:
    print("SIBO exists for each basis pair")
