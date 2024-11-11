import argparse
import itertools
import warnings
from pysat.solvers import Cadical195


def find_non_sibo_matroid(r):
    solver = Cadical195()
    warnings.simplefilter("ignore", SyntaxWarning)

    n = 2 * r
    E = list(range(n))
    id_to_subset = {
        i + 1: frozenset(subset)
        for i, subset in enumerate(itertools.combinations(E, r))
    }
    subset_to_id = {subset: id for id, subset in id_to_subset.items()}

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

    # A and B are bases with no SI-ordering
    A = list(range(r))
    B = list(range(r, n))
    solver.add_clause([subset_to_id[frozenset(A)]])
    solver.add_clause([subset_to_id[frozenset(B)]])
    for pA in itertools.permutations(range(r)):
        for pB in itertools.permutations(range(r)):
            clause = [
                -subset_to_id[
                    frozenset(
                        A[pA[k]] if (k <= i or k >= j + 1) else B[pB[k]]
                        for k in range(r)
                    )
                ]
                for i in range(-1, r)
                for j in range(i, r)
            ]
            solver.add_clause(clause)

    b = solver.solve()
    if b:
        while b:
            ans = solver.get_model()
            print(
                f"No SI-ordering for the basis pair {A}, {B} in the matroid with the following bases:",
                " ".join(str(set(id_to_subset[x])) for x in ans if x > 0),
            )
            # excluding all solutions isomorphic to the found one
            for pA in itertools.permutations(A):
                for pB in itertools.permutations(B):
                    clause = []
                    for X, id_X in subset_to_id.items():
                        Xp = {pA[x] if x in A else pB[x - r] for x in X}
                        clause.append(
                            -subset_to_id[frozenset(Xp)]
                            if ans[id_X - 1] > 0
                            else subset_to_id[frozenset(Xp)]
                        )
                    solver.add_clause(clause)
            b = solver.solve()
        print("SI-ordering exists for every other basis pair")
    else:
        print(f"Each matroid of rank {r} is SIBO")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--rank", type=int, required=True, help="Rank of matroid")
    args = parser.parse_args()

    find_non_sibo_matroid(args.rank)
