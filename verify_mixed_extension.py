#!/usr/bin/env python3
"""Exact semantic checks for the mixed rank-three equality classification.

This is a regression test for Theorem 5, not a substitute for its proof.
All arithmetic uses fractions.
"""

from fractions import Fraction
from itertools import combinations, permutations, product


N = 5
PERMS = tuple(permutations(range(N)))


def c_star(pair, signs=(1, 1, 1, 1, 1)):
    pair = frozenset(pair)
    triple = frozenset(range(N)) - pair
    matrix = [[Fraction(0) for _ in range(N)] for _ in range(N)]
    for i in range(N):
        matrix[i][i] = Fraction(1)
    a, b = sorted(pair)
    matrix[a][b] = matrix[b][a] = Fraction(1)
    for i, j in combinations(sorted(triple), 2):
        matrix[i][j] = matrix[j][i] = Fraction(-1, 2)
    return [
        [Fraction(signs[i] * signs[j]) * matrix[i][j] for j in range(N)]
        for i in range(N)
    ]


def monomial_vector(matrix):
    return tuple(
        product_fraction(matrix[i][sigma[i]] for i in range(N))
        for sigma in PERMS
    )


def product_fraction(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def permanent(matrix):
    return sum(monomial_vector(matrix), Fraction(0))


def hadamard(a, b):
    return [[a[i][j] * b[i][j] for j in range(N)] for i in range(N)]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    pairs = tuple(combinations(range(N), 2))
    vectors = {pair: monomial_vector(c_star(pair)) for pair in pairs}

    require(
        len(set(vectors.values())) == len(pairs),
        "distinct 2+3 partitions must have distinct monomial vectors",
    )
    require(
        all(sum(value != 0 for value in vector) == 12 for vector in vectors.values()),
        "every block support must contain exactly 12 permutations",
    )

    for pair in pairs:
        base = vectors[pair]
        for signs in product((-1, 1), repeat=N):
            require(
                monomial_vector(c_star(pair, signs)) == base,
                "diagonal sign switching changed a permutation-monomial vector",
            )

    constant = Fraction(13, 48)
    for pair in pairs:
        a = c_star(pair)
        require(
            permanent(hadamard(a, a)) == constant * permanent(a) ** 2,
            "a canonical equality matrix did not attain 13/48",
        )

    first, second = pairs[0], pairs[1]
    a, b = c_star(first), c_star(second)
    require(
        permanent(hadamard(a, b)) < constant * permanent(a) * permanent(b),
        "a mismatched block partition must not saturate the mixed inequality",
    )

    print("PASS: 10 distinct block supports, sign invariance, sharp matched pairs")
    print("PASS: mismatched block partition rejected as an equality case")


if __name__ == "__main__":
    main()
