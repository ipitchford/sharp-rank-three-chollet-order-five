#!/usr/bin/env python3
"""Independent exact checker for the rank-three self-Chollet certificate.

Python standard library only. No search routine, symbolic-algebra package,
or producer module is imported. Every proof-critical calculation is exact;
elapsed-time metadata alone uses a clock. The checker:
  (1) proves each rational Gram matrix positive definite using a rational
      Cholesky-like factor and strict diagonal dominance;
  (2) reconstructs the scalar/vector/tensor polynomial coefficients;
  (3) expands the full S_5 average of the residual, checking every coefficient;
  (4) verifies the exact sharpness witness.

Run: python3 verify.py
"""
from fractions import Fraction
from itertools import combinations_with_replacement, permutations
from pathlib import Path
from math import gcd, lcm
import argparse
import hashlib
import json
import time

PAIRS = tuple(combinations_with_replacement(range(5), 2))
PAIR_INDEX = {p: i for i, p in enumerate(PAIRS)}
PERMS = tuple(permutations(range(5)))
ZERO = (0,) * 15


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def add_term(poly, monomial, coefficient):
    if coefficient:
        poly[monomial] = poly.get(monomial, 0) + coefficient
        if not poly[monomial]:
            del poly[monomial]


def multiply_monomials(a, b):
    return tuple(x + y for x, y in zip(a, b))


def with_edges(m, edges):
    a = list(m)
    for i, j in edges:
        a[PAIR_INDEX[tuple(sorted((i, j)))]] += 1
    return tuple(a)


def target_polynomial():
    p = {}
    for sigma in PERMS:
        a = [0] * 15
        for i, j in enumerate(sigma):
            a[PAIR_INDEX[tuple(sorted((i, j)))]] += 1
        add_term(p, tuple(a), 1)
    result = {}
    for a, ca in p.items():
        add_term(result, tuple(2*x for x in a), -48*ca)
        for b, cb in p.items():
            add_term(result, multiply_monomials(a, b), 13*ca*cb)
    return result


def determinant_times_monomial(spec):
    rows, cols = spec['rows'], spec['cols']
    m = tuple(spec['monomial'])
    require(len(rows) == len(cols) == 4, 'A determinant is not 4 by 4.')
    require(len(set(rows)) == len(set(cols)) == 4, 'Repeated determinant index.')
    require(all(0 <= i < 5 for i in rows + cols), 'Determinant index out of range.')
    require(len(m) == 15 and all(isinstance(v, int) and v >= 0 for v in m), 'Invalid monomial.')
    poly = {}
    for sigma in permutations(range(4)):
        inversions = sum(sigma[i] > sigma[j] for i in range(4) for j in range(i+1, 4))
        a = with_edges(m, [(rows[i], cols[sigma[i]]) for i in range(4)])
        add_term(poly, a, -1 if inversions % 2 else 1)
    return poly


def check_positive_definite(block, qden):
    q = [[int(x) for x in row] for row in block['Q_numerators']]
    ell = [[int(x) for x in row] for row in block['L_numerators']]
    ld = int(block['L_denominator'])
    n = block['dimension']
    require(qden > 0 and ld > 0, 'Nonpositive denominator.')
    require(len(q) == len(ell) == n, 'Incorrect Gram dimension.')
    require(all(len(row) == n for row in q + ell), 'Nonsquare Gram matrix.')
    require(all(q[i][j] == q[j][i] for i in range(n) for j in range(n)), 'Nonsymmetric Q.')
    # W/(qden*ld^2) = Q - L L^T. Compute every entry in integers.
    diagonal = [0] * n
    off_sum = [0] * n
    for i in range(n):
        for j in range(i, n):
            product = sum(ell[i][k] * ell[j][k] for k in range(n))
            w = q[i][j] * ld * ld - qden * product
            if i == j:
                diagonal[i] = w
            else:
                off_sum[i] += abs(w)
                off_sum[j] += abs(w)
    margins = [diagonal[i] - off_sum[i] for i in range(n)]
    denominator = qden * ld * ld
    require(all(5*x > denominator for x in margins), 'Q >= (1/5)I not certified.')
    return q, Fraction(min(margins), denominator)


def validate_basis(data, dimension):
    spin = data['spin']
    require(spin in (0, 1, 2), 'Invalid covariant type.')
    require(data['coordinate_order'] == [list(p) for p in PAIRS], 'Coordinate order mismatch.')
    basis = data['basis']
    n = data['N']
    require(len(n) == len(basis), 'Basis/change-of-basis size mismatch.')
    require(all(len(row) == dimension for row in n), 'Incorrect covariant dimension.')
    require(all(isinstance(x, int) for row in n for x in row), 'Noninteger basis coefficient.')
    for item in basis:
        ij, a = item['pair'], item['powers']
        require(len(a) == 15 and all(isinstance(x, int) and x >= 0 for x in a), 'Invalid basis monomial.')
        degree = [0] * 5
        for (i, j), x in zip(PAIRS, a):
            degree[i] += x
            degree[j] += x
        if spin == 0:
            require(ij == [], 'Scalar basis has vector indices.')
        else:
            require(len(ij) == 2 and 0 <= ij[0] <= ij[1] < 5, 'Invalid covariant pair.')
            require(spin != 1 or ij[0] < ij[1], 'Invalid cross-product pair.')
            degree[ij[0]] += 1
            degree[ij[1]] += 1
        require(degree == [2]*5, 'Covariant has incorrect multidegree.')
    den = data['denominator']
    require(isinstance(den, int) and den > 0, 'Invalid basis denominator.')
    return basis, n, den


def expand_block(block, basis_data, q, common_basis_denominator):
    """Return 6*D^2*qden times the block's unaveraged SOS polynomial."""
    spin = block['spin']
    require(spin == basis_data['spin'], 'Spin mismatch.')
    basis, nmat, den = validate_basis(basis_data, len(q))
    require(common_basis_denominator % den == 0, 'Basis denominators incompatible.')
    scale = (common_basis_denominator // den)**2
    sparse = [[(k, x) for k, x in enumerate(row) if x] for row in nmat]
    # T = N Q_num. A = T N^T is used one entry at a time.
    t = [[sum(x*q[k][j] for k, x in row) for j in range(len(q))] for row in sparse]
    poly = {}
    for u, left in enumerate(basis):
        for v in range(u, len(basis)):
            a = sum(t[u][k]*x for k, x in sparse[v])
            if not a:
                continue
            a *= scale * (1 if u == v else 2)
            right = basis[v]
            m = multiply_monomials(left['powers'], right['powers'])
            if spin == 0:
                add_term(poly, m, 6*a)
            else:
                i, j = left['pair']
                k, ell = right['pair']
                first = with_edges(m, [(i, k), (j, ell)])
                second = with_edges(m, [(i, ell), (j, k)])
                if spin == 1:
                    add_term(poly, first, 6*a)
                    add_term(poly, second, -6*a)
                else:
                    add_term(poly, first, 3*a)
                    add_term(poly, second, 3*a)
                    add_term(poly, with_edges(m, [(i, j), (k, ell)]), -2*a)
    return poly


def permutation_average_numerator(poly):
    """Expand the sum over all 120 relabellings, without using orbit reduction."""
    total = {}
    maps = []
    for sigma in PERMS:
        maps.append([PAIR_INDEX[tuple(sorted((sigma[i], sigma[j])))] for i, j in PAIRS])
    for m, coefficient in poly.items():
        for relabel in maps:
            a = [0] * 15
            for old, new in enumerate(relabel):
                a[new] = m[old]
            add_term(total, tuple(a), coefficient)
    return total



def check_equality_implications(root, scalar_basis):
    """Check four exact scalar identities used in the equality classification."""
    data = json.loads((root/'equality_certificate.json').read_text())
    require(data['format'] == 'rank3-self-chollet-equality-v1', 'Unknown equality certificate.')
    from itertools import combinations
    def triangle(t):
        other = [i for i in range(5) if i not in t]
        return with_edges(ZERO, [(t[0],t[1]),(t[1],t[2]),(t[2],t[0])] + [(i,i) for i in other])
    expected = []
    p = {with_edges(ZERO, [(i,i) for i in range(5)]): 1}
    for t in combinations(range(5), 3):
        add_term(p, triangle(t), 8)
    expected.append(p)
    expected.append({triangle((2,3,4)): 1,
                     with_edges(ZERO, [(0,1),(0,1),(2,3),(3,4),(4,2)]): -1})
    expected.append({with_edges(ZERO, [(0,1),(1,2),(2,3),(3,0),(4,4)]): 1})
    p = {with_edges(ZERO, [(0,1),(0,1),(2,2),(3,3),(4,4)]): 1}
    add_term(p, triangle((2,3,4)), 8)
    for k in (2,3,4):
        add_term(p, triangle((0,1,k)), 2)
    expected.append(p)
    basis = scalar_basis['basis']
    nmat = scalar_basis['N']
    dimension = len(nmat[0])
    generators = []
    for j in range(dimension):
        poly = {}
        for item, row in zip(basis, nmat):
            add_term(poly, tuple(item['powers']), row[j])
        generators.append(poly)
    generators.extend(determinant_times_monomial(spec) for spec in data['relations'])
    require(len(data['identities']) == 4, 'Expected four equality implications.')
    for i, (identity, target) in enumerate(zip(data['identities'], expected)):
        require(identity['index'] == i, 'Equality identities out of order.')
        recorded = {tuple(term['powers']): term['coefficient'] for term in identity['target']}
        require(recorded == target, 'Equality target does not match the fixed statement.')
        coefficients = [Fraction(int(a), int(b)) for a,b in identity['coefficients']]
        require(len(coefficients) == len(generators), 'Incorrect equality coefficient count.')
        den = lcm(*(c.denominator for c in coefficients))
        residual = {m: -c*den for m,c in target.items()}
        for a, poly in zip(coefficients, generators):
            n = a.numerator*(den//a.denominator)
            for m,c in poly.items():
                add_term(residual, m, n*c)
        require(not residual, f'Equality implication {i} failed.')
    return 4

def permanent(a):
    total = Fraction(0)
    for sigma in permutations(range(len(a))):
        product = Fraction(1)
        for i, j in enumerate(sigma):
            product *= a[i][j]
        total += product
    return total


def check_sharpness():
    h = Fraction(1, 2)
    c = [[Fraction(x) for x in row] for row in
         [[1,1,0,0,0],[1,1,0,0,0],[0,0,1,-h,-h],
          [0,0,-h,1,-h],[0,0,-h,-h,1]]]
    # Rational Gram coordinates with a positive diagonal inner product.
    v = [[0,0,1],[0,0,1],[1,0,0],[-h,h,0],[-h,-h,0]]
    metric = [1,3,1]
    require(all(c[i][j] == sum(v[i][k]*metric[k]*v[j][k] for k in range(3))
                for i in range(5) for j in range(5)), 'Sharpness Gram factorisation failed.')
    # The rows 0,2,3 span a three-dimensional space.
    determinant = v[0][2] * (v[2][0]*v[3][1] - v[2][1]*v[3][0])
    require(determinant != 0, 'Sharpness rank is not three.')
    p, q = permanent(c), permanent([[x*x for x in row] for row in c])
    require(p == 3 and q == Fraction(39,16), 'Sharpness permanent calculation failed.')
    require(q/(p*p) == Fraction(13,48), 'Sharpness ratio failed.')
    return {'permanent': str(p), 'square_permanent': str(q), 'ratio': str(q/(p*p)), 'rank': 3}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--report', type=Path, default=None)
    args = parser.parse_args()
    root = args.directory
    start = time.time()
    cert = json.loads((root/'sos_certificate.json').read_text())
    require(cert['format'] == 'rank3-self-chollet-sos-v1', 'Unknown certificate format.')
    require(cert['constant'] == [13,48], 'Wrong target constant.')
    require(len(cert['blocks']) == 3, 'Expected three Gram blocks.')
    qden = int(cert['q_denominator'])
    bases = [json.loads((root/f'certificate_basis{s}.json').read_text()) for s in range(3)]
    common = lcm(*(b['denominator'] for b in bases))
    polynomial_denominator = 6*common*common*qden
    raw = {}
    margins = []
    for s, (block, basis) in enumerate(zip(cert['blocks'], bases)):
        require(block['spin'] == s, 'Gram blocks out of order.')
        q, margin = check_positive_definite(block, qden)
        margins.append(str(margin))
        print(f'PASS: rational Gram block {s}, dimension {len(q)}, lower bound > 1/5', flush=True)
        for m, coefficient in expand_block(block, basis, q, common).items():
            add_term(raw, m, coefficient)
    # Clear the determinant-term denominator exactly.
    ld = int(cert['lambda_denominator'])
    require(ld > 0, 'Invalid determinant coefficient denominator.')
    raw = {m: x*ld for m, x in raw.items()}
    relation_specs = json.loads((root/'relation_generators.json').read_text())
    numerators = [int(x) for x in cert['lambda_numerators']]
    require(len(relation_specs) == len(numerators) == 25, 'Wrong number of determinant terms.')
    for spec, numerator in zip(relation_specs, numerators):
        for m, coefficient in determinant_times_monomial(spec).items():
            add_term(raw, m, polynomial_denominator*numerator*coefficient)
    for m, coefficient in target_polynomial().items():
        add_term(raw, m, -polynomial_denominator*ld*coefficient)
    unaveraged_count = len(raw)
    residual = permutation_average_numerator(raw)
    require(not residual, f'Polynomial identity failed: {len(residual)} nonzero monomials.')
    print('PASS: full 120-permutation polynomial identity; zero nonzero residual coefficients', flush=True)
    equality_implications = check_equality_implications(root, bases[0])
    print('PASS: four exact scalar implications for the complete equality classification', flush=True)
    sharpness = check_sharpness()
    print('PASS: exact positive-semidefinite rank-three equality witness, ratio 13/48', flush=True)
    report = {'status':'PASS', 'arithmetic':'integer and rational; Python standard library only',
              'gram_dimensions':[b['dimension'] for b in cert['blocks']],
              'gram_lower_bound':'1/5', 'exact_gershgorin_margins':margins,
              'determinant_terms':len(relation_specs), 'permutations':120,
              'unaveraged_residual_monomials':unaveraged_count,
              'nonzero_averaged_residual_monomials':0, 'equality_implications':equality_implications, 'sharpness':sharpness,
              'seconds':time.time()-start,
              'certificate_sha256':hashlib.sha256((root/'sos_certificate.json').read_bytes()).hexdigest(),
              'input_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in
                ['sos_certificate.json','certificate_basis0.json','certificate_basis1.json',
                 'certificate_basis2.json','relation_generators.json','equality_certificate.json','verify.py']}}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+'\n')
    print('ALL EXACT CHECKS PASSED', flush=True)


if __name__ == '__main__':
    main()
