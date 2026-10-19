# Sharp rank-three Chollet inequalities in order five

**Creator:** Anonymous  
**Date:** 19 October 2026  
**Status:** exact computer-assisted proof, with complete self-product and mixed equality classifications; internally replayed, not externally peer reviewed or proof-assistant formalised.

## Result

For every real symmetric positive-semidefinite 5 × 5 matrix C with diagonal one and rank at most three,

    48 per(C ∘ C) ≤ 13 per(C)².

Equality holds exactly at signed simultaneous permutations of C★ = J₂ ⊕ E₃, where E₃ has diagonal one and off-diagonal entries −1/2. Consequently the sharp constant is 13/48.

For real positive-semidefinite 5 × 5 matrices A and B, each of rank at most three,

    per(A ∘ B) ≤ (13/48) per(A) per(B).

This constant is also sharp. With positive diagonals, equality holds exactly when the two diagonal normalizations have the same 2 + 3 block partition and each is a diagonal sign switching of the corresponding copy of C★. If either matrix has a zero diagonal entry, its zero row makes both sides vanish.

The theorem includes all rank drops. It is not a claim about complex correlation matrices or unrestricted rank.

## Read and verify

Read `proof.pdf` for the theorem, the certificate reduction, and the complete equality argument. The finite rational files are integral to the proof: the PDF alone does not reproduce their matrix entries.

Run in this directory:

```bash
python3 verify.py --report verification_report.json
python3 test_mutations.py
python3 verify_mixed_extension.py
```

The exact checker uses the Python standard library only. Its proof-critical operations are integer and rational arithmetic, with no numerical tolerance, optimiser, external data, or producer code. Python 3.9 or later supplies all imported standard-library functions. The recorded run took approximately 1.2 seconds; runtime is environment-dependent.

Expected final line:

```text
ALL EXACT CHECKS PASSED
```

## What the checker establishes

The checker reconstructs three covariant sum-of-squares expressions, with positive-definite rational Gram matrices of dimensions 52, 96 and 140. Rational factorisations and strict diagonal dominance certify Q ≥ (1/5)I for each matrix. It then expands a full 120-permutation polynomial identity. The difference between the target and the sum of squares is an explicit sum of monomial multiples of 4 × 4 Gram minors. Those minors vanish for rank at most three.

The checker also verifies four exact scalar identities used in the equality classification, and the rational witness values per(C★) = 3 and per(C★ ∘ C★) = 39/16. The written reduction from the four identities to the complete classification appears in Section 4 of the proof.

`test_mutations.py` deliberately corrupts a positivity certificate, a coefficient in the global identity, and an equality identity. All three corruptions must be rejected.

`verify_mixed_extension.py` is a supplementary semantic regression test for the analytic mixed theorem. It checks all ten 2 + 3 block partitions, all diagonal sign switches, and a deliberately mismatched partition. It is not a substitute for the proof in the manuscript.

## Files

- `proof.pdf`, `proof.tex`: the mathematical argument and its source.
- `verify.py`: the separate, standard-library-only exact checker.
- `verify_mixed_extension.py`: exact regression checks for the mixed equality classification.
- `sos_certificate.json`: the three rational Gram matrices, rational positivity certificates, and determinant coefficients.
- `certificate_basis0.json`, `certificate_basis1.json`, `certificate_basis2.json`: explicit scalar, vector and traceless-tensor covariants and integer basis maps.
- `relation_generators.json`: 25 monomial multiples of 4 × 4 determinants.
- `equality_certificate.json`: exact scalar consequences that classify equality.
- `verification_report.json`, `mutation_test_report.json`: results from actual executions.
- `STATUS.json`, `SHA256SUMS`: scope, verification status and file integrity.
- `discovery/`: a numerical multistart search and its recorded output, neither used in the proof.

The `discovery/search.py` script additionally needs NumPy, SciPy and Numba. Certificate verification does not.

## Provenance and limitations

The initial research target is equation (7) in Zeng, Liu and Ratnavelu's August 2026 preprint. The preprint supplied the candidate constant and extremiser; the global certificate here does not assume the correctness of its structured-case or local proofs.

A numerical invariant sum-of-squares search supplied a candidate, followed by rational reconstruction and exact repair. The package does not contain the generator, solver transcript or repair trace, so it does not reproduce certificate construction. Numerical discovery and search coverage are not evidence of global validity; the exact polynomial identity and rational positivity checks carry the proof.

The retained discovery log contains 500 starts, of which 79 carry a successful termination flag. These runs are not proof evidence.

The public-literature check was targeted, not exhaustive. No independent external mathematical review or proof-assistant formalisation has been completed. No claim is made here about the general Chollet conjecture.

## Sources

Zeng, Z., Liu, H., & Ratnavelu, K. (2026, August 18). *Sharp self-Chollet inequalities in order five: Rank two and structured rank three* [Preprint, version 1]. Preprints.org. https://doi.org/10.20944/preprints202608.1240.v1

Li, Q., Zhao, Q., & Luo, W. (2026, August 31). *Chollet's permanent conjecture through order six via border induction* [Preprint, version 1]. Preprints.org. https://www.preprints.org/manuscript/202608.2196

Gatermann, K., & Parrilo, P. A. (2004). Symmetry groups, semidefinite programs, and sums of squares. *Journal of Pure and Applied Algebra, 192*(1–3), 95–128. https://doi.org/10.1016/j.jpaa.2003.12.011

Peyrl, H., & Parrilo, P. A. (2008). Computing sum of squares decompositions with rational coefficients. *Theoretical Computer Science, 409*(2), 269–281. https://doi.org/10.1016/j.tcs.2008.09.025
