# AI index

## Core result

- **Natural-language statement:** [proof.pdf](proof.pdf), Theorems 1 and 5
- **LaTeX source:** [proof.tex](proof.tex)
- **Scope:** real positive-semidefinite 5 × 5 matrices, rank at most three
- **Sharp coefficient:** 13/48
- **Self-product equality:** signed simultaneous permutations of `J2 direct-sum E3`
- **Mixed equality:** after positive-diagonal normalization, a common 2 + 3 block partition with independent diagonal sign switches; zero-row cases are degenerate equalities

## Proof objects

- **Rational sum-of-squares data:** [sos_certificate.json](sos_certificate.json)
- **Covariant bases:** [certificate_basis0.json](certificate_basis0.json), [certificate_basis1.json](certificate_basis1.json), [certificate_basis2.json](certificate_basis2.json)
- **Rank-three ideal terms:** [relation_generators.json](relation_generators.json)
- **Equality implications:** [equality_certificate.json](equality_certificate.json)
- **Exact checker:** [verify.py](verify.py)
- **Corruption tests:** [test_mutations.py](test_mutations.py)
- **Mixed-extension regression:** [verify_mixed_extension.py](verify_mixed_extension.py)

## Trust boundaries

- `verify.py` checks the supplied finite identity and rational positivity data.
- `test_mutations.py` checks sensitivity to three proof-relevant corruptions.
- `verify_mixed_extension.py` checks the finite block-support consequences of the analytic mixed theorem.
- The written Gram reduction and equality deductions remain mathematical prose requiring inspection.
- No proof-assistant formalization, external peer review or exhaustive priority audit is claimed.
- Certificate construction is not reproducible because the generator and repair trace are absent.

## Entry points

- Human overview and commands: [README.md](README.md)
- Assurance boundary: [ASSURANCE.md](ASSURANCE.md)
- Machine scope: [STATUS.json](STATUS.json)
- Citation and related work: [proof.tex](proof.tex) and [proof.pdf](proof.pdf)
- File integrity: [SHA256SUMS](SHA256SUMS)
