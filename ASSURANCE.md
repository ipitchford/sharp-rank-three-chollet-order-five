# Assurance statement

## What is established

The package contains a written reduction and an explicit rational certificate for the sharp real order-five rank-at-most-three self-Chollet inequality. The supplied standard-library checker verifies the complete polynomial identity, positivity of the three rational blocks, the four equality implications and the sharp witness. A separate implementation in the private review record reproduced the decisive exact calculations with a different algebraic organization.

The mixed two-matrix theorem and its equality classification are analytic consequences of the self-product theorem and Cauchy–Schwarz on permutation monomials. `verify_mixed_extension.py` supplies an exact semantic regression check for the block-support argument.

## What replay means

Successful replay establishes that the supplied finite rational data satisfy the identities checked by the programs in the recorded environment. The mutation suite shows that three plausible corruptions are rejected. Replay does not by itself establish historical priority, external peer review, independent reproduction, formal proof-assistant verification or correctness of arguments that appear only in prose.

## Important limits

- The certificate generator, numerical solver transcript and rational-repair trace are not included. The finished certificate can be checked, but its construction cannot be regenerated from this package.
- The 500-start numerical search is discovery evidence only; 79 starts carry successful termination flags.
- The novelty search is bounded and dated 19 October 2026.
- The theorem is real, order five and rank at most three. It does not assert the complex or unrestricted-rank Chollet conjecture.
- Internal AI-assisted review is not external specialist review or journal peer review.

## Reproduction

Run, from a fresh extraction:

```bash
python3 -B verify.py --report verification_report.json
python3 -B test_mutations.py
python3 -B verify_mixed_extension.py
PYTHONDONTWRITEBYTECODE=1 python3 -O -B verify.py --report verification_report_optimized.json
PYTHONDONTWRITEBYTECODE=1 python3 -O -B test_mutations.py
PYTHONDONTWRITEBYTECODE=1 python3 -O -B verify_mixed_extension.py
```

All proof-critical checks use exact integer or rational arithmetic.
