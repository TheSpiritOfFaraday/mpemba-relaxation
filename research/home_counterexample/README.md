# A finite strict counterexample above 8%

[Project home](../../README.md) · [Research files](../README.md) · [Research map](../../docs/RESEARCH_MAP.md)

This later example specializes the quarter-limit family to a finite reversible
four-state generator. Its exact certificate verifies strict triplet failure and
encloses the asymptotic slow-mode amplitude advantage between approximately
8.1718% and 8.1823%. This is a mathematical amplitude comparison, not a certified
cooling-time percentage or a hardware result.

| File | Role |
| --- | --- |
| [certify_intermediate.py](certify_intermediate.py) | Exact rational acceptance checks with a stated spectral perturbation bound |
| [intermediate_certificate.json](intermediate_certificate.json) | Saved generator, strict triplet residuals, gap bounds, and amplitude enclosure |
| [check_intermediate.py](check_intermediate.py) | Exact structural checks plus 160-digit numerical thermal screening |
| [intermediate_results.json](intermediate_results.json) | Saved screening results, including numerical finite-time and crossing estimates |

The screening results are not interval certificates for crossing time, uniqueness,
or physical tolerance. The exact certificate and the numerical screen make
different claims. Both use the archived
[quarter-limit family](../quarter_limit_certificate.json).

From the package root, in a disposable copy with the documented dependencies:

```sh
python3 research/home_counterexample/certify_intermediate.py
python3 research/home_counterexample/check_intermediate.py
```

These commands overwrite their corresponding JSON receipts. The certificate also
imports a helper that requires `mpmath`. See [REPRODUCING.md](../../REPRODUCING.md).

[PACKAGE_CHECKS.json](../../PACKAGE_CHECKS.json) records a fresh successful replay
of the exact certificate with unchanged receipt bytes. It does not record a new
full screening run or physical validation.
