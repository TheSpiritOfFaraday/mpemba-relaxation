# Historical audit records

[Project home](../../README.md) · [Research files](../README.md) · [Documentation](../../docs/README.md)

This directory preserves research-era reviews, command receipts, logs, visual
checks, and isolated copies used for earlier rebuilds. Date-stamped directories
describe those past runs. Their presence is not a claim that all checks were
rerun when this GitHub package was prepared.

| Start with | What it contains |
| --- | --- |
| [Final research review](final_completion.json) | Scope, requirements, evidence pointers, and limitations of the recorded review |
| [Claim reviews](claim_reviews.json) | Review of mathematical statements and dependencies |
| [Command adjudication](command_adjudication.json) | Recorded command outcomes and exceptions |
| [Source correspondence](source_correspondence.json) | Historical correspondence between claims and cited sources |
| [Core replay](20261013T011946Z_core/) | Recorded core-check run |
| [Extension replay](20261013T012015Z_extensions/) | Recorded extension-check run |
| [Physical replay](20261013T021448Z_physical/) | Recorded physical-model checks |
| [Formal-command receipt](20261013_formal_command/receipt.json) | Historical Lean command evidence; see the [formal coverage guide](../../formal/README.md) |
| [Final document review](final_delivery/verification.json) | Recorded document and teaching-guide verification |

For checks freshly performed while preparing this release, use the root-level
[PACKAGE_CHECKS.json](../../PACKAGE_CHECKS.json) and
[PRIVACY_REPORT.json](../../PRIVACY_REPORT.json). The first reports bounded
scientific smoke checks; the second reports privacy verification. They do not
replace the analytic arguments or external peer review.

Historical hashes bind historical bytes. PDF metadata removal and other privacy
transformations can change a hash without changing the underlying mathematics.
The publication's [ANONYMIZATION.json](../../ANONYMIZATION.json) explains this
distinction. Downloaded third-party full text is omitted; retrieval references
remain in [THIRD_PARTY_SOURCES.json](../../THIRD_PARTY_SOURCES.json).

Some `isolated/` directories contain old copies of papers and programs. Use the
top-level manuscripts and [research entry points](../README.md) for the assembled
release; retain these snapshots as historical evidence.
