# Research files

[Project home](../README.md) · [Documentation](../docs/README.md) · [Research map](../docs/RESEARCH_MAP.md)

Start with the [focused paper](../output/pdf/focused_paper.pdf) for the argument,
then use this directory to inspect its calculations and supporting evidence.
The complete command guide is [RUNME.txt](../RUNME.txt); the shorter
[reproduction guide](../REPRODUCING.md) explains dependencies and safe replay.

| I want to… | Open |
| --- | --- |
| Check the primary rational counterexample | [Standalone exact checker](standalone_primary_check.py) and [rejection tests](test_standalone_primary_check.py) |
| Check the focused paper's different rational witness | [Walkthrough checker](check_freshman_walkthrough.py) and [negative controls](test_freshman_walkthrough.py) |
| Explore the finite witness with an amplitude enclosure above 8% | [Intermediate-counterexample guide](home_counterexample/README.md) |
| Follow the quarter-limit construction | [Certificate program](quarter_limit_certificate.py) and [saved certificate](quarter_limit_certificate.json) |
| Replay groups of checks | [Replay runner](recheck_all.py) and [reproduction instructions](../REPRODUCING.md) |
| Inspect the paper build | [Build script](build_papers.sh) |
| Follow the teaching-guide build | [Guide builder](build_freshman_guide.py) and [finished guide](../freshman_guide.html) |
| Examine research-era audits and snapshots | [Audit guide](detailed_audit/README.md) |
| Read formal proof coverage | [Lean guide](../formal/README.md) |
| Browse figures and their data | [Figure guide](../output/figures/README.md) |

Scripts, saved outputs, exploratory searches, and historical drafts coexist here.
A saved result is evidence of its recorded run, not a fresh execution. The
publication's bounded fresh checks are in [PACKAGE_CHECKS.json](../PACKAGE_CHECKS.json);
its privacy audit is in [PRIVACY_REPORT.json](../PRIVACY_REPORT.json). Neither
establishes every analytic claim or replaces external peer review.

Run research programs in a disposable copy: many overwrite receipts, plots, or
manuscripts. Preserve the original package and its manifest first.
