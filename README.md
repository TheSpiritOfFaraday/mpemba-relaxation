# Mpemba relaxation beyond triplet tests

**Can a hotter preparation relax closer to equilibrium even when every triplet test fails?**
This research gives a reversible four-state counterexample and a sharp bound on
its asymptotic slow-mode amplitude advantage.

[Read the paper](output/pdf/focused_paper.pdf) ·
[Learn the ideas](output/pdf/freshman_guide.pdf) ·
[Explore the research](docs/RESEARCH_MAP.md) ·
[Reproduce a result](REPRODUCING.md)

## Choose a starting point

| I want to… | Start here | Next step |
| --- | --- | --- |
| Understand the main result | [Focused paper — PDF](output/pdf/focused_paper.pdf) | [Technical appendices](output/pdf/focused_appendices.pdf) |
| Learn the background step by step | [Beginner's guide — PDF](output/pdf/freshman_guide.pdf) | Download the repository and open [the interactive guide](freshman_guide.html) locally |
| See what was proved and what remains open | [Results and open questions](RESULTS.txt) | [Research map by topic](docs/RESEARCH_MAP.md) |
| Run a small exact check | [Quick reproduction instructions](REPRODUCING.md#small-exact-checks) | [Full command guide](RUNME.txt) |
| Inspect the formal proofs | [Lean proof guide](formal/README.md) | [Theorem axiom audit](formal/Audit.lean) |
| Find a particular file | [Browse by topic](docs/RESEARCH_MAP.md) | [Complete file catalog](docs/catalog/README.md) |

## The result, in brief

The ordinary-cooling triplet implication fails for a strictly positive reversible
four-state system. Under the full triplet-failure premise, the asymptotic
slow-amplitude advantage has a sharp supremum of **25%**. This is an amplitude
comparison; it does not measure a refrigeration gain or a universal reduction
in cooling time.

The project also studies finite classifications, higher-dimensional examples,
finite-time distances, and a separate idealized Arrhenius model. The
[research map](docs/RESEARCH_MAP.md) keeps each result beside its hypotheses,
source, and supporting checks. A later [finite strict witness above 8%](research/home_counterexample/README.md)
is included with its exact checker and separately labeled numerical screening.

## Read, inspect, or run

| Collection | What you will find |
| --- | --- |
| [Papers and guides](output/pdf/README.md) | Focused paper, appendices, comprehensive report, teaching guide |
| [Research source](research/README.md) | Topic routes to programs, mathematical fragments, and saved results |
| [Figures and data](output/figures/README.md) | Figure previews, vector PDFs, plotting data, and provenance |
| [Lean formalization](formal/README.md) | Pinned dependencies, theorem coverage, and build instructions |
| [Audit evidence](research/detailed_audit/README.md) | Historical review records and snapshots, with scope explained |
| [Reference literature](docs/REFERENCES.md) | Links to the cited third-party papers; their downloads are omitted |
| [Repository guide](docs/README.md) | Reading routes, the folder map, and the complete file catalog |

## A first check

From the repository root, with Python 3 installed:

```sh
python3 -I -S research/standalone_primary_check.py
python3 -I -S research/test_standalone_primary_check.py
```

These standard-library commands check the **primary counterexample** and its
rejection controls. The focused paper uses a different rational witness; its
checks and the rest of the reproduction routes are explained in
[REPRODUCING.md](REPRODUCING.md). Use a disposable copy for broader replay because
many research scripts regenerate their output files.

## How to read the evidence

The work includes analytic arguments, executable certificates, numerical studies,
and Lean proofs with different scopes. Lean covers the specified rational
companion, its eventual Pearson chi-square reversal, and the sharp four-state
amplitude supremum. Broader KL, physical-model, classification, and
higher-dimensional claims extend beyond that formal coverage. The work has not
received external human peer review.

[Recent package checks](PACKAGE_CHECKS.json) describe the bounded computations
rerun for this release. [Historical research records](research/detailed_audit/README.md)
remain historical evidence. Anonymization changes some file hashes without
changing mathematical content; see [the transformation record](ANONYMIZATION.json).

## Release information

[Reproduction](REPRODUCING.md) · [Privacy checks](PRIVACY_REPORT.json) ·
[Publication instructions](PUBLICATION.md) · [Rights and attribution](RIGHTS.md) ·
[File checksums](PACKAGE_MANIFEST.json)

This is an anonymous publication package with no existing Git history. The
publishing account and previously disclosed research can still create identity
links. No reuse license has been selected automatically.
