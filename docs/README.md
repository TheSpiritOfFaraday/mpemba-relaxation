# Repository guide

[Home](../README.md) · [Research map](RESEARCH_MAP.md) · [All files](catalog/README.md)

Use these routes to follow the research without opening the full historical
archive. All original research files remain at their established paths so that
programs, manuscript inputs, and evidence links continue to work.

## Three reading routes

| Route | Suggested order |
| --- | --- |
| Learn the ideas | [Beginner's guide](../output/pdf/freshman_guide.pdf) → [focused paper](../output/pdf/focused_paper.pdf) → [results and open questions](../RESULTS.txt) |
| Review the mathematics | [Focused paper](../output/pdf/focused_paper.pdf) → [appendices](../output/pdf/focused_appendices.pdf) → [research map](RESEARCH_MAP.md) → [comprehensive report](../output/pdf/research_report.pdf) |
| Check a result yourself | [Small exact checks](../REPRODUCING.md#small-exact-checks) → [topic-specific checkers](RESEARCH_MAP.md) → [Lean guide](../formal/README.md) → [full reproduction commands](../RUNME.txt) |

The interactive [HTML teaching guide](../freshman_guide.html) is intended to be
opened locally after downloading the repository. Its [PDF version](../output/pdf/freshman_guide.pdf)
is the convenient reading copy when browsing files.

## Folder map

| Folder | Purpose | Human-readable guide |
| --- | --- | --- |
| Root | Manuscripts, principal results, reproduction and release instructions | [Home](../README.md) |
| `docs/` | Topic map, references, and navigation | This page |
| `research/` | Scientific programs, mathematical fragments, and saved outputs | [Research guide](../research/README.md) |
| `research/home_counterexample/` | Later finite strict counterexample | [Witness guide](../research/home_counterexample/README.md) |
| `research/detailed_audit/` | Historical reviews, receipts, and isolated snapshots | [Audit guide](../research/detailed_audit/README.md) |
| `research/archive/` | Earlier discovery artifacts | [Archived files](catalog/README.md#historical-discovery-archive) |
| `formal/` | Lean sources and pinned dependency manifests | [Formal proof guide](../formal/README.md) |
| `output/pdf/` | Exported papers and guides | [Document shelf](../output/pdf/README.md) |
| `output/figures/` | Figures, data tables, and plotting receipts | [Figure guide](../output/figures/README.md) |
| `output/playwright/` | Historical screenshots of the teaching guide | [Output guide](../output/README.md) |
| `tools/` | Package integrity, anonymous Git setup, and catalog generation | [Tool guide](../tools/README.md) |

## Find something specific

- **By result:** the [research map](RESEARCH_MAP.md) connects each topic to its
  explanation, implementation, and evidence.
- **By filename or file type:** the [file catalog](catalog/README.md) lists the
  package in short pages, with links and file sizes.
- **By cited paper:** [reference literature](REFERENCES.md) groups the source links
  used in the research.
- **By command:** [REPRODUCING.md](../REPRODUCING.md) gives the main routes;
  [RUNME.txt](../RUNME.txt) retains the detailed command guide.

## Evidence labels

| Label | How to interpret it |
| --- | --- |
| Exact executable check | A program checks specified arithmetic or certificate conditions; read its analytic assumptions |
| Lean formal proof | A theorem within the stated formal coverage, checked by the pinned Lean environment |
| Numerical study | Computational exploration or a plotted illustration; its precision does not make it an exact proof |
| Historical receipt | A record of an earlier run, bound to the historical inputs and hashes |
| Package check | A bounded check rerun while assembling this release; see [the recorded scope](../PACKAGE_CHECKS.json) |

The [privacy report](../PRIVACY_REPORT.json) concerns the publication package.
It does not extend the mathematical claims or establish external peer review.
