# Reproducing the research

Run commands from the package root unless stated otherwise. Use a disposable copy
for reproduction: many original programs overwrite generated JSON, logs, figures,
or manuscripts. Retain the downloaded package and its manifest as the baseline.

## Small exact checks

These checks need only Python's standard library:

```sh
python3 -I -S research/standalone_primary_check.py
python3 -I -S research/test_standalone_primary_check.py
python3 -I -S research/check_freshman_walkthrough.py
python3 -I -S research/test_freshman_walkthrough.py
```

The first pair checks the primary counterexample and rejection controls. The
second pair reconstructs the focused paper's different rational witness and its
negative controls; it writes receipts under `research/detailed_audit/`. These
checks do not establish every analytic claim in the manuscripts.

## Python environment

Python 3.11 is recommended for a new environment. Historical replay receipts used
Python 3.9.13. The five scientific libraries in `requirements.txt` are pinned to
versions documented in `RUNME.txt`; that is not a claim that every Python/platform
combination has been tested. Additional document-inspection dependencies are pinned
to the packaging environment's versions and identified separately in that file.

```sh
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
unset PYTHONOPTIMIZE
```

Do not use Python `-O` for the broader certificate suite: several programs rely on
assertions. The standalone checker above deliberately keeps its checks active in
optimized mode.

## Broader replay

```sh
python3 research/recheck_all.py --group core
python3 research/recheck_all.py --group extensions
python3 research/recheck_all.py --group physical
```

These groups execute proof-supporting computations and create timestamped audit
directories. The physical group can be lengthy; it replays the continuous
robustness obligations. Successful execution does not independently validate
every analytic reduction on which a certificate relies. Discovery optimizers
are separate and may follow different paths with different libraries or BLAS.

The later finite witness can be recomputed with:

```sh
python3 research/home_counterexample/certify_intermediate.py
python3 research/home_counterexample/check_intermediate.py
```

The certificate uses exact rational arithmetic and a stated spectral perturbation
bound; the screening script uses high-precision numerics. Neither is hardware
validation or a physical tolerance-box certificate.

## Lean

Install elan, Git, Python 3, and a POSIX shell. Initial dependency downloads require
network access. Preserve `formal/lean-toolchain`, `formal/lakefile.toml`, and
`formal/lake-manifest.json`; do not run `lake update` to reproduce this dependency
set.

```sh
cd formal
lake exe cache get
./check.sh
cd ..
```

The pinned versions are Lean 4.19.0 and mathlib revision
`c44e0c8ee63ca166450922a373c7409c5d26b00b`. The expected successful final line is:

```text
PASS: 43 theorem axiom audits; standard logical axioms only.
```

[formal/README.md](formal/README.md) documents the precise theorem scope and an
interpreter fallback for the mathlib cache downloader.

## Documents and historical checks

With the Python dependencies and Tectonic installed:

```sh
sh research/build_papers.sh
```

This regenerates the three manuscripts and PDFs. The interactive guide can be
rebuilt with `python3 research/build_freshman_guide.py`. Preserve the baseline
release first. Rebuilt PDFs and logs can acquire new metadata and local paths;
rerun the publication privacy checks before distributing regenerated files.

Some historical checks, including `research/check_publication_consistency.py`,
require original PDF hashes. Metadata sanitization changes those hashes. Such
historical checks are not unconditional acceptance tests for this release.
`ANONYMIZATION.json` records published hashes and transformation categories, not a
mapping to private originals. Do not replace historical expected hashes merely to
make a check pass or report historical checks as freshly passed.

Downloaded literature files are omitted and listed in `THIRD_PARTY_SOURCES.json`.
Checks that inspect those historical downloads require obtaining the cited
materials separately. The exact standalone checks above do not need them.

[RUNME.txt](RUNME.txt) retains the complete research command guide, hypotheses,
distinctions between numerical and exact evidence, and known open questions.
