#!/bin/sh
# Rebuild the three distinct documents without replacing the archived original.
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
python3 research/physical_rate_plots.py
python3 research/revise_report.py
python3 research/build_focused_appendices.py
python3 research/focused_overview.py
if command -v tectonic >/dev/null 2>&1; then
    task_tex=tectonic
elif [ -x tmp/pdfs/tectonic ]; then
    task_tex="$PWD/tmp/pdfs/tectonic"
else
    echo "Install Tectonic, or compile the three .tex files with a standard LaTeX engine." >&2
    exit 1
fi
mkdir -p output/pdf tmp/pdfs
for task_doc in focused_paper focused_appendices research_report; do
    "$task_tex" --keep-logs --outdir output/pdf "$task_doc.tex" >"tmp/pdfs/$task_doc-build.log" 2>&1
    echo "Built output/pdf/$task_doc.pdf"
done
