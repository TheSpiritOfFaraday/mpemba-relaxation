#!/bin/sh
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
if [ -x .tools/elan/bin/lake ]; then
  ELAN_HOME="$PWD/.tools/elan"
  export ELAN_HOME
  task_lake="$PWD/.tools/elan/bin/lake"
else
  task_lake=lake
fi
"$task_lake" build
task_audit="$(mktemp)"
trap 'test ! -e "$task_audit" || rm -- "$task_audit"' EXIT HUP INT TERM
"$task_lake" env lean -DwarningAsError=true Audit.lean > "$task_audit"
cat "$task_audit"
python3 - "$task_audit" <<'PY'
import re
import sys
from pathlib import Path

output = Path(sys.argv[1]).read_text()
audits = re.findall(r"'([^']+)' depends on axioms: \[([^\]]*)\]", output)
if len(audits) != 43:
    raise SystemExit(f"Expected 43 audited theorems; found {len(audits)}")
declared = re.findall(r"^#print axioms (\S+)$", Path("Audit.lean").read_text(), re.M)
reported = [theorem for theorem, _ in audits]
if len(set(declared)) != 43 or set(reported) != set(declared):
    raise SystemExit("The axiom reports do not match the 43 distinct declared audits")
required = {
    "Mpemba.appendixB_counterexample",
    "Mpemba.overlap_implies_eventual_distance_order",
    "Mpemba.master_equation_distance_reversal",
    "Mpemba.Quarter.sharp_quarter_isLUB",
    "Mpemba.Quarter.sharp_quarter_not_attained",
    "Mpemba.Quarter.Sharp.eventually_family_attainable",
}
if not required <= set(reported):
    raise SystemExit(f"Missing requested-result audits: {sorted(required-set(reported))}")
allowed = {"propext", "Classical.choice", "Quot.sound"}
for theorem, raw in audits:
    used = {name.strip() for name in raw.split(",") if name.strip()}
    extra = used - allowed
    if extra:
        raise SystemExit(f"Unexpected axioms in {theorem}: {sorted(extra)}")
print(f"PASS: {len(audits)} theorem axiom audits; standard logical axioms only.")
PY
