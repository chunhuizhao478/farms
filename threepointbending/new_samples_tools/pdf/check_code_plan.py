#!/usr/bin/env python3
"""
check_code_plan.py -- the /code-plan skill's mechanical rules, checked on a plan-format Markdown file.

WHAT
    ASCII source; inline `$` balanced on every line; every `(#anchor)` link resolves to a heading id; table cells of at
    most CELL_WORDS words; the prose word budgets of the skill per section; at most 6 constraints, 5 risks, 6 phases.

USAGE
    python3 check_code_plan.py DOC.md          exit status = number of rule failures
"""
import re
import sys
from pathlib import Path

CELL_WORDS = 15
BUDGETS = {"Summary": 150, "Formulation": 120, "Approach": 100, "Testing Strategy": 80}
PHASE_BUDGET = 200


def auto_id(title: str) -> str:
    """pandoc's auto_identifiers for a heading without an explicit {#id}."""
    t = re.sub(r"\\eqref\{[^}]*\}|`", "", title).lower()
    t = re.sub(r"[^a-z0-9 _.-]", "", t)
    return re.sub(r"\s+", "-", t.strip())


def strip_code_spans(line: str) -> str:
    return re.sub(r"`[^`]*`", "", line)


def prose_words(lines):
    """Words of narration: no tables, display math, code, checkboxes, numbered or bulleted spec lists, provenance tags."""
    n, in_eq, in_code = 0, False, False
    for line in lines:
        s = line.strip()
        if s.startswith("```"):
            in_code = not in_code
            continue
        if in_code:
            continue
        if s.startswith("\\begin{equation}"):
            in_eq = True
        if in_eq:
            in_eq = not s.startswith("\\end{equation}")
            continue
        if (not s or s.startswith("|") or s.startswith("- [") or s.startswith("$$") or re.match(r"^(\d+\.|-)\s", s)
                or s.startswith("*From") or s.startswith("*Derived") or s.startswith("*Definition")
                or s.startswith("*Modeling") or re.match(r"^\s{2,}", line)):
            continue
        n += len(strip_code_spans(s).split())
    return n


def main(path: str) -> int:
    text = Path(path).read_text()
    lines = text.split("\n")
    fails = 0

    bad = [(i + 1, c) for i, l in enumerate(lines) for c in l if ord(c) > 127]
    print(f"non-ASCII characters: {len(bad)}" + (f"  first at line {bad[0][0]}: {bad[0][1]!r}" if bad else ""))
    fails += bool(bad)

    in_code, dollar_bad = False, []
    for i, l in enumerate(lines):
        if l.strip().startswith("```"):
            in_code = not in_code
            continue
        if not in_code and strip_code_spans(l).replace("\\$", "").count("$") % 2:
            dollar_bad.append(i + 1)
    print(f"lines with an odd number of $: {len(dollar_bad)}" + (f"  {dollar_bad[:8]}" if dollar_bad else ""))
    fails += bool(dollar_bad)

    ids = set()
    for l in lines:
        m = re.match(r"^#{1,6}\s+(.*?)\s*(\{#([^}]+)\})?\s*$", l)
        if m:
            ids.add(m.group(3) if m.group(3) else auto_id(m.group(1)))
    links = set(re.findall(r"\]\(#([^)]+)\)", text))
    missing = sorted(links - ids)
    print(f"links: {len(links)} distinct anchors, unresolved: {missing}")
    fails += bool(missing)

    long_cells = []
    for i, l in enumerate(lines):
        if l.startswith("|") and not re.match(r"^\|[\s:|-]+\|$", l):
            for cell in l.strip().strip("|").split("|"):
                if len(cell.split()) > CELL_WORDS:
                    long_cells.append((i + 1, len(cell.split()), cell.strip()[:60]))
    print(f"table cells over {CELL_WORDS} words: {len(long_cells)}")
    for row in long_cells[:10]:
        print("   line %d: %d words: %s" % row)
    fails += bool(long_cells)

    sections, cur = {}, None
    for l in lines:
        m = re.match(r"^##\s+(.*?)\s*(\{#[^}]+\})?\s*$", l)
        if m:
            cur = m.group(1)
            sections[cur] = []
        elif cur:
            sections[cur].append(l)
    print("prose words per section (budget):")
    phases = 0
    for name, body in sections.items():
        words = prose_words(body)
        budget = BUDGETS.get(name, PHASE_BUDGET if name.startswith("Phase ") and ":" in name else None)
        phases += name.startswith("Phase ") and ":" in name
        over = budget is not None and words > budget
        fails += over
        if budget is not None or words:
            print(f"   {name[:48]:<48} {words:>5}" + (f" / {budget}" if budget else "") + ("   OVER" if over else ""))
    n_constraints = sum(1 for l in sections.get("Constraints", []) if l.startswith("- "))
    n_risks = sum(1 for l in sections.get("Risks", []) if l.startswith("|")) - 2
    print(f"constraints {n_constraints} (<= 6), risks {n_risks} (<= 5), phases {phases} (<= 6)")
    fails += (n_constraints > 6) + (n_risks > 5) + (phases > 6)
    print(f"FAILURES: {fails}")
    return fails


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
