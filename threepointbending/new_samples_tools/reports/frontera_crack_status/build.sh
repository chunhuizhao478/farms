#!/bin/bash
# build.sh -- REPORT.md -> REPORT.pdf (house /code-plan document build, same tools as ../../pdf/build.sh).
# pandoc + xelatex, TOC to depth 2, RoyalBlue links, ../../pdf/code_plan_header.tex (fvextra wrapping, footnotesize
# tables) and the ../../pdf/codebreak.lua filter; ../../pdf/check_code_plan.py runs first (ASCII, $ balance, links,
# cell and word budgets).  Figures come from render_figures.py (run it first).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PDF_TOOLS="$(cd "$HERE/../../pdf" && pwd)"
DOC=REPORT
export PATH=/opt/homebrew/bin:/Library/TeX/texbin:$PATH
python3 "$PDF_TOOLS/check_code_plan.py" "$HERE/$DOC.md"
cd "$HERE"
pandoc $DOC.md -s -o $DOC.tex --pdf-engine=xelatex --toc --toc-depth=2 \
  -V colorlinks=true -V linkcolor=RoyalBlue -V toccolor=black -H "$PDF_TOOLS/code_plan_header.tex" \
  -V monofont="Menlo" -V monofontoptions="Scale=0.84" \
  --lua-filter="$PDF_TOOLS/codebreak.lua"
for pass in 1 2 3; do xelatex -interaction=nonstopmode $DOC.tex > /dev/null 2>&1 || true; done
echo "$HERE/$DOC.pdf: $(grep -a -o 'Output written on [^ ]* ([0-9]* pages)' $DOC.log | sed 's/.*(//; s/)//')"
echo "errors: $(grep -c '^!' $DOC.log || true)   missing chars: $(grep -c 'Missing character' $DOC.log || true)"
echo "undefined references: $(grep -c 'undefined' $DOC.log || true)"
echo "overfull > 15 pt:"; grep -o 'Overfull \\hbox ([0-9.]*pt' $DOC.log | sed 's/[^0-9.]//g' | awk '$1>15' || true
[ -n "${KEEP_LOG:-}" ] || rm -f $DOC.aux $DOC.out $DOC.toc $DOC.tex $DOC.log
