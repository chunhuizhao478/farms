#!/bin/bash
# build.sh -- new_samples_tools/README.md -> README.pdf (house /code-plan document build).
# pandoc + xelatex, TOC to depth 2, RoyalBlue links, code_plan_header.tex (fvextra wrapping,
# footnotesize tables) and the codebreak.lua filter, all copied into this folder so the build is
# self-contained.  check_code_plan.py runs first (ASCII, $ balance, links, cell and word budgets).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DOC_DIR="$(dirname "$HERE")"
DOC=README
export PATH=/opt/homebrew/bin:/Library/TeX/texbin:$PATH
python3 "$HERE/check_code_plan.py" "$DOC_DIR/$DOC.md"
cd "$DOC_DIR"
pandoc $DOC.md -s -o $DOC.tex --pdf-engine=xelatex --toc --toc-depth=2 \
  -V colorlinks=true -V linkcolor=RoyalBlue -V toccolor=black -H "$HERE/code_plan_header.tex" \
  -V monofont="Menlo" -V monofontoptions="Scale=0.84" \
  --lua-filter="$HERE/codebreak.lua"
for pass in 1 2 3; do xelatex -interaction=nonstopmode $DOC.tex > /dev/null 2>&1 || true; done
echo "$DOC_DIR/$DOC.pdf: $(grep -a -o 'Output written on [^ ]* ([0-9]* pages)' $DOC.log | sed 's/.*(//; s/)//')"
echo "errors: $(grep -c '^!' $DOC.log || true)   missing chars: $(grep -c 'Missing character' $DOC.log || true)"
echo "undefined references: $(grep -c 'undefined' $DOC.log || true)"
echo "overfull > 15 pt:"; grep -o 'Overfull \\hbox ([0-9.]*pt' $DOC.log | sed 's/[^0-9.]//g' | awk '$1>15' || true
[ -n "${KEEP_LOG:-}" ] || rm -f $DOC.aux $DOC.out $DOC.toc $DOC.tex $DOC.log
