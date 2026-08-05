#!/usr/bin/env python3
"""Generate + run a standalone equivalence test for computeAlphaCr.

Extracts the polynomial expressions from three sources -- pre-change main, this
branch's patched file, and cdbm@99b41e9c -- and compiles them into one program so
they can be compared without building MOOSE.

Asserts:
  soft path  == pre-change main      EXACTLY (bit-for-bit) -> protects existing users
  stiff path == cdbm@99b41e9c        within 1 ulp          -> gives the paper its physics
"""
import pathlib, re, subprocess, sys, tempfile

WT = pathlib.Path("/Users/chunhuizhao/projects/farms_benchmark_poroelastic/farms/.claude/worktrees/damage-paper")
REL = "src/kernels/cdbm/BreakageVarForcingFuncDevOld.C"

def git_show(ref):
    return subprocess.run(["git", "show", f"{ref}:{REL}"], cwd=WT,
                          capture_output=True, text=True, check=True).stdout

def live_exprs(text):
    """Return the two live 'alphacr = ...;' right-hand sides, in file order."""
    out = []
    for l in text.split("\n"):
        s = l.strip()
        if s.startswith("alphacr = ") and not s.startswith("//") and "1.0;" not in s:
            rhs = s[len("alphacr = "):].rstrip(";")
            rhs = re.sub(r"^alphacr = ", "", rhs)   # published doubled assignment
            out.append(rhs)
    return out

old_main = live_exprs(git_show("origin/main"))
published = live_exprs(git_show("99b41e9c"))
patched_txt = (WT / REL).read_text()

# patched file has 4 live expressions: stiff1, soft1, stiff2, soft2
patched = live_exprs(patched_txt)
if len(old_main) != 2 or len(published) != 2 or len(patched) != 4:
    sys.exit(f"unexpected expression counts: main={len(old_main)} pub={len(published)} patched={len(patched)}")
stiff1, soft1, stiff2, soft2 = patched

cpp = f"""
#include <cstdio>
#include <cmath>
#include <cstdlib>
using std::pow; using std::sqrt; using std::fabs;

static double old_main_b1(double xi){{ return {old_main[0]}; }}
static double old_main_b2(double xi){{ return {old_main[1]}; }}
static double published_b1(double xi){{ return {published[0]}; }}
static double published_b2(double xi){{ return {published[1]}; }}
static double soft_b1(double xi){{ return {soft1}; }}
static double soft_b2(double xi){{ return {soft2}; }}
static double stiff_b1(double xi){{ return {stiff1}; }}
static double stiff_b2(double xi){{ return {stiff2}; }}

static int ulp_diff(double a, double b) {{
  if (a == b) return 0;
  double m = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
  double eps = std::nextafter(m, m*2.0) - m;
  return (int)(fabs(a-b)/eps + 0.5);
}}

int main() {{
  const double xi_0 = -0.8, xi_1 = 0.8248, xi_max = 1.5;
  const double step = 1e-4;
  long n1 = 0, n2 = 0;
  long soft_exact_fail = 0, stiff_ulp_fail = 0; int worst_ulp = 0;

  for (double xi = xi_0 + step; xi <= xi_1; xi += step) {{
    ++n1;
    if (soft_b1(xi) != old_main_b1(xi)) ++soft_exact_fail;
    int u = ulp_diff(stiff_b1(xi), published_b1(xi));
    if (u > worst_ulp) worst_ulp = u;
    if (u > 1) ++stiff_ulp_fail;
  }}
  for (double xi = xi_1 + step; xi <= xi_max; xi += step) {{
    ++n2;
    if (soft_b2(xi) != old_main_b2(xi)) ++soft_exact_fail;
    int u = ulp_diff(stiff_b2(xi), published_b2(xi));
    if (u > worst_ulp) worst_ulp = u;
    if (u > 1) ++stiff_ulp_fail;
  }}

  printf("samples: branch1=%ld branch2=%ld (step %g)\\n", n1, n2, step);
  printf("soft  vs pre-change main : %s (%ld bit-level mismatches)\\n",
         soft_exact_fail ? "FAIL" : "PASS exact", soft_exact_fail);
  printf("stiff vs cdbm@99b41e9c   : %s (%ld beyond 1 ulp, worst %d ulp)\\n",
         stiff_ulp_fail ? "FAIL" : "PASS", stiff_ulp_fail, worst_ulp);

  // Sanity: the two calibrations must actually differ, else the selector is a no-op.
  double d = fabs(stiff_b1(0.0) - soft_b1(0.0));
  printf("calibrations differ at xi=0: stiff=%.9g soft=%.9g |diff|=%.9g\\n",
         stiff_b1(0.0), soft_b1(0.0), d);
  if (d < 1e-12) {{ printf("FAIL: selector is a no-op\\n"); return 1; }}

  return (soft_exact_fail || stiff_ulp_fail) ? 1 : 0;
}}
"""

tmp = pathlib.Path(tempfile.mkdtemp())
srcf = tmp / "t.cpp"
srcf.write_text(cpp)
exe = tmp / "t"
r = subprocess.run(["g++", "-O0", "-std=c++17", str(srcf), "-o", str(exe)],
                   capture_output=True, text=True)
if r.returncode:
    print(r.stderr[:3000]); sys.exit("compile failed")
print(subprocess.run([str(exe)], capture_output=True, text=True).stdout)
sys.exit(subprocess.run([str(exe)]).returncode)
