# GitHub-first guard: the checkout must be clean and its HEAD must be on GitHub. Sourced by every LS6 job script.
cd $FARMS_DIR || exit 2
git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "ERROR: uncommitted changes in $FARMS_DIR"; git status --short; exit 2; }
[ -n "$(git branch -r --contains HEAD)" ] || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
echo "farms commit $(git rev-parse HEAD)"
