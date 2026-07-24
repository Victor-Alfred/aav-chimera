#!/usr/bin/env bash
#
# Rename the project (repo name, Python package, CLI entry points) in one pass.
#
# Usage:
#   ./scripts/rename_project.sh <new-repo-name> <new_package_name> <new-cli-name>
#
# Example — switching to "chimerabench":
#   ./scripts/rename_project.sh chimerabench chimerabench chimerabench
#
# Rewrites pyproject.toml, README.md, docs/, tests/, examples/, CITATION.cff,
# CONTRIBUTING.md, workflows, and the src/ package directory.
#
set -euo pipefail

# macOS ships BSD sed (needs -i ''), Linux ships GNU sed (needs -i).
if sed --version >/dev/null 2>&1; then SEDI=(-i); else SEDI=(-i ""); fi

NEW_REPO="${1:?Usage: rename_project.sh <new-repo-name> <new_package_name> <new-cli-name>}"
NEW_PKG="${2:?Missing package name (must be a valid Python identifier)}"
NEW_CLI="${3:?Missing CLI name}"

OLD_REPO="aav-chimera"
OLD_PKG="aav_chimera"
OLD_CLI="aav-chimera"
OLD_CLI_SIM="aav-chimera-sim"
NEW_CLI_SIM="${NEW_CLI}-sim"

if [[ ! "$NEW_PKG" =~ ^[a-z_][a-z0-9_]*$ ]]; then
  echo "ERROR: package name '$NEW_PKG' is not a valid Python identifier." >&2
  exit 1
fi

echo ">> Renaming package directory src/$OLD_PKG -> src/$NEW_PKG"
if [ -d "src/$OLD_PKG" ]; then
  git mv "src/$OLD_PKG" "src/$NEW_PKG" 2>/dev/null || mv "src/$OLD_PKG" "src/$NEW_PKG"
fi

echo ">> Rewriting references"
FILES=$(grep -rl "$OLD_REPO\|$OLD_PKG\|$OLD_CLI" \
  --include="*.py" --include="*.md" --include="*.toml" \
  --include="*.yml" --include="*.cff" --include="*.sh" --include="*.txt" \
  . 2>/dev/null || true)

for f in $FILES; do
  # Order matters: longest / most specific patterns first.
  sed "${SEDI[@]}" "s/${OLD_CLI_SIM}/${NEW_CLI_SIM}/g" "$f"
  sed "${SEDI[@]}" "s/${OLD_REPO}/${NEW_REPO}/g" "$f"
  sed "${SEDI[@]}" "s/${OLD_PKG}/${NEW_PKG}/g" "$f"
  sed "${SEDI[@]}" "s/${OLD_CLI}/${NEW_CLI}/g" "$f"
done

echo ">> Remaining references to the old name (should be empty):"
grep -rn "$OLD_REPO\|$OLD_PKG" \
  --include="*.py" --include="*.md" --include="*.toml" \
  --include="*.yml" --include="*.cff" --include="*.sh" . 2>/dev/null \
  || echo "   none — clean."

cat <<MSG

Done. Next steps:
  1. Rename the containing directory:  cd .. && mv $OLD_REPO $NEW_REPO
  2. Reinstall and verify:             pip install -e ".[dev]" && ruff check src tests scripts && pytest
  3. Re-read the README title/tagline by hand — prose isn't safely automatable.
MSG
