#!/usr/bin/env python3
"""
One-shot publish prep: fill in every placeholder, then init git and commit.

Pure Python and standard library only, so it behaves identically on macOS,
Linux and Windows (no sed/BSD-vs-GNU differences).

Usage
-----
    python3 scripts/setup_repo.py \
        --username Victor-Alfred \
        --first Victor --last Alfred \
        --email vicalf_1@yahoo.com

Optional:
    --version 0.1.0          package version (default: keep current)
    --date 2026-07-24        release date for CITATION.cff (default: today)
    --repo-name aav-chimera  GitHub repo name if different from aav-chimera
    --no-git                 only substitute placeholders, skip git init/commit

Safe to run more than once.
"""

import argparse
import datetime as _dt
import re
import subprocess
import sys
from pathlib import Path

TEXT_SUFFIXES = {".md", ".toml", ".cff", ".yml", ".yaml", ".py", ".sh", ".txt"}
EXTRA_FILES = {"LICENSE"}
SKIP_DIRS = {".git", "__pycache__", ".pytest_cache", ".ruff_cache", "results"}
THIS_FILE = Path(__file__).name


def find_repo_root() -> Path:
    """Locate the repo root by walking up until pyproject.toml is found."""
    here = Path(__file__).resolve().parent
    for candidate in (here.parent, here, Path.cwd()):
        if (candidate / "pyproject.toml").exists():
            return candidate
    sys.exit(
        "ERROR: could not find pyproject.toml.\n"
        "Run this from inside the aav-chimera folder:\n"
        "    python3 scripts/setup_repo.py --help"
    )


def iter_text_files(root: Path):
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.relative_to(root).parts):
            continue
        if path.name == THIS_FILE:
            continue
        if path.suffix in TEXT_SUFFIXES or path.name in EXTRA_FILES:
            yield path


def substitute(root: Path, mapping: dict) -> int:
    """Apply literal string replacements across all text files."""
    changed = 0
    for path in iter_text_files(root):
        try:
            original = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        updated = original
        for old, new in mapping.items():
            updated = updated.replace(old, new)
        if updated != original:
            path.write_text(updated, encoding="utf-8")
            changed += 1
            print(f"   updated {path.relative_to(root)}")
    return changed


def set_version(root: Path, version: str) -> None:
    """Update the version string in pyproject.toml, CITATION.cff and __init__.py."""
    targets = [
        (root / "pyproject.toml", r'^version = ".*"$', f'version = "{version}"'),
        (root / "CITATION.cff", r"^version: .*$", f"version: {version}"),
        (root / "src" / "aav_chimera" / "__init__.py",
         r'^__version__ = ".*"$', f'__version__ = "{version}"'),
    ]
    for path, pattern, replacement in targets:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        new_text, n = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
        if n:
            path.write_text(new_text, encoding="utf-8")
            print(f"   version -> {version} in {path.relative_to(root)}")


def set_dates(root: Path, release_date: str) -> None:
    """Set CITATION.cff release date and the LICENSE copyright year."""
    cff = root / "CITATION.cff"
    if cff.exists():
        text = cff.read_text(encoding="utf-8")
        text, n = re.subn(r"^date-released: .*$",
                          f"date-released: {release_date}", text,
                          count=1, flags=re.MULTILINE)
        if n:
            cff.write_text(text, encoding="utf-8")
            print(f"   release date -> {release_date}")

    lic = root / "LICENSE"
    if lic.exists():
        year = release_date.split("-")[0]
        text = lic.read_text(encoding="utf-8")
        text, n = re.subn(r"Copyright \(c\) \d{4}",
                          f"Copyright (c) {year}", text, count=1)
        if n:
            lic.write_text(text, encoding="utf-8")
            print(f"   copyright year -> {year}")


def strip_publish_note(root: Path) -> None:
    """Remove the pre-publication reminder block from the README."""
    readme = root / "README.md"
    if not readme.exists():
        return
    text = readme.read_text(encoding="utf-8")
    cleaned = re.sub(
        r"\n*---\n*\n*<sub><b>Before publishing:.*?</sub>\s*$",
        "\n",
        text,
        flags=re.DOTALL,
    )
    if cleaned != text:
        readme.write_text(cleaned.rstrip() + "\n", encoding="utf-8")
        print("   removed pre-publication note from README.md")


def check_clean(root: Path) -> bool:
    """Report any placeholders that survived."""
    tokens = ["yourusername", "YOUR NAME", "YOUR FIRST NAME",
              "YOUR LAST NAME", "you@example.com"]
    found = []
    for path in iter_text_files(root):
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for token in tokens:
            if token in text:
                found.append(f"{path.relative_to(root)}: {token}")
    if found:
        print("\n   WARNING — placeholders still present:")
        for item in found:
            print(f"     {item}")
        return False
    print("   no placeholders remain.")
    return True


def run(cmd, cwd, check=True):
    return subprocess.run(cmd, cwd=cwd, check=check,
                          capture_output=True, text=True)


def git_setup(root: Path, name: str, email: str) -> bool:
    """Initialise the repo on branch main and create the initial commit."""
    if not (root / ".git").exists():
        print(">> git init (branch: main)")
        try:
            run(["git", "init", "-b", "main"], root)
        except subprocess.CalledProcessError:
            # Older git without -b support
            run(["git", "init"], root)
            run(["git", "checkout", "-b", "main"], root, check=False)
    else:
        print(">> existing .git found — reusing it")

    # Ensure the branch is named main regardless of the user's git defaults.
    branch = run(["git", "rev-parse", "--abbrev-ref", "HEAD"], root,
                 check=False).stdout.strip()
    if branch and branch != "main":
        run(["git", "branch", "-M", "main"], root, check=False)
        print(f"   renamed branch {branch} -> main")

    run(["git", "config", "user.name", name], root)
    run(["git", "config", "user.email", email], root)
    run(["git", "add", "-A"], root)

    staged = run(["git", "diff", "--cached", "--name-only"], root,
                 check=False).stdout.strip()
    if not staged:
        print("   nothing new to commit (already committed?)")
        return True

    message = (
        "Initial commit: benchmark-validated AAV chimeric read detection\n\n"
        "Detection pipeline, ground-truth Nanopore read simulator, and the\n"
        "benchmark harness that scores one against the other.\n\n"
        "Chimeric detection: F1 0.973, precision 1.000, recall 0.947 on 5,071\n"
        "reads with known ground truth (0 false positives).\n\n"
        "Includes unit tests, ruff linting, CI across Python 3.9/3.11/3.12,\n"
        "and a benchmark workflow that regenerates figures on every push."
    )
    run(["git", "commit", "-m", message], root)
    count = len(staged.splitlines())
    print(f"   committed {count} files on branch main")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Fill placeholders and prepare the repo for its first push.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--username", required=True, help="GitHub username")
    parser.add_argument("--first", required=True, help="Your first name")
    parser.add_argument("--last", required=True, help="Your last name")
    parser.add_argument("--email", required=True, help="Contact email")
    parser.add_argument("--version", default=None,
                        help="Package version (default: leave unchanged)")
    parser.add_argument("--date", default=None,
                        help="Release date YYYY-MM-DD (default: today)")
    parser.add_argument("--repo-name", default="aav-chimera",
                        help="GitHub repo name (default: aav-chimera)")
    parser.add_argument("--no-git", action="store_true",
                        help="Only substitute placeholders; skip git init/commit")
    args = parser.parse_args()

    release_date = args.date or _dt.date.today().isoformat()
    try:
        _dt.date.fromisoformat(release_date)
    except ValueError:
        sys.exit(f"ERROR: --date must be YYYY-MM-DD, got '{release_date}'")

    root = find_repo_root()
    full_name = f"{args.first} {args.last}"
    print(f">> Repo root: {root}")
    print(">> Substituting placeholders")

    mapping = {
        "YOUR FIRST NAME": args.first,
        "YOUR LAST NAME": args.last,
        "YOUR NAME": full_name,
        "you@example.com": args.email,
        "yourusername": args.username,
    }
    if args.repo_name != "aav-chimera":
        mapping[f"{args.username}/aav-chimera"] = f"{args.username}/{args.repo_name}"

    substitute(root, mapping)

    if args.version:
        set_version(root, args.version)
    set_dates(root, release_date)
    strip_publish_note(root)

    print(">> Verifying")
    check_clean(root)

    if args.no_git:
        print("\nSkipped git (--no-git).")
        return

    git_setup(root, full_name, args.email)

    print(f"""
Done. Now create the GitHub repo and push:

  gh repo create {args.repo_name} --public --source=. --remote=origin --push

If that fails, create an EMPTY repo named '{args.repo_name}' at github.com/new
(no README, no license, no .gitignore) and then:

  git remote add origin https://github.com/{args.username}/{args.repo_name}.git
  git push -u origin main
""")


if __name__ == "__main__":
    main()
