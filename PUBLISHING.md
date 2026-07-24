# Publishing this repo

Four commands. Run them from inside the `aav-chimera` folder.

## 0. Check you're in the right folder

```bash
ls
```

You should see `README.md`, `pyproject.toml`, `src`, `tests`, `scripts`.

If you instead see a single `aav-chimera` folder, the archive unzipped one level
deep — `cd aav-chimera` and check again.

## 1. Fill in your details, init git, commit

```bash
python3 scripts/setup_repo.py \
  --username Victor-Alfred \
  --first Victor \
  --last Alfred \
  --email you@email.com
```

Optional flags: `--version 1.0.0`, `--date 2026-07-24`, `--repo-name other-name`.

This is pure Python (no `sed`), so it works the same on macOS and Linux. It's safe
to run more than once, and it fixes a stray `master` branch if you already ran
`git init` by hand.

Expected output ends with `committed 26 files on branch main`.

## 2. Authenticate with GitHub (once per machine)

```bash
gh auth status || gh auth login
```

If it prompts: **GitHub.com** → **HTTPS** → **Yes** → **Login with a web browser**.
Copy the 8-character code from the terminal *before* pressing Enter — that code is
what the browser page asks for. Nothing is emailed to you.

## 3. Create the repo and push

```bash
gh repo create aav-chimera --public --source=. --remote=origin --push
```

No `gh`? Create an **empty** repo at <https://github.com/new> named `aav-chimera`
— no README, no license, no .gitignore, since those would collide — then:

```bash
git remote add origin https://github.com/Victor-Alfred/aav-chimera.git
git push -u origin main
```

## 4. Confirm

```bash
gh repo view --web
```

Check the README renders with figures, then open the **Actions** tab. CI (lint +
30 tests on Python 3.9/3.11/3.12) and Benchmark (regenerates figures, publishes
metrics) should be running.

The repo must be **public** for the dynamic F1 badge to work — shields.io reads
`results/benchmark_results.json` over raw.githubusercontent, which 404s on private
repos. Badges stay grey until the first workflow run completes (~2 minutes).

---

## Troubleshooting

**`no such file or directory: ./scripts/setup_repo.py`**
You're in the wrong folder, or the archive nested one level. Run `ls` (step 0), or
`find . -name setup_repo.py` to locate it.

**`--push enabled but no commits found`**
Step 1 hasn't run successfully yet. Run it, confirm with `git log --oneline`.

**`Permission denied`**
Zip archives don't preserve the executable bit. Either `chmod +x scripts/*.sh`, or
invoke through the interpreter: `python3 scripts/setup_repo.py ...`.

**Commits don't appear on your contribution graph**
The email must be verified on your GitHub account (Settings → Emails). Otherwise use
GitHub's noreply address:

```bash
git config user.email "Victor-Alfred@users.noreply.github.com"
git commit --amend --reset-author --no-edit
```

**Repo name other than `aav-chimera`**
Pass `--repo-name your-name` in step 1 so the badge URLs are retargeted. To change
the *package* and CLI names too, see `scripts/rename_project.sh`.
