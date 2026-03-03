# Publishing rnaseqmut-py to PyPI

## 1) One-time setup

1. Create a PyPI account: https://pypi.org/account/register/
2. Create project ownership for `rnaseqmut-py` (first upload creates project).
3. In GitHub repo settings:
   - Settings → Secrets and variables → Actions
   - Add repository secret `PYPI_API_TOKEN` (if you use token-based upload)
4. Optional but recommended: configure PyPI Trusted Publisher for this GitHub repo/workflow.

## 2) Build locally

```bash
cd python
python3 -m pip install --upgrade build twine
python3 -m build
python3 -m twine check dist/*
```

## 3) Upload manually (token)

```bash
python3 -m twine upload dist/*
```

Use `__token__` as username and your PyPI token as password.

## 4) Upload via GitHub Actions (recommended)

- Workflow: `.github/workflows/pypi-publish.yml`
- Trigger by publishing a GitHub Release.
- The workflow builds from `python/` and uploads to PyPI.

## 5) Version bump before next release

Edit `python/pyproject.toml`:

```toml
version = "0.1.1"
```

Then tag and release again.
