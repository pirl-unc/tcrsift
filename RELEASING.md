# Releasing TCRsift

Merging a release PR to `main` and publishing it to PyPI are performed together
in one pass so `main` never runs ahead of the latest published package:

1. Complete all changes — including the version bump in `tcrsift/version.py`
   and changelog entry in `docs/release-notes.md` — on the PR branch.
2. Wait until all CI checks (`test` 3.9–3.12, `lint`, `build`, docs `build`)
   are green.
3. Squash-merge the PR with commit title `X.Y.Z: <summary> (#PR)`.
4. Immediately tag `vX.Y.Z` (GitHub release) and upload the built distributions:

```bash
rm -rf dist/ build/
python -m build
twine check dist/tcrsift-X.Y.Z*
twine upload dist/tcrsift-X.Y.Z*
```
