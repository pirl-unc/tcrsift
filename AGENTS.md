# Working agreements for tcrsift

## Release workflow (hard rule)

**Merging a PR to `main` and releasing it are ONE coupled action. Never merge and then wonder whether it's ready to release.**

- Do *all* work — including the version bump in `tcrsift/version.py` — on the feature branch / open PR. Keep iterating there (CI green) until it is genuinely release-ready.
- Only then: merge (squash, with a `X.Y.Z: … (#PR)` title), **and in the same pass** tag `vX.Y.Z` and publish to PyPI (`python -m build` + `twine upload dist/tcrsift-X.Y.Z*`).
- Do **not** merge to "make progress", to bank CI-green work, or while still deciding whether to release. If it's not ready to ship, it stays on the branch.
- Rationale: the pilot installs from PyPI. `main` running ahead of the last published version creates a catch-up gap that grows and gets riskier.

## Release mechanics

- Per-PR version bump; squash-merge title is `X.Y.Z: <summary> (#PR)`.
- Tag `vX.Y.Z` after merge; publish wheel + sdist to PyPI (`~/.pypirc` is configured).
- CI (tests 3.9–3.12, lint, build) must be green. Coveralls is informational (no branch protection).
- `gh pr edit`/`gh issue view` hit a deprecated projects-classic GraphQL bug — use `gh api` (REST) for PR/issue edits.
