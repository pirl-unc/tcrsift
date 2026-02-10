# til-select implementation plan

## Goal
Add `tcrsift til-select` to support TIL-only prioritization from 10x VDJ + GEX timepoints, with CLI inputs/outputs aligned to:

`~/code/pfo-analysis/pfo004/full-length-tcrs-in-TILs/v2/harmonize_abtcr_timepoints.py`

## Scope and constraints
- Keep existing `til-clonotype` and `match-til` behavior unchanged.
- Reuse/factor existing annotation and utility functionality where possible.
- Add tests first (TDD) for new functionality.
- Keep CLI backward-compatible with existing commands.

## Task list

1. Baseline contracts
- [ ] Capture CLI contract for v2-compatible inputs and default outputs.
- [ ] Define minimal required output artifact list to match v2 workflow.

2. Tests first (fail initially)
- [ ] Add parser tests for new `til-select` command and key arguments.
- [ ] Add unit tests for sample/config parsing and input discovery.
- [ ] Add unit tests for consensus/clonotype loading and harmonization.
- [ ] Add unit tests for marker scoring from contigs + 10x h5.
- [ ] Add unit tests for selection masks and subset generation.
- [ ] Add CLI integration-style test for end-to-end small synthetic dataset.

3. Implementation
- [ ] Add new module `tcrsift/til_select.py` with:
  - input parsing/discovery
  - harmonization
  - marker scoring
  - selection pipeline
  - output writing and basic plots/report
- [ ] Reuse existing annotation loading/matching (`tcrsift.annotate`) for DB matching.
- [ ] Reuse validation helpers for path and numeric checks where appropriate.
- [ ] Keep output filenames aligned with v2 script defaults/names.

4. CLI wiring
- [ ] Add `cmd_til_select` in `tcrsift/cli.py`.
- [ ] Add `til-select` parser section with v2-compatible flags and aliases.
- [ ] Export public helpers from `tcrsift/__init__.py` as appropriate.

5. Documentation
- [ ] Update `README.md` with TIL-only selection workflow examples.
- [ ] Update `docs/user-guide/cli.md` with full `til-select` command docs.
- [ ] Update `docs/user-guide/pipeline.md` for TIL-only analysis guidance.
- [ ] Add/extend API docs for `tcrsift.til_select`.

6. Validation
- [ ] Run targeted tests for new module and CLI.
- [ ] Run lint (`bash lint.sh`).
- [ ] Run full tests (`bash test.sh`).
- [ ] Summarize parity and any known deltas vs v2 behavior.
