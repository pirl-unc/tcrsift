# PYTHON_ARGCOMPLETE_OK
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Command-line interface subpackage for TCRsift.

Each workflow domain owns both its parser registration and its command
handlers in a dedicated submodule:

- :mod:`tcrsift.cli.pipeline`: ``run``, ``generate-config``
- :mod:`tcrsift.cli.load`: ``load``, ``load-sct``, ``phenotype``, ``clonotype``
- :mod:`tcrsift.cli.cells`: ``cells qc|embed|annotate``
- :mod:`tcrsift.cli.filter`: ``filter``, ``select``, ``unify``, ``cohort``
- :mod:`tcrsift.cli.til`: ``prioritize``, ``til-prioritize``, ``til-select``,
  ``match-til``, ``til-clonotype``
- :mod:`tcrsift.cli.annotate`: ``annotate``, ``annotate-gex``, ``mnemonic``,
  ``log-pgen``, ``pgen fetch|train|annotate|status``
- :mod:`tcrsift.cli.assemble`: ``assemble``, ``audit-alleles``
- :mod:`tcrsift.cli.report`: ``report bundle|cohort-figures|selected|...``
- :mod:`tcrsift.cli.data`: ``data list|clear|download``
"""

from __future__ import annotations

import argparse
import logging
import sys

try:
    from datacache.cli_help import (
        HelpArgumentParser,
        HelpFormatter,
        enable_argcomplete,
    )
except ImportError:
    HelpArgumentParser = argparse.ArgumentParser
    HelpFormatter = argparse.RawDescriptionHelpFormatter

    def enable_argcomplete(parser, **_kwargs):
        try:
            import argcomplete
        except ImportError:
            return False
        argcomplete.autocomplete(parser)
        return True

from ..version import __version__
from . import annotate, assemble, cells, data, filter, load, pipeline, report, til
from .annotate import (  # noqa: F401
    cmd_annotate,
    cmd_annotate_gex,
    cmd_log_pgen,
    cmd_mnemonic,
    cmd_pgen_annotate,
    cmd_pgen_fetch,
    cmd_pgen_status,
    cmd_pgen_train,
)
from .assemble import cmd_assemble, cmd_audit_alleles  # noqa: F401
from .cells import (  # noqa: F401
    _CELLS_QC_THRESHOLDS,
    _UNSET_CLI,
    _cells_threshold,
    cmd_cells_annotate,
    cmd_cells_embed,
    cmd_cells_qc,
)
from .common import LEADER_CHOICES, _format_size, setup_logging  # noqa: F401
from .data import cmd_data_clear, cmd_data_download, cmd_data_list  # noqa: F401
from .filter import cmd_cohort, cmd_filter, cmd_select, cmd_unify  # noqa: F401
from .load import cmd_clonotype, cmd_load, cmd_load_sct, cmd_phenotype  # noqa: F401
from .pipeline import cmd_generate_config, cmd_run  # noqa: F401
from .report import (  # noqa: F401
    cmd_report_bundle,
    cmd_report_cohort,
    cmd_report_cohort_figures,
    cmd_report_polished,
    cmd_report_selected,
    cmd_report_selected_combine,
    cmd_report_sequence,
)
from .til import (  # noqa: F401
    _parse_per_lineage_min,
    cmd_match_til,
    cmd_til_clonotype,
    cmd_til_prioritize,
    cmd_til_select,
)

_PRIMARY_WORKFLOW_COMMANDS = ("prioritize", "run", "assemble", "data")
_EXPLORE_ANNOTATE_COMMANDS = ("cells", "annotate", "report")


class _TieredHelpFormatter(HelpFormatter):
    """Group top-level ``tcrsift --help`` subcommands by workflow role."""

    def _format_action(self, action):
        if not (isinstance(action, argparse._SubParsersAction) and action.dest == "command"):
            return super()._format_action(action)

        help_by_cmd = {
            sub.dest.split(" (", 1)[0]: (sub.help or "")
            for sub in action._get_subactions()
        }
        shown: set[str] = set()
        sections = [
            ("Primary workflows", [c for c in _PRIMARY_WORKFLOW_COMMANDS if c in help_by_cmd]),
            ("Exploration & annotation", [c for c in _EXPLORE_ANNOTATE_COMMANDS if c in help_by_cmd]),
        ]
        shown.update(c for _, cmds in sections for c in cmds)
        plumbing = [c for c in help_by_cmd if c not in shown]
        if plumbing:
            sections.append(("Modular pipeline steps & legacy commands", plumbing))

        lines: list[str] = []
        for idx, (heading, cmds) in enumerate(sections):
            if idx > 0:
                lines.append("")
            lines.append(f"{heading}:")
            for cmd in cmds:
                desc = " ".join(help_by_cmd[cmd].split())
                lines.append(f"  {cmd:<18} {desc}")
        lines.append("")
        return "\n".join(lines)


def create_parser() -> argparse.ArgumentParser:
    """Create the top-level TCRsift argument parser."""
    parser = HelpArgumentParser(
        prog="tcrsift",
        formatter_class=_TieredHelpFormatter,
        description=(
            "TCRsift: prioritize TCR clonotypes from fresh tissue/blood or "
            "antigen-enrichment cultures and assemble synthesis-ready constructs."
        ),
        epilog=(
            "Run 'tcrsift <command> --help' for command-specific flags.\n"
            "Quick start:\n"
            "  Ex vivo TIL / blood:  tcrsift prioritize samples.yaml -o out/ --context solid-tumor --assemble\n"
            "  Peptide / culture:    tcrsift run --sample-sheet samples.yaml -o out/\n"
            "  Reference cache:      tcrsift data download --db vdjdb --group germline"
        ),
    )
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser._positionals.title = None

    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    load.register_load(subparsers)
    cells.register(subparsers)
    load.register_phenotype(subparsers)
    load.register_clonotype(subparsers)
    filter.register_filter(subparsers)
    filter.register_select(subparsers)
    report.register(subparsers)
    annotate.register_annotate(subparsers)
    data.register(subparsers)
    til.register(subparsers)
    assemble.register_assemble(subparsers)
    assemble.register_audit_alleles(subparsers)
    annotate.register_pgen(subparsers)
    pipeline.register_run(subparsers)
    load.register_load_sct(subparsers)
    annotate.register_annotate_gex(subparsers)
    filter.register_unify(subparsers)
    annotate.register_mnemonic(subparsers)
    filter.register_cohort(subparsers)
    pipeline.register_generate_config(subparsers)

    return parser


def main(args=None):
    """Main CLI entry point."""
    parser = create_parser()
    enable_argcomplete(parser)
    args = parser.parse_args(args)

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    try:
        rc = args.func(args)
    except Exception as e:
        logging.error(f"Error: {e}")
        if hasattr(args, "verbose") and args.verbose:
            raise
        sys.exit(1)
    if isinstance(rc, int) and not isinstance(rc, bool) and rc != 0:
        sys.exit(rc)


_LAZY_PIPELINE_EXPORTS = frozenset({"_common_contig_root", "run_pipeline"})


def __getattr__(name: str):
    if name in _LAZY_PIPELINE_EXPORTS:
        from .. import pipeline as _pipeline_mod

        value = getattr(_pipeline_mod, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__) | set(_LAZY_PIPELINE_EXPORTS))


__all__ = [
    "LEADER_CHOICES",
    "cmd_annotate",
    "cmd_annotate_gex",
    "cmd_assemble",
    "cmd_audit_alleles",
    "cmd_cells_annotate",
    "cmd_cells_embed",
    "cmd_cells_qc",
    "cmd_clonotype",
    "cmd_cohort",
    "cmd_data_clear",
    "cmd_data_download",
    "cmd_data_list",
    "cmd_filter",
    "cmd_generate_config",
    "cmd_load",
    "cmd_load_sct",
    "cmd_log_pgen",
    "cmd_match_til",
    "cmd_mnemonic",
    "cmd_pgen_annotate",
    "cmd_pgen_fetch",
    "cmd_pgen_status",
    "cmd_pgen_train",
    "cmd_phenotype",
    "cmd_report_bundle",
    "cmd_report_cohort",
    "cmd_report_cohort_figures",
    "cmd_report_polished",
    "cmd_report_selected",
    "cmd_report_selected_combine",
    "cmd_report_sequence",
    "cmd_run",
    "cmd_select",
    "cmd_til_clonotype",
    "cmd_til_prioritize",
    "cmd_til_select",
    "cmd_unify",
    "create_parser",
    "main",
    "run_pipeline",  # pylint: disable=undefined-all-variable
    "setup_logging",
]
