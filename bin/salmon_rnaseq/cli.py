from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path
from subprocess import CompletedProcess
from typing import Any

import click


class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath: str):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def pathing(path: str, new: bool = False, overwrite: bool = True) -> Path:
    """Guarantees correct expansion rules for pathing.
    :param Union[str, Path] path: path of folder or file you wish to expand.
    :param bool new: will check if distination exists if new  (will check parent path regardless).
    :return: A pathlib.Path object.
    >>> pathing('~/Desktop/folderofgoodstuffs/')
    /home/user/Desktop/folderofgoodstuffs
    """
    path = Path(path)
    # Expand shortened path
    if str(path)[0] == "~":
        path = path.expanduser()
    # Exand local path
    if str(path)[0] == ".":
        path = path.resolve()
    else:
        path = path.absolute()
    # Making sure new paths don't exist while also making sure existing paths actually exist.
    if new:
        if not path.parent.exists():
            raise ValueError(f"ERROR ::: Parent directory of {path} does not exist.")
        if path.exists() and not overwrite:
            raise ValueError(f"ERROR ::: {path} already exists!")
    else:
        if not path.exists():
            raise ValueError(f"ERROR ::: Path {path} does not exist.")
    return path


class Shell:
    """Shell wrapper for cwtool pipeline.cwl"""

    @staticmethod
    def run(args: list[str], capture_output: bool = False) -> CompletedProcess[Any]:

        cmd = [
            "cwltool",
            str(Path(__file__).parent.parent.parent / "pipeline.cwl"),
            *args,
        ]
        proc: CompletedProcess[Any] = subprocess.run(
            cmd, capture_output=capture_output, encoding="utf-8"
        )
        return proc

    def help(self) -> str:
        args = ["--help"]
        help_str = self.run(args, capture_output=True).stdout
        # TODO: should be a direct change not a replace here
        # help_str = "usage: salmon-rnaseq\n" + "\n".join(help_str.split("\n")[1:])
        return help_str


class RichHelp(click.Command):
    """Override the default help command to use rich markup with CWLTOOL message."""

    def format_help(self, ctx: click.Context, formatter: click.HelpFormatter) -> None:
        click.echo(Shell().help())


@click.command(
    cls=RichHelp,
    context_settings=dict(
        ignore_unknown_options=True, allow_extra_args=True, help_option_names=["-h", "--help"]
    ),
)
@click.pass_context
@click.option(
    "-o",
    "--outdir",
    required=False,
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
    ),
    default=".",
)
def main(ctx: click.Context, outdir: str) -> None:
    """Run CWLTOOL pipeline.cwl with rich markup help.

    Parameters
    ----------
    ctx : click.Context
        Click context object.
    """
    if not ctx.args:
        click.echo(Shell().help())
        ctx.exit()

    Path(outdir).mkdir(exist_ok=True)
    for i, arg in enumerate(ctx.args):
        if Path(arg).exists():
            ctx.args[i] = str(pathing(arg))
    outdir = str(pathing(outdir))

    with cd(outdir):
        Shell.run(ctx.args)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
