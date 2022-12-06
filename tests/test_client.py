from click.testing import CliRunner

from salmon_rnaseq import cli


def test_command_line_interface(fixture_setup):
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    # should be able to handle no args and print help
    assert result.exit_code == 0
    assert "--help" in result.output
    # direct help message call
    help_result = runner.invoke(cli.main, ["--help"])
    assert help_result.exit_code == 0
    assert "--help" in help_result.output


def test_command_line_10x_v2_sn(fixture_setup):
    runner = CliRunner()
    dir_10x_v2_sn = fixture_setup.get_dir_10x_v2_sn()
    proc = runner.invoke(
        cli.main,
        [
            "--assay",
            "10x_v2_sn",
            "--fastq_dir",
            dir_10x_v2_sn + "/input",
            "-o",
            dir_10x_v2_sn + "/output",
        ],
    )
    print(proc)
