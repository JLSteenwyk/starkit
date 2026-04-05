import sys
import textwrap

from argparse import (
    ArgumentParser,
    SUPPRESS,
    RawDescriptionHelpFormatter,
)

from .version import __version__


def create_parser() -> ArgumentParser:
    parser = ArgumentParser(
        add_help=False,
        formatter_class=RawDescriptionHelpFormatter,
        usage=SUPPRESS,
        description=textwrap.dedent(
            f"""\
           _____ _             _  _______ _______
          / ____| |           | |/ /_   _|__   __|
         | (___ | |_ __ _ _ __| ' /  | |    | |
          \\___ \\| __/ _` | '__|  <   | |    | |
          ____) | || (_| | |  | . \\ _| |_   | |
         |_____/ \\__\\__,_|_|  |_|\\_\\_____|  |_|


        Version: {__version__}
        StarKIT: Starship prediction in fungal genomes.

        Usage: starkit <input> [optional arguments]
        """  # noqa
        ),
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()

    # required arguments
    required = parser.add_argument_group(
        "required arguments",
        description=textwrap.dedent(
            """\
        <input>                                     input genome file
                                                    (GenBank or FASTA; must be the first argument)
        """
        ),
    )

    required.add_argument("input", type=str, help=SUPPRESS)

    # optional arguments
    optional = parser.add_argument_group(
        "optional arguments",
        description=textwrap.dedent(
            """\
        -g, --gff <gff3_file>                      GFF3 annotation file
                                                    (required if input is FASTA)

        -o, --output <output_prefix>                output file prefix
                                                    (default: input file with '.starkit' suffix)

        -e, --evalue <threshold>                    e-value threshold for captain HMM search
                                                    (default: 1e-10)

        --min-size <bp>                             minimum predicted Starship size in bp
                                                    (default: 15000)

        --max-size <bp>                             maximum predicted Starship size in bp
                                                    (default: 700000)

        --evidence <level>                          minimum evidence level to report
                                                    high, medium, low, all (default: all)

        -l, --log                                   creates a log file

        -q, --quiet                                 suppress stdout messages

        -h, --help                                  help message
        -v, --version                               print version
        """
        ),
    )

    optional.add_argument(
        "-g",
        "--gff",
        type=str,
        required=False,
        help=SUPPRESS,
        metavar="GFF3 file",
    )

    optional.add_argument("-o", "--output", help=SUPPRESS, metavar="output")

    optional.add_argument(
        "-e",
        "--evalue",
        type=float,
        required=False,
        help=SUPPRESS,
        metavar="e-value threshold",
    )

    optional.add_argument(
        "--min-size",
        type=int,
        required=False,
        help=SUPPRESS,
        metavar="minimum size",
    )

    optional.add_argument(
        "--max-size",
        type=int,
        required=False,
        help=SUPPRESS,
        metavar="maximum size",
    )

    evidence_choices = ["high", "medium", "low", "all"]
    optional.add_argument(
        "--evidence",
        type=str,
        required=False,
        choices=evidence_choices,
        help=SUPPRESS,
        metavar="evidence level",
    )

    optional.add_argument(
        "-l",
        "--log",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        required=False,
        help=SUPPRESS,
    )

    optional.add_argument(
        "-h",
        "--help",
        action="help",
        help=SUPPRESS,
    )

    optional.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"starkit {__version__}",
        help=SUPPRESS,
    )

    return parser
