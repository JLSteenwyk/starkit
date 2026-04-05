import logging
import os.path
import sys

from .settings import DEFAULT_EVALUE, DEFAULT_MIN_SIZE, DEFAULT_MAX_SIZE

logger = logging.getLogger(__name__)


def process_args(args) -> dict:
    """Process args from argparser and set defaults."""
    input_file = args.input
    output_prefix = args.output or f"{input_file}.starkit"
    gff_file = args.gff

    if not os.path.isfile(input_file):
        logger.warning("Input file does not exist")
        sys.exit()

    if gff_file and not os.path.isfile(gff_file):
        logger.warning("GFF3 file does not exist")
        sys.exit()

    evalue = float(args.evalue) if args.evalue is not None else DEFAULT_EVALUE
    min_size = int(args.min_size) if args.min_size is not None else DEFAULT_MIN_SIZE
    max_size = int(args.max_size) if args.max_size is not None else DEFAULT_MAX_SIZE
    evidence = args.evidence or "all"
    use_log = args.log or False
    quiet = args.quiet or False

    return dict(
        input_file=input_file,
        output_prefix=output_prefix,
        gff_file=gff_file,
        evalue=evalue,
        min_size=min_size,
        max_size=max_size,
        evidence=evidence,
        use_log=use_log,
        quiet=quiet,
    )
