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
    if evalue <= 0 or evalue > 1:
        logger.warning("E-value must be between 0 (exclusive) and 1. Using default.")
        evalue = DEFAULT_EVALUE

    min_size = int(args.min_size) if args.min_size is not None else DEFAULT_MIN_SIZE
    max_size = int(args.max_size) if args.max_size is not None else DEFAULT_MAX_SIZE
    if min_size < 0:
        min_size = 0
    if max_size < 0:
        max_size = 0
    evidence = args.evidence or "all"
    use_log = args.log or False
    quiet = args.quiet or False
    no_homology = getattr(args, "no_homology", False) or False
    library = getattr(args, "library", None)
    if library and not os.path.isfile(library):
        logger.warning(f"Library file does not exist: {library}")
        library = None
    run_mode = getattr(args, "mode", "relaxed") or "relaxed"

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
        no_homology=no_homology,
        library=library,
        run_mode=run_mode,
    )
