"""starkit.__main__: executed when starkit is called as script"""
import sys

from .starkit import main

main(sys.argv[1:])
