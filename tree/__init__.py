"""Tree data structures and the Newick parser."""

from .tree_parser import (
    ParseError, TNode, TTree,
    allocate_sequences, dispose_tree, parse_trees)

__all__ = [
    "TNode", "TTree", "ParseError",
    "parse_trees",
    "allocate_sequences", "dispose_tree",
]
