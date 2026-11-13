"""Alignment and rate output writers: strict PHYLIP, relaxed PHYLIP, NEXUS, and FASTA.

Maps integer state buffers to ASCII sequences via a cached byte translation
table and formats tip or ancestral alignments, plus the per-site relative rate
table for ``-wr``.
"""

import numpy as np

_TABLES = {}


def _state_table(state_chars):
    """Byte translation table cached per alphabet (states are 0..n-1)."""
    table = _TABLES.get(state_chars)
    if table is None:
        head = state_chars.encode("ascii")
        table = head + bytes(range(len(head), 256))
        _TABLES[state_chars] = table
    return table


def _encode(seq, num_sites, state_chars):
    """Characters for the first *num_sites* states of *seq*."""
    return seq[:num_sites].astype(np.uint8).tobytes().translate(
        _state_table(state_chars)).decode("ascii")


def _datatype(state_chars):
    return "DNA" if len(state_chars) == 4 else "PROTEIN"


def write_sequences(f, fmt, tree_set, partition_lengths, num_taxa,
                    total_sites, state_chars, tree_no=0, dataset_no=0,
                    is_nuc=True):
    """Write the tip sequences of every partition in *fmt*."""
    if fmt == "n":
        title = ""
        if tree_no and dataset_no:
            title = f"Tree {tree_no}, Dataset {dataset_no}"
        elif tree_no:
            title = f"Tree {tree_no}"
        elif dataset_no:
            title = f"Dataset {dataset_no}"
        rows = [(_name(tree_set, i),
                 "".join(_encode(tree.tips[i].sequence, partition_lengths[k],
                                 state_chars)
                         for k, tree in enumerate(tree_set)))
                for i in range(num_taxa)]
        _write_nexus_block(f, title, rows, "DNA" if is_nuc else "PROTEIN",
                           total_sites)
        return

    if fmt == "f":
        for i in range(num_taxa):
            f.write(f">{_name(tree_set, i)}\n")
            seq = "".join(_encode(tree.tips[i].sequence, partition_lengths[k],
                                  state_chars)
                          for k, tree in enumerate(tree_set))
            for j in range(0, len(seq), 72):
                f.write(seq[j:j + 72] + "\n")
        return

    relaxed = fmt == "r"
    f.write(f"{num_taxa} {total_sites}\n")
    for i in range(num_taxa):
        name = _name(tree_set, i)
        f.write(f"{name} " if relaxed else f"{name[:10]:<10}")
        for k, tree in enumerate(tree_set):
            f.write(_encode(tree.tips[i].sequence, partition_lengths[k],
                            state_chars))
        f.write("\n")


def _name(tree_set, tip_no):
    return tree_set[0].names[tip_no]


def _write_nexus_block(f, title, rows, datatype, nchar):
    """Emit a NEXUS DATA block from ``(name, sequence)`` pairs."""
    f.write(f"Begin DATA; [{title}]\n" if title else "Begin DATA;\n")
    f.write(f"  Dimensions NTAX={len(rows)} NCHAR={nchar};\n")
    f.write(f"  Format MISSING=? GAP=- DATATYPE={datatype};\n")
    f.write("  Matrix\n")
    width = max(len(name) for name, _ in rows)
    for name, seq in rows:
        f.write(name + " " * (width - len(name) + 1) + seq + "\n")
    f.write("  ;\nEnd;\n\n")


def write_ancestral(f, fmt, tree, state_chars):
    """Write the sequences of every node, tips and ancestors alike.

    Internal nodes are numbered from ``num_tips + 1`` upwards in pre-order.
    """
    rows = _collect_ancestral(tree, state_chars)
    num_sites = len(tree.root.sequence)

    if fmt == "f":
        for name, seq in rows:
            f.write(f">{name}\n")
            for i in range(0, len(seq), 72):
                f.write(seq[i:i + 72] + "\n")
    elif fmt == "n":
        _write_nexus_block(f, "Ancestral sequences", rows,
                           _datatype(state_chars), num_sites)
    else:
        f.write(f"{len(rows)} {num_sites}\n")
        relaxed = fmt == "r"
        for name, seq in rows:
            f.write(f"{name} {seq}\n" if relaxed else f"{name:<10}{seq}\n")


def _collect_ancestral(tree, state_chars):
    """``(label, sequence)`` for every node, root first, in pre-order."""
    rows = []
    next_number = [tree.num_tips]

    def visit(node):
        if node is None:
            return
        if node.tip_no == -1:
            next_number[0] += 1
            label = str(next_number[0])
        else:
            label = tree.names[node.tip_no]
        rows.append((label, _encode(node.sequence, len(node.sequence),
                                    state_chars)))
        if node.tip_no == -1:
            visit(node.branch1)
            visit(node.branch2)

    visit(tree.root)
    if not tree.rooted:
        visit(tree.root.branch0)
    return rows


def write_rates(f, site_rates, num_sites):
    """The per-site relative rate table behind ``-wr``."""
    f.write("Relative rates for each site:\n")
    for i in range(num_sites):
        f.write(f"{i + 1}\t{site_rates[i]:.6f}\n")
    f.write("\n")
