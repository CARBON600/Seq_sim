"""Alignment writers: strict PHYLIP, relaxed PHYLIP, NEXUS and FASTA.
Integer state buffers are mapped to ASCII through a cached byte translation
table; the per-site rate table behind -wr is written here too.
input: per-taxon state buffers, partition lengths, tip names, format code
output: text on the open handle f
"""

import re

import numpy as np

# static setting
_NEXUS_PLAIN = re.compile(r"^[A-Za-z0-9_.-]+$")
_TABLES = {}


###function definition

def _nexus_label(name):
    """Quote a label NEXUS cannot read back without quotes."""
    if _NEXUS_PLAIN.match(name):
        return name
    return "'" + name.replace("'", "''") + "'"


def _relaxed_label(name):
    """Reject a name that relaxed PHYLIP cannot delimit from the sequence."""
    if re.search(r"\s", name):
        raise ValueError(
            "relaxed PHYLIP separates the name from the sequence with a "
            "space, so the name " + repr(name) + " cannot be written "
            "unambiguously; use -o f (FASTA) or rename the taxon")
    return name


def _state_table(state_chars):
    """Byte translation table cached per alphabet (states are 0..n-1)."""
    table = _TABLES.get(state_chars)
    if table is None:
        head = state_chars.encode("ascii")
        table = head + bytes(range(len(head), 256))
        _TABLES[state_chars] = table
    return table


def _encode(seq, num_sites, state_chars):
    """Characters for the first num_sites states of seq."""
    return seq[:num_sites].astype(np.uint8).tobytes().translate(
        _state_table(state_chars)).decode("ascii")


def _datatype(state_chars):
    return "DNA" if len(state_chars) == 4 else "PROTEIN"


def write_sequences(f, fmt, tree_set, partition_lengths, num_taxa,
                    total_sites, state_chars, tree_no=0, dataset_no=0):
    """Write the tip sequences of every partition in the format code fmt."""
    names, rows = _rows_by_taxon(tree_set, partition_lengths, state_chars,
                                 num_taxa)
    if fmt == "n":
        title = ""
        if tree_no and dataset_no:
            title = f"Tree {tree_no}, Dataset {dataset_no}"
        elif tree_no:
            title = f"Tree {tree_no}"
        elif dataset_no:
            title = f"Dataset {dataset_no}"
        nexus_rows = list(zip(names, rows))
        _write_nexus_block(f, title, nexus_rows, _datatype(state_chars), total_sites)
        return

    if fmt == "f":
        for name, seq in zip(names, rows):
            f.write(f">{name}\n")
            for j in range(0, len(seq), 72):
                f.write(seq[j:j + 72] + "\n")
        return

    relaxed = fmt == "r"
    labels = None if relaxed else _phylip_labels(names)
    f.write(f"{num_taxa} {total_sites}\n")
    for idx, (name, seq) in enumerate(zip(names, rows)):
        f.write(f"{_relaxed_label(name)} " if relaxed else f"{labels[idx]:<10}")
        f.write(seq)
        f.write("\n")


def _require_unique(names, part_no):
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise ValueError(f"partition {part_no} has duplicate tip name(s) {dup}")


def _phylip_labels(names):
    """10-character labels for strict PHYLIP, refusing to merge two taxa."""
    clipped = [name[:10] for name in names]
    dup = sorted({n for n in clipped if clipped.count(n) > 1})
    if dup:
        raise ValueError(f"strict PHYLIP keeps 10 characters and these names would collide: {dup}; use -o r (relaxed) instead")
    return clipped


def _rows_by_taxon(tree_set, partition_lengths, state_chars, num_taxa):
    """Concatenated characters per taxon, matched by tip name.

    Partitions may list their tips in a different order, so the tip name is
    the only stable key for the concatenation.
    """
    names = [tree_set[0].names[i] for i in range(num_taxa)]
    _require_unique(names, 1)
    indexes = []
    for k, tree in enumerate(tree_set, start=1):
        tree_names = [tree.names[i] for i in range(tree.num_tips)]
        _require_unique(tree_names, k)
        missing = [n for n in names if n not in tree_names]
        extra = [n for n in tree_names if n not in names]
        if missing or extra:
            raise ValueError(
                f"partition {k} taxa do not match partition 1 "
                f"(missing {missing}, extra {extra}); every partition must "
                f"carry exactly the same taxa")
        indexes.append({n: i for i, n in enumerate(tree_names)})
    rows = []
    for name in names:
        rows.append("".join(
            _encode(tree.tips[index[name]].sequence, partition_lengths[k],
                    state_chars)
            for k, (tree, index) in enumerate(zip(tree_set, indexes))))
    return names, rows


def _write_nexus_block(f, title, rows, datatype, nchar):
    """Emit a NEXUS DATA block from (name, sequence) pairs."""
    f.write(f"Begin DATA; [{title}]\n" if title else "Begin DATA;\n")
    f.write(f"  Dimensions NTAX={len(rows)} NCHAR={nchar};\n")
    f.write(f"  Format MISSING=? GAP=- DATATYPE={datatype};\n")
    f.write("  Matrix\n")
    labels = [_nexus_label(n) for n, _ in rows]
    width = max(len(x) for x in labels)
    for label, (_, seq) in zip(labels, rows):
        f.write(label + " " * (width - len(label) + 1) + seq + "\n")
    f.write("  ;\nEnd;\n\n")


###----ancestral----###

def write_ancestral(f, fmt, tree, state_chars):
    """Write the sequences of every node, tips and ancestors alike.

    Internal nodes are numbered from num_tips + 1 upwards in pre-order.
    """
    write_ancestral_multi(f, fmt, [tree], state_chars)


def _node_signature(tree):
    # (is_tip, name) in the pre-order _collect_ancestral uses; auto-numbered
    # internal labels are left out because they collide with numeric tip names
    sig = []

    def visit(node):
        if node is None:
            return
        if node.tip_no == -1:
            sig.append((False, ""))
            visit(node.branch1)
            visit(node.branch2)
        else:
            sig.append((True, tree.names[node.tip_no]))

    visit(tree.root)
    if not tree.rooted:
        visit(tree.root.branch0)
    return sig


def write_ancestral_multi(f, fmt, trees, state_chars):
    """Ancestral output for one dataset, with partitions concatenated.

    Rows are concatenated in the same order the tip alignment uses. Every
    partition must carry the same topology and tip names, else this raises.
    """
    rows = _collect_ancestral(trees[0], state_chars)
    if len(trees) > 1:
        merged = [[label, seq] for label, seq in rows]
        base_sig = _node_signature(trees[0])
        for k, tree in enumerate(trees[1:], start=2):
            if _node_signature(tree) != base_sig:
                raise ValueError(
                    "ancestral output with multiple partitions needs the "
                    f"same topology in every partition; partition {k} "
                    "differs from partition 1")
            part = _collect_ancestral(tree, state_chars)
            for j, (_, seq) in enumerate(part):
                merged[j][1] += seq
        rows = [(label, seq) for label, seq in merged]
    num_sites = len(rows[0][1]) if rows else 0

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
        labels = None if relaxed else _phylip_labels([n for n, _ in rows])
        for idx, (name, seq) in enumerate(rows):
            f.write(f"{_relaxed_label(name)} {seq}\n" if relaxed
                    else f"{labels[idx]:<10}{seq}\n")


def _collect_ancestral(tree, state_chars):
    """(label, sequence) for every node, root first, in pre-order."""
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
    labels = [label for label, _ in rows]
    dup = sorted({l for l in labels if labels.count(l) > 1})
    if dup:
        raise ValueError(
            f"ancestral node labels collide ({dup}): internal nodes are numbered from num_tips+1, so a tip with that numeric name clashes")
    return rows


###----rate table----###

def write_rates(f, site_rates, num_sites):
    """The per-site relative rate table behind -wr."""
    f.write("Relative rates for each site:\n")
    for i in range(num_sites):
        f.write(f"{i + 1}\t{site_rates[i]:.6f}\n")
    f.write("\n")
