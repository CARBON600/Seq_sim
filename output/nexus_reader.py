"""Reader for BEAST/MrBayes-style NEXUS tree files.

Pulls out the Newick strings, the TRANSLATE map and the [&R]/[&U] rootedness
hints, while letting seq_sim's own [f=N]/[r=N] annotations survive.
"""

import re

_FR_TOKEN = re.compile(r"^[fr]\s*=\s*-?\d+$", re.IGNORECASE)
_TREE_CMD = re.compile(r"(?i)(?:TREE|UTREE)\s+\*?\s*(\S+)\s*=")


def is_nexus(text):
    """Whether *text* looks like NEXUS rather than plain Newick."""
    return bool(re.search(r"(?im)^\s*#NEXUS|^\s*BEGIN\s+\w+\s*;", text))


def _is_fr_annotation(content):
    """Whether *content* is nothing but f=N / r=N tokens."""
    tokens = [t.strip() for t in content.split(",") if t.strip()]
    if not tokens:
        return False
    return all(_FR_TOKEN.match(t) for t in tokens)


def _strip_bracket_content(text, keep_fr=False):
    """Remove [...] comments, tracking nesting.

    With *keep_fr*, a top-level annotation made only of f=/r= tokens (``[f=1]``,
    ``[f=1,r=2]``) is kept so per-branch model indices survive; everything else
    (``[&R]``, ``[&height=...]``) is dropped.
    """
    result = []
    depth = 0
    buf = []
    for ch in text:
        if ch == "[":
            if depth == 0:
                buf = []
            depth += 1
        elif ch == "]":
            if depth > 0:
                depth -= 1
                if depth == 0 and keep_fr and _is_fr_annotation("".join(buf)):
                    result.append("[" + "".join(buf) + "]")
        elif depth == 0:
            result.append(ch)
        elif depth == 1:
            buf.append(ch)
    return "".join(result)


def _parse_translate(block):
    """The TRANSLATE command as {int_id: name}."""
    translate = {}
    m = re.search(r"(?is)TRANSLATE\s+(.+?);", block)
    if not m:
        return translate
    body = _strip_bracket_content(m.group(1))
    for item in body.split(","):
        item = item.strip()
        if not item:
            continue
        parts = item.split(maxsplit=1)
        if len(parts) == 2:
            try:
                translate[int(parts[0])] = parts[1]
            except ValueError:
                pass
    return translate


def _extract_tree_newicks(block):
    """Every (newick, rooted) pair in a TREES block."""
    results = []
    clean = _strip_bracket_content(block, keep_fr=True)

    for m in _TREE_CMD.finditer(clean):
        tree_name = m.group(1)
        pos = m.end()

        paren_start = clean.find("(", pos)
        if paren_start == -1:
            continue

        depth = 0
        paren_end = None
        for i in range(paren_start, len(clean)):
            if clean[i] == "(":
                depth += 1
            elif clean[i] == ")":
                depth -= 1
                if depth == 0:
                    paren_end = i + 1
                    break
        if paren_end is None:
            continue

        semi_pos = clean.find(";", paren_end)
        if semi_pos == -1:
            continue

        results.append((clean[paren_start:semi_pos + 1],
                        _detect_rooted(block, tree_name)))
    return results


def _detect_rooted(original_block, tree_name):
    """True/False for an [&R]/[&U] sitting between '=' and the first '('."""
    pattern = r"(?i)(?:TREE|UTREE)\s+\*?\s*" + re.escape(tree_name) + r"\s*="
    m = re.search(pattern, original_block)
    if not m:
        return None

    pos = m.end()
    depth = 0
    raw = []
    while pos < len(original_block):
        ch = original_block[pos]
        if ch == "[":
            depth += 1
        elif ch == "]":
            if depth > 0:
                depth -= 1
        elif ch == "(" and depth == 0:
            break
        raw.append(ch)
        pos += 1

    prefix = "".join(raw)
    if "[&R]" in prefix:
        return True
    if "[&U]" in prefix:
        return False
    return None


def parse_nexus_trees(text):
    """Parse NEXUS text into (translate_map, [(newick, rooted), ...]).

    Raises ValueError when the text holds no trees.
    """
    translate = {}
    all_trees = []

    blocks = re.findall(r"(?is)BEGIN\s+TREES\s*;(.*?)(?:END|ENDBLOCK)\s*;", text)
    for block in blocks:
        translate.update(_parse_translate(block))
        all_trees.extend(_extract_tree_newicks(block))

    if not all_trees:
        raise ValueError("No trees found in NEXUS file")
    return translate, all_trees
