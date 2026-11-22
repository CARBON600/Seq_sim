"""Reader for BEAST/MrBayes-style NEXUS tree files.
Pulls out the Newick strings, the TRANSLATE map and the [&R]/[&U] rootedness
hints, while seq_sim's own [f=N]/[r=N] annotations survive.
input: NEXUS text
output: Newick strings, a TRANSLATE map and a rootedness flag per tree
"""

import re

_FR_TOKEN = re.compile(r"^[fr]\s*=\s*[+-]?\d+$", re.IGNORECASE)
# A bare number in brackets is seq_sim's per-branch rate multiplier in Newick
# ("A:0.1[2.0]"), so it survives the NEXUS strip.
_NUM_TOKEN = re.compile(r"^[+-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|nan|inf|infinity)$", re.IGNORECASE)
_TREE_CMD = re.compile(r"(?i)(?:TREE|UTREE)\s+\*?\s*(\S+)\s*=")


###function definition

def _quote_mask(text):
    """True for characters outside a single-quoted NEXUS label."""
    mask = []
    in_quote = False
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "'":
            if in_quote and i + 1 < len(text) and text[i + 1] == "'":
                mask.extend([False, False])
                i += 2
                continue
            in_quote = not in_quote
            mask.append(False)
            i += 1
            continue
        mask.append(not in_quote)
        i += 1
    return mask

def _unquote(text):
    """Strip one layer of NEXUS quoting from a label."""
    t = text.strip()
    if len(t) >= 2 and t[0] == "'" and t[-1] == "'":
        return t[1:-1].replace("''", "'")
    return t

def _split_quoted(text, sep=","):
    """Split on sep only where it is outside a quoted label."""
    mask = _quote_mask(text)
    parts = []
    start = 0
    for idx, ch in enumerate(text):
        if ch == sep and mask[idx]:
            parts.append(text[start:idx])
            start = idx + 1
    parts.append(text[start:])
    return parts


def is_nexus(text):
    """Whether text looks like NEXUS rather than plain Newick."""
    return bool(re.search(r"(?im)^\s*#NEXUS|^\s*BEGIN\s+\w+\s*;", text))


def _is_kept_annotation(content):
    """Whether every token is something seq_sim's Newick parser understands.

    A bracket may mix them: [f=1,2.0] sets composition 1 and a rate multiplier
    of 2.0, so the whole bracket is kept rather than dropped.
    """
    tokens = [t.strip() for t in content.split(",") if t.strip()]
    if not tokens:
        return False
    return all(_FR_TOKEN.match(t) or _NUM_TOKEN.match(t) for t in tokens)


def _strip_bracket_content(text, keep_fr=False):
    """Remove [...] comments, tracking nesting and quoted labels.

    A '[' inside a quoted NEXUS label ('A[west]') is part of the name, not a
    comment, so quotes are tracked while scanning.
    """
    mask = _quote_mask(text)
    result = []
    depth = 0
    buf = []
    for idx, ch in enumerate(text):
        if mask[idx] and (ch == "[" or ch == "]"):
            if ch == "[":
                if depth == 0:
                    buf = []
                depth += 1
            elif depth > 0:
                depth -= 1
                if depth == 0 and keep_fr and _is_kept_annotation("".join(buf)):
                    result.append("[" + "".join(buf) + "]")
            continue
        if depth == 0:
            result.append(ch)
        elif depth == 1:
            buf.append(ch)
    return "".join(result)

def _statements(text):
    """Split NEXUS text at top-level semicolons, ignoring quoted labels."""
    mask = _quote_mask(text)
    out = []
    start = 0
    depth = 0
    for i, ch in enumerate(text):
        if not mask[i]:
            continue
        if ch == "[":
            depth += 1
        elif ch == "]":
            if depth:
                depth -= 1
        elif ch == ";" and depth == 0:
            out.append(text[start:i])
            start = i + 1
    tail = text[start:]
    if tail.strip():
        out.append(tail)
    return out


def _trees_blocks(text):
    """Every BEGIN TREES ... END block, as a ';'-joined statement string."""
    blocks = []
    cur = None
    for stmt in _statements(text):
        words = [w.upper() for w in stmt.strip().split()]
        if "BEGIN" in words and words[words.index("BEGIN") + 1:][:1] == ["TREES"]:
            cur = []
            continue
        head = words[0] if words else ""
        if cur is not None and head in ("END", "ENDBLOCK"):
            blocks.append(";".join(cur))
            cur = None
            continue
        if cur is not None:
            cur.append(stmt)
    if cur is not None:
        blocks.append(";".join(cur))
    return blocks


def _parse_translate(block):
    """The TRANSLATE command as {int_id: name}, quote aware."""
    translate = {}
    for stmt in _statements(block):
        if not stmt.strip().upper().startswith("TRANSLATE"):
            continue
        body = stmt.strip()[len("TRANSLATE"):]
        for item in _split_quoted(body, ","):
            item = item.strip()
            if not item:
                continue
            parts = item.split(None, 1)
            if len(parts) == 2:
                try:
                    translate[int(parts[0])] = _unquote(parts[1])
                except ValueError:
                    pass
    return translate

def _extract_tree_newicks(block):
    """Every (newick, rooted) pair in a TREES block."""
    results = []
    clean = _strip_bracket_content(block, keep_fr=True)
    mask = _quote_mask(clean)

    for m in _TREE_CMD.finditer(clean):
        tree_name = m.group(1)
        pos = m.end()

        paren_start = next((i for i in range(pos, len(clean)) if clean[i] == "(" and mask[i]), -1)
        if paren_start == -1:
            continue

        depth = 0
        paren_end = None
        for i in range(paren_start, len(clean)):
            if not mask[i]:
                continue
            if clean[i] == "(":
                depth += 1
            elif clean[i] == ")":
                depth -= 1
                if depth == 0:
                    paren_end = i + 1
                    break
        if paren_end is None:
            continue

        semi_pos = next((i for i in range(paren_end, len(clean)) if clean[i] == ";" and mask[i]), len(clean) - 1)

        # Keep a partition prefix ("[100,1](...)") that sits just before the
        # Newick; _strip_bracket_content has already removed [&R] comments.
        prefix_start = paren_start
        j = paren_start - 1
        while j >= 0 and clean[j].isspace():
            j -= 1
        if j >= 0 and clean[j] == "]":
            depth = 0
            k = j
            while k >= 0:
                if clean[k] == "]":
                    depth += 1
                elif clean[k] == "[":
                    depth -= 1
                    if depth == 0:
                        break
                k -= 1
            if k >= 0:
                prefix_start = k
        text = clean[prefix_start:semi_pos + 1].rstrip()
        if not text.endswith(";"):
            text += ";"
        results.append((text, _detect_rooted(block, tree_name)))
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
    mask = _quote_mask(original_block)
    while pos < len(original_block):
        ch = original_block[pos]
        if not mask[pos]:
            raw.append(ch)
            pos += 1
            continue
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


def parse_nexus_tree_entries(text):
    # One entry per TREE: (newick, rooted, translate-of-its-own-block).
    # TRANSLATE is block-scoped, so each entry carries its own map.
    entries = []
    blocks = _trees_blocks(text)
    for block in blocks:
        block_translate = _parse_translate(block)
        for newick, rooted in _extract_tree_newicks(block):
            entries.append((newick, rooted, block_translate))
    if not entries:
        raise ValueError("No trees found in NEXUS file")
    return entries

def parse_nexus_trees(text):
    """Parse NEXUS text into (translate_map, [(newick, rooted), ...]).

    Raises ValueError when the text holds no trees.
    """
    translate = {}
    all_trees = []

    blocks = _trees_blocks(text)
    for block in blocks:
        translate.update(_parse_translate(block))
        all_trees.extend(_extract_tree_newicks(block))

    if not all_trees:
        raise ValueError("No trees found in NEXUS file")
    return translate, all_trees
