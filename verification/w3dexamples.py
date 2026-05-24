#!/usr/bin/env python3
"""
w3dexamples.py

Search and catalog WARP3D verification examples.

This script scans a WARP3D verification directory, identifies top-level
verification input files, records files referenced by

    *input from filename

and builds a searchable JSON index plus a Markdown catalog.

Important WARP3D policy used by this tool
-----------------------------------------

Files referenced by "*input from ..." are treated as support/include files.
They are recorded in the index, but they are not treated as independent
verification examples and are not searched by default.

Feature terms
-------------

Feature names come from a required simple plain-text term file, for example:

    @element
    q3disop | hex20
    l3disop | hex8

    @material
    mises | j2 plasticity | von mises

The first item on each line is the canonical term. Items after "|" are aliases.
The generated feature tag is:

    category_canonical_term

Examples:

    q3disop       -> element_q3disop
    mises         -> material_mises
    tied contact  -> constraint_tie_mesh_constraints, if defined that way

Typical use
-----------

    python3 w3dexamples.py index /path/to/verification --terms z_features.txt

    python3 w3dexamples.py find "mises q3disop domain"

    python3 w3dexamples.py find --debug "tie-mesh bilinear"

    python3 w3dexamples.py show path_fragment

    python3 w3dexamples.py tags

    python3 w3dexamples.py copy path_fragment ./scratch_example

Python version
--------------

This script targets Python 3.8 or newer.
"""

import argparse
import json
import os
import re
import shlex
import shutil
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_INDEX_FILE = "w3dexamples_index.json"
DEFAULT_CATALOG_FILE = "w3dexamples_catalog.md"

# Candidate top-level verification file names.
# These are deliberately broad. A later pass removes files referenced
# by "*input from ..." so include files do not become search results.
MAIN_FILE_PATTERNS = [
    "test*",
    "*.inp",
    "*.in",
    "*.dat",
    "*.w3d",
    "*.warp3d",
    "input*",
]

# Files to ignore even if their names match MAIN_FILE_PATTERNS.
EXCLUDE_FILE_NAMES = {
    "driver.inp",
}

# Filename prefixes to ignore.
EXCLUDE_FILE_PREFIXES = (
    "clean",
)

# Directory names to skip during the tree walk.
SKIP_DIR_NAMES = {
    ".git",
    ".svn",
    "__pycache__",
    "build",
    "tmp",
    "temp",
    "results",
    "output",
    "outputs",
}

# Include syntax recognized by this tool:
#
#     *input from filename
#     *input from "filename with spaces"
#     *input from 'filename with spaces'
#
INCLUDE_RE = re.compile(
    r"^\s*\*\s*input\s+from\s+(?P<name>.+?)\s*$",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Small utility functions
# ---------------------------------------------------------------------------

def read_text_file(path):
    """Read a text file using UTF-8 first, then Latin-1 as a fallback."""
    encodings = ["utf-8", "latin-1"]

    for enc in encodings:
        try:
            return path.read_text(encoding=enc).splitlines()
        except UnicodeDecodeError:
            pass

    return path.read_text(errors="replace").splitlines()


def slugify(text):
    """
    Convert a term to a safe lowercase tag component.

    Examples:
        "tie-mesh constraints" -> "tie_mesh_constraints"
        "J-integral"           -> "j_integral"
    """
    s = text.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    if not s:
        s = "unknown"
    return s


def term_to_regex(term):
    """
    Convert a literal search term into a regex.

    Hyphen, underscore, and whitespace are treated as equivalent separators.
    For example, "tied contact" matches:

        tied contact
        tied-contact
        tied_contact
    """
    pieces = re.split(r"[\s_-]+", term.strip())
    pieces = [p for p in pieces if p]

    if not pieces:
        return r"a^"       # regex that never matches

    escaped = [re.escape(p) for p in pieces]
    body = r"[\s_-]+".join(escaped)

    # These guards prevent "tet" from matching "tet10".
    return r"(?<![A-Za-z0-9])" + body + r"(?![A-Za-z0-9])"


def strip_quotes(s):
    """Remove matching single or double quotes around a string."""
    s = s.strip()

    if len(s) >= 2:
        if s[0] == '"' and s[-1] == '"':
            return s[1:-1]
        if s[0] == "'" and s[-1] == "'":
            return s[1:-1]

    return s


def is_comment_line(line):
    """
    WARP3D verification-file comment convention used by this tool:

        column 1 is c, C, #, or !
        column 2 is a blank space
    """
    if len(line) < 2:
        return False

    return line[0] in ("c", "C", "#", "!") and line[1] == " "


def comment_text(line):
    """Return the text after the two-character comment prefix."""
    return line[2:].strip()


# ---------------------------------------------------------------------------
# Feature-term file handling
# ---------------------------------------------------------------------------

def load_feature_terms(term_file):
    """
    Load feature terms from a plain-text file.

    Format:

        @category
        canonical term | alias 1 | alias 2

    The returned list contains dictionaries with these fields:

        tag
        family
        canonical
        aliases
        patterns
    """
    if term_file is None:
        raise RuntimeError("--terms feature file is required")

    term_file = Path(term_file)

    if not term_file.exists():
        raise RuntimeError("feature term file not found: {}".format(term_file))

    terms = []
    current_family = None

    for line_number, raw in enumerate(read_text_file(term_file), start=1):
        line = raw.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        if line.startswith("@"):
            current_family = slugify(line[1:])
            continue

        if current_family is None:
            raise RuntimeError(
                "{}:{}: feature term appears before first @category".format(
                    term_file, line_number
                )
            )

        pieces = [p.strip() for p in line.split("|") if p.strip()]
        if not pieces:
            continue

        canonical = pieces[0]
        aliases = pieces[1:]
        tag = current_family + "_" + slugify(canonical)

        patterns = []
        for piece in pieces:
            patterns.append(term_to_regex(piece))

        terms.append(
            {
                "tag": tag,
                "family": current_family,
                "canonical": canonical,
                "aliases": aliases,
                "patterns": patterns,
            }
        )

    return terms


# ---------------------------------------------------------------------------
# Verification-file discovery and include scanning
# ---------------------------------------------------------------------------

def should_skip_dir(path):
    """Return True if any component of path is a directory we want to skip."""
    for part in path.parts:
        if part in SKIP_DIR_NAMES:
            return True

    return False


def looks_like_main_file(path):
    """
    Return True if path looks like a possible top-level verification input.

    This is only the first pass. Files included by another input file are
    removed later.
    """
    if not path.is_file():
        return False

    name = path.name.lower()

    if name.startswith("."):
        return False

    if name in EXCLUDE_FILE_NAMES:
        return False

    for prefix in EXCLUDE_FILE_PREFIXES:
        if name.startswith(prefix):
            return False

    bad_suffixes = {
        ".out",
        ".log",
        ".msg",
        ".sta",
        ".bak",
        ".old",
        ".orig",
        ".pyc",
    }

    if path.suffix.lower() in bad_suffixes:
        return False

    for pattern in MAIN_FILE_PATTERNS:
        if path.match(pattern):
            return True

    return False


def parse_include_filename(line):
    """
    Return the filename in a WARP3D '*input from ...' line.

    Return None if the line is not an include command.
    """
    match = INCLUDE_RE.match(line)
    if not match:
        return None

    name = strip_quotes(match.group("name"))

    # Remove simple trailing shell-style comments.
    name = re.split(r"\s+[!#].*$", name)[0].strip()

    if name:
        return name

    return None


def resolve_include(include_name, current_file, root):
    """
    Resolve an include file path.

    Try:
        1. absolute path if supplied
        2. relative to the file containing the include command
        3. relative to the verification root

    If the file does not exist, return the most likely relative path.
    """
    include_path = Path(include_name)

    if include_path.is_absolute():
        return include_path

    candidate = current_file.parent / include_path
    if candidate.exists():
        return candidate

    candidate = root / include_path
    if candidate.exists():
        return candidate

    return current_file.parent / include_path


def discover_candidate_main_files(root):
    """Find files whose names look like verification inputs."""
    files = []

    for dirpath, dirnames, filenames in os.walk(root):
        directory = Path(dirpath)

        # Prevent os.walk from descending into skipped directories.
        keep_dirnames = []
        for d in dirnames:
            if d not in SKIP_DIR_NAMES:
                keep_dirnames.append(d)
        dirnames[:] = keep_dirnames

        if should_skip_dir(directory):
            continue

        for filename in filenames:
            path = directory / filename
            if looks_like_main_file(path):
                files.append(path.resolve())

    return sorted(set(files))


def collect_included_paths(path, root, included_paths, visited_stack=None):
    """
    Recursively collect files referenced by '*input from ...'.

    Any file reached by an include command is a support file. It should not
    appear as an independent search result.
    """
    if visited_stack is None:
        visited_stack = []

    path = path.resolve()

    if path in visited_stack:
        return

    if not path.exists():
        return

    new_stack = visited_stack + [path]

    for line in read_text_file(path):
        include_name = parse_include_filename(line)
        if include_name is None:
            continue

        include_path = resolve_include(include_name, path, root).resolve()
        included_paths.add(include_path)

        if include_path.exists():
            collect_included_paths(
                include_path,
                root,
                included_paths,
                visited_stack=new_stack,
            )


def remove_included_files_from_candidates(candidates, root):
    """
    Remove any candidate file referenced by '*input from ...'.

    This prevents support files such as get_j.inp from appearing in find
    results as if they were independent verification problems.
    """
    included_paths = set()

    for candidate in candidates:
        collect_included_paths(candidate, root, included_paths)

    main_files = []
    for candidate in candidates:
        if candidate.resolve() not in included_paths:
            main_files.append(candidate.resolve())

    return sorted(main_files)


def scan_file_for_includes(path, root, include_records, visited_stack=None, depth=0):
    """
    Record include commands reachable from path.

    Include-file contents are not added to the searchable/indexed text.
    """
    if visited_stack is None:
        visited_stack = []

    path = path.resolve()

    if path in visited_stack:
        include_records.append(
            {
                "include_file": str(path),
                "included_from": str(visited_stack[-1]) if visited_stack else "",
                "line_number": 0,
                "depth": depth,
                "status": "cycle-detected",
            }
        )
        return

    if not path.exists():
        include_records.append(
            {
                "include_file": str(path),
                "included_from": str(visited_stack[-1]) if visited_stack else "",
                "line_number": 0,
                "depth": depth,
                "status": "missing",
            }
        )
        return

    lines = read_text_file(path)
    new_stack = visited_stack + [path]

    for line_number, line in enumerate(lines, start=1):
        include_name = parse_include_filename(line)
        if include_name is None:
            continue

        include_path = resolve_include(include_name, path, root).resolve()

        if include_path.exists():
            status = "ok"
        else:
            status = "missing"

        include_records.append(
            {
                "include_file": str(include_path),
                "included_from": str(path),
                "line_number": line_number,
                "depth": depth + 1,
                "status": status,
            }
        )

        scan_file_for_includes(
            include_path,
            root,
            include_records,
            visited_stack=new_stack,
            depth=depth + 1,
        )


# ---------------------------------------------------------------------------
# Index creation and feature detection
# ---------------------------------------------------------------------------

def read_main_file_lines(path):
    """
    Return records for the top-level file only.

    Each record is:

        {
            "source_file": Path(...),
            "line_number": integer,
            "text": line text
        }
    """
    records = []

    for line_number, line in enumerate(read_text_file(path), start=1):
        records.append(
            {
                "source_file": path.resolve(),
                "line_number": line_number,
                "text": line,
            }
        )

    return records


def infer_title(indexed_lines, fallback):
    """Infer a short description from early WARP3D comment lines."""
    candidates = []

    for rec in indexed_lines[:80]:
        line = rec["text"]

        if not line.strip():
            continue

        if is_comment_line(line):
            text = comment_text(line)
            if len(text) >= 12:
                candidates.append(text)

        elif len(candidates) >= 1:
            break

    if candidates:
        return candidates[0][:120]

    return fallback


def detect_features(indexed_lines, root, feature_terms):
    """
    Detect features in the top-level file text.

    Return:
        features          sorted list of feature tags
        feature_families  dictionary: family -> list of feature tags
        evidence          list of evidence dictionaries
    """
    compiled = {}
    tag_to_family = {}

    for term in feature_terms:
        compiled_patterns = []
        for pattern in term["patterns"]:
            compiled_patterns.append(re.compile(pattern, re.IGNORECASE))

        compiled[term["tag"]] = compiled_patterns
        tag_to_family[term["tag"]] = term["family"]

    features = set()
    family_sets = {}
    evidence = []
    evidence_count = {}

    max_evidence_per_pattern = 5

    for rec in indexed_lines:
        source_file = rec["source_file"]
        line_number = rec["line_number"]
        line = rec["text"]

        for feature in compiled:
            family = tag_to_family.get(feature, "unknown")

            for pattern in compiled[feature]:
                if not pattern.search(line):
                    continue

                features.add(feature)

                if family not in family_sets:
                    family_sets[family] = set()
                family_sets[family].add(feature)

                key = (feature, pattern.pattern)
                count = evidence_count.get(key, 0)

                if count >= max_evidence_per_pattern:
                    continue

                try:
                    rel_source = str(source_file.resolve().relative_to(root.resolve()))
                except ValueError:
                    rel_source = str(source_file)

                evidence.append(
                    {
                        "feature": feature,
                        "family": family,
                        "pattern": pattern.pattern,
                        "source_file": rel_source,
                        "line_number": line_number,
                        "line_text": line.strip()[:240],
                    }
                )

                evidence_count[key] = count + 1

    feature_families = {}
    for family in sorted(family_sets):
        feature_families[family] = sorted(family_sets[family])

    return sorted(features), feature_families, evidence


def build_index(root, feature_terms):
    """Build the full verification-problem index."""
    root = root.resolve()

    candidates = discover_candidate_main_files(root)
    main_files = remove_included_files_from_candidates(candidates, root)

    records = []

    for main_file in main_files:
        indexed_lines = read_main_file_lines(main_file)

        include_records = []
        scan_file_for_includes(main_file, root, include_records)

        features, feature_families, evidence = detect_features(
            indexed_lines,
            root,
            feature_terms,
        )

        indexed_text = "\n".join(rec["text"] for rec in indexed_lines)
        rel_main = str(main_file.resolve().relative_to(root))
        directory = str(main_file.parent.resolve().relative_to(root))
        title = infer_title(indexed_lines, fallback=main_file.stem)

        records.append(
            {
                "main_file": str(main_file.resolve()),
                "rel_main_file": rel_main,
                "directory": directory,
                "title": title,
                "includes": include_records,
                "features": features,
                "feature_families": feature_families,
                "evidence": evidence,
                "indexed_line_count": len(indexed_lines),
                "indexed_text": indexed_text,
            }
        )

    return records


# ---------------------------------------------------------------------------
# Index reading/writing and Markdown catalog
# ---------------------------------------------------------------------------

def write_json_index(records, index_file, root, feature_terms):
    """Write the JSON index."""
    payload = {
        "schema_version": 3,
        "root": str(root.resolve()),
        "n_problems": len(records),
        "feature_terms": feature_terms,
        "problems": records,
    }

    index_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def read_json_index(index_file):
    """Read a JSON index created by this tool."""
    payload = json.loads(index_file.read_text(encoding="utf-8"))

    root = Path(payload["root"])
    records = payload.get("problems", [])

    return root, records


def write_markdown_catalog(records, catalog_file):
    """Write a human-readable Markdown catalog."""
    lines = []
    lines.append("# WARP3D Verification Problem Catalog")
    lines.append("")
    lines.append("Number of indexed problems: **{}**".format(len(records)))
    lines.append("")

    by_dir = {}
    for record in records:
        directory = record.get("directory", ".")
        if directory not in by_dir:
            by_dir[directory] = []
        by_dir[directory].append(record)

    for directory in sorted(by_dir):
        lines.append("## {}".format(directory))
        lines.append("")
        lines.append("| Problem | Title / inferred description | Features | Includes |")
        lines.append("|---|---|---|---:|")

        records_in_dir = sorted(
            by_dir[directory],
            key=lambda r: r.get("rel_main_file", ""),
        )

        for record in records_in_dir:
            features = record.get("features", [])
            if features:
                feature_text = ", ".join(features)
            else:
                feature_text = "—"

            include_count = 0
            for inc in record.get("includes", []):
                if inc.get("status") == "ok":
                    include_count += 1

            title = record.get("title", "").replace("|", "\\|")
            rel_main = record.get("rel_main_file", "")

            lines.append(
                "| `{}` | {} | {} | {} |".format(
                    rel_main,
                    title,
                    feature_text,
                    include_count,
                )
            )

        lines.append("")

    lines.append("# Feature Index")
    lines.append("")

    feature_map = {}
    family_map = {}

    for record in records:
        for feature in record.get("features", []):
            if feature not in feature_map:
                feature_map[feature] = []
            feature_map[feature].append(record.get("rel_main_file", ""))

        for family, tags in record.get("feature_families", {}).items():
            if family not in family_map:
                family_map[family] = set()
            for tag in tags:
                family_map[family].add(tag)

    for family in sorted(family_map):
        lines.append("## {}".format(family))
        lines.append("")

        for feature in sorted(family_map[family]):
            lines.append("### {}".format(feature))
            lines.append("")

            for problem in sorted(feature_map.get(feature, [])):
                lines.append("- `{}`".format(problem))

            lines.append("")

    catalog_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Search and display
# ---------------------------------------------------------------------------

def split_query(query):
    """Split a user query while respecting quotes."""
    if not query:
        return []

    try:
        return shlex.split(query)
    except ValueError:
        # Fall back to simple splitting if the user has unmatched quotes.
        return query.split()


def query_term_matches_text(term, text):
    """Return True if term appears in text, with separator normalization."""
    pattern = re.compile(term_to_regex(term), re.IGNORECASE)
    return pattern.search(text) is not None


def count_query_term_matches(term, text):
    """Count occurrences of term in text, with separator normalization."""
    pattern = re.compile(term_to_regex(term), re.IGNORECASE)
    return len(list(pattern.finditer(text)))


def record_has_required_features(record, required_features):
    """Check explicit --feature filters."""
    feature_set = set(record.get("features", []))

    for required in required_features:
        required_slug = slugify(required)

        matched = False

        if required in feature_set:
            matched = True

        if required_slug in feature_set:
            matched = True

        for feature in feature_set:
            if feature.endswith("_" + required_slug):
                matched = True

        if not matched:
            return False

    return True


def score_record(record, query_terms, required_features):
    """
    Score a record for a find query.

    Important search rule:
        Each query term must be supported by the top-level file path, inferred
        title, or indexed top-level file text.

    The search does not match a problem merely because the query text appears
    inside a generated feature tag.
    """
    if not record_has_required_features(record, required_features):
        return -1

    score = 0
    score += 30 * len(required_features)

    path_text = record.get("rel_main_file", "")
    title_text = record.get("title", "")
    indexed_text = record.get("indexed_text", "")

    for term in query_terms:
        term_score = 0

        if query_term_matches_text(term, path_text):
            term_score += 12

        if query_term_matches_text(term, title_text):
            term_score += 10

        count = count_query_term_matches(term, indexed_text)
        if count > 0:
            term_score += min(count, 10)

        # All query terms are required. This keeps searches such as
        # "tie-mesh bilinear" tight and predictable.
        if term_score == 0:
            return -1

        score += term_score

    return score


def evidence_matches_query(ev, query_terms, required_features):
    """Return True if an evidence line helps explain the current query."""
    feature = ev.get("feature", "")
    feature_slug = slugify(feature)

    for required in required_features:
        required_slug = slugify(required)
        if feature == required or feature_slug.endswith("_" + required_slug):
            return True

    line = ev.get("line_text", "")
    for term in query_terms:
        if query_term_matches_text(term, line):
            return True

    return False


def matched_feature_tags(record, query_terms, required_features):
    """
    Return features directly associated with the current find request.

    Normal find output should not list every feature detected in a broad
    verification example. It should list only the features that explain why
    this result matched.
    """
    matched = set()

    for ev in record.get("evidence", []):
        if evidence_matches_query(ev, query_terms, required_features):
            matched.add(ev.get("feature", ""))

    # Explicit --feature requests should appear even if no evidence line was
    # kept due to evidence limits.
    feature_set = set(record.get("features", []))
    for required in required_features:
        required_slug = slugify(required)
        for feature in feature_set:
            if feature == required or feature.endswith("_" + required_slug):
                matched.add(feature)

    return sorted([m for m in matched if m])


def format_wrapped_list(label, values, per_line=4, indent="   "):
    """Return formatted lines for a label plus a wrapped comma-separated list."""
    if not values:
        return [indent + label + ": —"]

    lines = []
    first_prefix = indent + label + ": "
    continuation_prefix = " " * len(first_prefix)

    for i in range(0, len(values), per_line):
        chunk = values[i:i + per_line]
        text = ", ".join(chunk)

        if i == 0:
            lines.append(first_prefix + text)
        else:
            lines.append(continuation_prefix + text)

    return lines


def select_evidence(record, query_terms, required_features):
    """Select evidence lines to show in --debug output."""
    selected = []

    for ev in record.get("evidence", []):
        if evidence_matches_query(ev, query_terms, required_features):
            selected.append(ev)

    if not selected:
        selected = record.get("evidence", [])[:5]

    return selected


def print_find_results(records, query, required_features, limit, debug=False):
    """Print user-facing find results."""
    query_terms = split_query(query)

    scored = []
    for record in records:
        score = score_record(record, query_terms, required_features)

        if score >= 0:
            if score > 0 or required_features or not query_terms:
                scored.append((score, record))

    scored.sort(key=lambda pair: (-pair[0], pair[1].get("rel_main_file", "")))
    scored = scored[:limit]

    if not scored:
        print("No matching verification problems found.")
        return

    if len(scored) == 1:
        noun = "problem"
    else:
        noun = "problems"

    print("Found {} matching verification {}.\n".format(len(scored), noun))

    for i, pair in enumerate(scored, start=1):
        score, record = pair

        print("{}. Location: {}".format(i, record.get("rel_main_file", "")))

        title = record.get("title", "")
        if title:
            print("   Description: {}".format(title))

        matched_features = matched_feature_tags(
            record,
            query_terms,
            required_features,
        )

        for line in format_wrapped_list("Matched features", matched_features):
            print(line)

        if debug:
            all_features = record.get("features", [])
            for line in format_wrapped_list("All detected features", all_features):
                print(line)

            print("   Score: {}".format(score))

            ok_count = 0
            bad_count = 0
            for inc in record.get("includes", []):
                if inc.get("status") == "ok":
                    ok_count += 1
                else:
                    bad_count += 1

            if bad_count:
                print("   Includes: {} ok, {} missing/problematic".format(ok_count, bad_count))
            else:
                print("   Includes: {} ok".format(ok_count))

            evidence = select_evidence(record, query_terms, required_features)
            if evidence:
                print("   Evidence:")
                for ev in evidence[:8]:
                    print(
                        "      {}:{}: [{}] {}".format(
                            ev.get("source_file", ""),
                            ev.get("line_number", ""),
                            ev.get("feature", ""),
                            ev.get("line_text", ""),
                        )
                    )

        print("")


def show_record(records, target):
    """Show a detailed card for one indexed problem."""
    target_lower = target.lower()
    matches = []

    for record in records:
        rel_path = record.get("rel_main_file", "").lower()
        full_path = record.get("main_file", "").lower()

        if target_lower in rel_path or target_lower in full_path:
            matches.append(record)

    if not matches:
        print("No indexed problem matches: {}".format(target))
        return

    if len(matches) > 1:
        print("Multiple matches for {}:\n".format(target))
        for record in matches[:25]:
            print("  {}".format(record.get("rel_main_file", "")))
        return

    record = matches[0]

    print("Problem: {}".format(record.get("rel_main_file", "")))
    print("Title:   {}".format(record.get("title", "")))
    print("Dir:     {}".format(record.get("directory", "")))
    print("Indexed top-level lines: {}".format(record.get("indexed_line_count", 0)))

    for line in format_wrapped_list(
        "Features",
        record.get("features", []),
        per_line=4,
        indent="",
    ):
        print(line)

    if record.get("feature_families"):
        print("")
        print("Feature families:")
        for family in sorted(record["feature_families"]):
            tags = record["feature_families"][family]
            for line in format_wrapped_list(family, tags, per_line=4, indent="  "):
                print(line)

    print("")
    print("Includes:")
    includes = record.get("includes", [])

    if not includes:
        print("  —")
    else:
        for inc in includes:
            print(
                "  {:14s} depth={:<2d} from {}:{} -> {}".format(
                    inc.get("status", ""),
                    inc.get("depth", 0),
                    inc.get("included_from", ""),
                    inc.get("line_number", 0),
                    inc.get("include_file", ""),
                )
            )

    print("")
    print("Evidence:")
    evidence = record.get("evidence", [])
    if not evidence:
        print("  —")
    else:
        for ev in evidence[:60]:
            print(
                "  {}:{}: [{}/{}] {}".format(
                    ev.get("source_file", ""),
                    ev.get("line_number", ""),
                    ev.get("family", ""),
                    ev.get("feature", ""),
                    ev.get("line_text", ""),
                )
            )


def print_tags(records):
    """List detected feature tags and counts."""
    counts = {}
    families = {}

    for record in records:
        for family, tags in record.get("feature_families", {}).items():
            for tag in tags:
                families[tag] = family

        for feature in record.get("features", []):
            counts[feature] = counts.get(feature, 0) + 1

    if not counts:
        print("No feature tags found.")
        return

    width = max(len(tag) for tag in counts)

    sorted_tags = sorted(
        counts.keys(),
        key=lambda tag: (families.get(tag, ""), tag),
    )

    for tag in sorted_tags:
        family = families.get(tag, "unknown")
        print("{:<{width}}  {:5d}  {}".format(
            tag,
            counts[tag],
            family,
            width=width,
        ))


def copy_problem(records, target, dest):
    """Copy a problem and all of its recorded include files to dest."""
    target_lower = target.lower()
    matches = []

    for record in records:
        rel_path = record.get("rel_main_file", "").lower()
        full_path = record.get("main_file", "").lower()

        if target_lower in rel_path or target_lower in full_path:
            matches.append(record)

    if not matches:
        print("No indexed problem matches: {}".format(target))
        return

    if len(matches) > 1:
        print("Multiple matches for {}; use a more specific path.".format(target))
        for record in matches[:25]:
            print("  {}".format(record.get("rel_main_file", "")))
        return

    record = matches[0]
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)

    files_to_copy = set()
    files_to_copy.add(Path(record["main_file"]))

    for inc in record.get("includes", []):
        if inc.get("status") == "ok":
            files_to_copy.add(Path(inc["include_file"]))

    for source_file in sorted(files_to_copy):
        target_path = dest / source_file.name
        shutil.copy2(source_file, target_path)
        print("copied {} -> {}".format(source_file, target_path))


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def build_arg_parser():
    parser = argparse.ArgumentParser(
        prog="w3dexamples.py",
        description="Index and search WARP3D verification examples.",
    )

    sub = parser.add_subparsers(dest="command")

    p_index = sub.add_parser("index", help="Build verification problem index.")
    p_index.add_argument("root", type=Path, help="Verification directory root.")
    p_index.add_argument("--terms", type=Path, required=True,
                         help="Plain-text WARP3D feature term file.")
    p_index.add_argument("--index", type=Path, default=Path(DEFAULT_INDEX_FILE))
    p_index.add_argument("--catalog", type=Path, default=Path(DEFAULT_CATALOG_FILE))
    p_index.add_argument("--no-catalog", action="store_true")

    p_find = sub.add_parser("find", help="Search indexed verification problems.")
    p_find.add_argument("query", nargs="?", default="", help="Keyword or phrase query.")
    p_find.add_argument("--index", type=Path, default=Path(DEFAULT_INDEX_FILE))
    p_find.add_argument("--feature", action="append", default=[],
                        help="Required feature tag.")
    p_find.add_argument("--limit", type=int, default=20)
    p_find.add_argument("--debug", action="store_true",
                        help="Show scores, include counts, and evidence lines.")

    p_show = sub.add_parser("show", help="Show one indexed problem card.")
    p_show.add_argument("target", help="Path fragment or problem file.")
    p_show.add_argument("--index", type=Path, default=Path(DEFAULT_INDEX_FILE))

    p_tags = sub.add_parser("tags", help="List detected feature tags.")
    p_tags.add_argument("--index", type=Path, default=Path(DEFAULT_INDEX_FILE))

    p_copy = sub.add_parser("copy", help="Copy a problem and its included files.")
    p_copy.add_argument("target", help="Path fragment or problem file.")
    p_copy.add_argument("dest", type=Path, help="Destination directory.")
    p_copy.add_argument("--index", type=Path, default=Path(DEFAULT_INDEX_FILE))

    return parser


def load_records_for_command(index_file):
    """Read records for commands that operate on an existing index."""
    if not index_file.exists():
        print("error: index file not found: {}".format(index_file), file=sys.stderr)
        print("run: w3dexamples.py index /path/to/verification", file=sys.stderr)
        return None

    root, records = read_json_index(index_file)
    return root, records


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 2

    if args.command == "index":
        root = args.root

        if not root.exists() or not root.is_dir():
            print("error: verification root is not a directory: {}".format(root),
                  file=sys.stderr)
            return 2

        try:
            feature_terms = load_feature_terms(args.terms)
        except RuntimeError as exc:
            print("error reading feature terms: {}".format(exc), file=sys.stderr)
            return 2

        records = build_index(root, feature_terms)
        write_json_index(records, args.index, root, feature_terms)

        print("loaded {} feature term(s)".format(len(feature_terms)))
        print("indexed {} verification problem(s)".format(len(records)))
        print("wrote {}".format(args.index))

        if not args.no_catalog:
            write_markdown_catalog(records, args.catalog)
            print("wrote {}".format(args.catalog))

        return 0

    loaded = load_records_for_command(args.index)
    if loaded is None:
        return 2

    _root, records = loaded

    if args.command == "find":
        print_find_results(
            records,
            query=args.query,
            required_features=args.feature,
            limit=args.limit,
            debug=args.debug,
        )
        return 0

    if args.command == "show":
        show_record(records, args.target)
        return 0

    if args.command == "tags":
        print_tags(records)
        return 0

    if args.command == "copy":
        copy_problem(records, args.target, args.dest)
        return 0

    parser.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
