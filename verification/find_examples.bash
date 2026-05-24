#!/usr/bin/env bash
#
# find_examples
#
# Convenience wrapper for:
#
#     python3 w3dexamples.py find "search terms"
#
# Usage:
#
#     find_examples "mises q3disop domain"
#     find_examples --debug "tie-mesh bilinear"
#     find_examples --feature material_mises --feature element_q3disop
#
# Notes:
#   - This script finds a usable Python 3.8+ command automatically.
#   - It assumes w3dexamples.py is in the same directory as this script.
#   - All command-line arguments are passed directly to:
#
#         w3dexamples.py find
#

set -e

# Directory containing this wrapper script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

W3DEXAMPLES="${SCRIPT_DIR}/w3dexamples.py"

if [[ ! -f "$W3DEXAMPLES" ]]; then
    echo "Error: cannot find w3dexamples.py in:" >&2
    echo "  $SCRIPT_DIR" >&2
    exit 1
fi

# Return success only if the supplied Python command is version 3.8 or newer.
check_python_38()
{
    "$@" -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 8) else 1)'         >/dev/null 2>&1
}

# Print the version of the supplied Python command.
python_version()
{
    "$@" -c 'import sys; print("{}.{}.{}".format(*sys.version_info[:3]))'         2>/dev/null || true
}

# Find a usable Python 3.8+ command.
PY=()

if command -v python3 >/dev/null 2>&1 && check_python_38 python3; then
    PY=(python3)
elif command -v python >/dev/null 2>&1 && check_python_38 python; then
    PY=(python)
elif command -v py >/dev/null 2>&1 && check_python_38 py -3; then
    # Windows Python launcher, useful in Git Bash/MSYS contexts.
    PY=(py -3)
fi

if [[ ${#PY[@]} -eq 0 ]]; then
    echo "Error: could not find Python 3.8 or newer." >&2
    echo "Tried: python3, python, py -3" >&2
    echo "" >&2

    if command -v python3 >/dev/null 2>&1; then
        echo "python3 version: $(python_version python3)" >&2
    fi
    if command -v python >/dev/null 2>&1; then
        echo "python version:  $(python_version python)" >&2
    fi
    if command -v py >/dev/null 2>&1; then
        echo "py -3 version:   $(python_version py -3)" >&2
    fi

    exit 1
fi
echo " "
exec "${PY[@]}" "$W3DEXAMPLES" find "$@"
