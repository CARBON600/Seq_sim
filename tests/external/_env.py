"""Directory lookups shared by the external check scripts.

Nothing here falls back to a developer's own layout: a script that needs a
directory either finds the variable set or stops with a message naming it.
"""
import os
import sys


def directory(name, hint):
    """Return os.environ[name], or exit with a message saying what to set."""
    value = os.environ.get(name)
    if not value:
        sys.exit("set %s to %s" % (name, hint))
    return value
