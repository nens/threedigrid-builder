import os
import sys

# Prefer the in-repo package over any installed version in site-packages.
# This ensures the compiled extension in the workspace (fgrid .so) is used.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)
