"""Test configuration: never contact PyPI from the test suite."""
import os

os.environ.setdefault("PANDAMAP_NO_UPDATE_CHECK", "1")
