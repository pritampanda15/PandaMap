"""
Tell users when a newer PandaMap is available.

The check runs in a background thread, asks PyPI at most once a day (the
answer is cached under ``~/.cache/pandamap``), never raises, and prints one
notice per process to stderr.  Set ``PANDAMAP_NO_UPDATE_CHECK=1`` to turn it
off.
"""
import json
import os
import re
import sys
import threading
import time
from importlib.metadata import version

ENV_DISABLE = "PANDAMAP_NO_UPDATE_CHECK"
CHECK_INTERVAL = 24 * 3600           # seconds between questions to PyPI
_CACHE = os.path.join(os.path.expanduser("~"), ".cache", "pandamap", "update_check.json")
_lock = threading.Lock()
_notified = False


def _version_tuple(text):
    """Leading numeric components of a version string: '4.3.0rc1' -> (4, 3, 0)."""
    return tuple(int(n) for n in re.findall(r"\d+", text.split("+")[0])[:4])


def update_message(latest, current, package_name="pandamap"):
    return (
        f"\nPandaMap {latest} is available (you have {current}).\n"
        f"Please use the latest version for improved utilization of PandaMap:\n"
        f"    pip install --upgrade {package_name}\n"
        f"(Set {ENV_DISABLE}=1 to hide this notice.)\n"
    )


def _latest_version(package_name):
    """Latest version on PyPI, from the daily cache when it is fresh."""
    try:
        with open(_CACHE, encoding="utf-8") as fh:
            cached = json.load(fh)
        if time.time() - float(cached["checked"]) < CHECK_INTERVAL:
            return cached["latest"]
    except Exception:
        pass
    import requests
    response = requests.get(f"https://pypi.org/pypi/{package_name}/json", timeout=2)
    latest = response.json()["info"]["version"]
    try:
        os.makedirs(os.path.dirname(_CACHE), exist_ok=True)
        with open(_CACHE, "w", encoding="utf-8") as fh:
            json.dump({"checked": time.time(), "latest": latest}, fh)
    except OSError:
        pass
    return latest


def check_for_updates(package_name="pandamap", env_var_disable=ENV_DISABLE):
    """
    Print a notice to stderr when PyPI has a newer version than the one
    installed.  Returns the newer version string, or None.

    Never raises: any failure (offline, proxy, missing metadata) is ignored.
    """
    global _notified
    if os.getenv(env_var_disable):
        return None
    try:
        current = version(package_name)
        latest = _latest_version(package_name)
        if _version_tuple(latest) > _version_tuple(current):
            with _lock:
                if not _notified:
                    _notified = True
                    print(update_message(latest, current, package_name), file=sys.stderr)
            return latest
    except Exception:
        pass  # Silently ignore errors
    return None


def check_for_updates_in_background():
    """Start the check without delaying the caller."""
    if os.getenv(ENV_DISABLE):
        return None
    thread = threading.Thread(target=check_for_updates, daemon=True)
    thread.start()
    return thread
