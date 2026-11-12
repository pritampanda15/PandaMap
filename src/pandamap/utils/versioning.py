import os
import re
import sys
from importlib.metadata import version


def _version_tuple(text):
    """Leading numeric components of a version string: '4.3.0rc1' -> (4, 3, 0)."""
    return tuple(int(n) for n in re.findall(r"\d+", text.split("+")[0])[:4])


def check_for_updates(package_name="pandamap", env_var_disable="PANDAMAP_NO_UPDATE_CHECK"):
    """
    Check PyPI for a newer version and print a one-line notice to stderr.

    Never raises: any failure (offline, proxy, missing metadata) is ignored.
    Set the environment variable ``PANDAMAP_NO_UPDATE_CHECK=1`` to skip it.
    """
    if os.getenv(env_var_disable):
        return

    try:
        import requests

        current_version = version(package_name)
        response = requests.get(
            f"https://pypi.org/pypi/{package_name}/json", timeout=2
        )
        latest_version = response.json()["info"]["version"]

        if _version_tuple(latest_version) > _version_tuple(current_version):
            print(
                f"{package_name} {latest_version} is available "
                f"(you have {current_version}). "
                f"Update with: pip install --upgrade {package_name}  "
                f"[set {env_var_disable}=1 to silence this notice]",
                file=sys.stderr,
            )
    except Exception:
        pass  # Silently ignore errors
