import warnings
import requests
import shutil
from packaging.version import Version
from importlib.metadata import version


def check_for_updates(package_name="pandamap", env_var_disable="PANDAMAP_NO_UPDATE_CHECK"):
    """
    Check PyPI for a newer version and notify the user.

    Parameters:
    - package_name (str): The name of the package to check.
    - env_var_disable (str): Environment variable to skip version check.
    """
    import os

    if os.getenv(env_var_disable):
        return

    try:
        current_version = version(package_name)
        response = requests.get(
            f"https://pypi.org/pypi/{package_name}/json", timeout=2
        )
        latest_version = response.json()["info"]["version"]

        if Version(latest_version) > Version(current_version):
            message = (
                f"\n\U0001F6A8 [bold red]{package_name} {latest_version} is available![/bold red] "
                f"[dim](you have {current_version})[/dim]\n\n"
                f"[yellow]Update with:[/yellow] [green]pip install --upgrade {package_name}[/green]\n"
                f"[dim]To disable update checks, set: {env_var_disable}=1[/dim]\n"
            )

            if shutil.which("rich"):
                try:
                    from rich.console import Console
                    console = Console()
                    console.print(message)
                except ImportError:
                    warnings.warn(message, UserWarning, stacklevel=2)
            else:
                warnings.warn(message, UserWarning, stacklevel=2)

    except Exception:
        pass  # Silently ignore errors
