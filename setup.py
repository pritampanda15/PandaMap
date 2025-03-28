from setuptools import setup, find_packages

# This is a shim for compatibility with older tools
# pyproject.toml is the preferred way to configure package metadata
setup(
    package_dir={"": "src"},
    packages=find_packages(where="src"),
)
