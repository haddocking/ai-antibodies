"""Setup aiabs."""
from setuptools import find_packages, setup

with open("requirements.txt") as f:
    required = f.read().splitlines()

setup(
    name="aiabs",
    license="Apache License 2.0",
    version="0.0.0",
    author="",
    description="",
    author_email="",
    include_package_data=True,
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[],
    python_requires=">=3.10, <4",
    install_requires=required,
    entry_points={
        "console_scripts": [
            "aiabs=aiabs.cli:maincli",
            "aiabs-analysis=aiabs.cli_analysis:maincli"
        ],
    },
)