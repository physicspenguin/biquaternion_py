#!/usr/bin/env python3

from setuptools import setup

me = "Daren A. Thimm"
my_email = "daren.thimm@uibk.ac.at"
setup(
    name="biquaternion-py",
    version="0.1.1",
    packages=[
        "biquat_py",
    ],
    license="GPL v3",
    license_files="LICENSE",
    url="https://git.uibk.ac.at/geometrie-vermessung/biquaternion_py",
    author=me,
    author_email=my_email,
    maintainer=me,
    maintainer_email=my_email,
    description="Implementation of general biquaternion algebras.",
    long_description=open("./README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "numpy==1.24.3",
        "setuptools==67.7.0",
        "sympy==1.12",
    ],
    extras_require={
        "dev": [
            "flake8 == 6.0.0",
            "pytest == 7.2.0",
            "pre-commit == 3.3.3",
            "pytest-cov == 4.0.0",
        ]
    },
)
