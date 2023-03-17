from setuptools import find_packages
from setuptools import setup


INSTALL_REQUIRES = [
    "scanpy==1.9.3",
    "pyyaml==6.0",
    "doubletdetection==4.2",
    "harmonypy==0.0.9",
    "openpyxl==3.1.1",
    "anticor-features==0.2.0",
    "numpy==1.23",
    "gin-config==0.5.0",
]

TESTS_REQUIRE = ["pytest"]

with open("README.md", "r", encoding="utf-8") as file:
    long_description = file.read()

setup(
    name="moresca",
    version="0.1.0",
    description="Reproducible boilerplate-free workflow management for Scanpy-based scRNA-seq analysis.",
    license="AGPL-3.0 license",
    packages=find_packages(),
    author="Matthias Bruhns",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author_email="matthias.bruhns@posteo.de",
    url="https://github.com/claassenlab/MORESCA",
    keywords=["scRNA-seq", "reproducibility", "workflow management"],
    install_requires=INSTALL_REQUIRES,
    test_suite="tests",
    tests_require=TESTS_REQUIRE,
    extras_require={"dev": ["black", "flake8"]},
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
