import os

from setuptools import setup, Extension
#from distutils.core import setup
#from Cython.Build import cythonize


HERE = os.path.abspath(os.path.dirname(__file__))
about = {}
version = os.path.join(HERE, "hmseekr", "__version__.py")
with open(version, "r", encoding="utf-8") as f:
    exec(f.read(), about)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Check if the Cython-generated C file exists
c_source = "./hmseekr/kmersc.c"
if os.path.exists(c_source):
    ext_modules = [Extension("kmersc", sources=[c_source])]
else:
    raise FileNotFoundError(f"Required source file not found: {c_source}")


requirements = [
    "tqdm",
    "numpy",
    "pandas",
    "scipy",
]

# ensure matplotlib to be greater than 3.5.3 which is compatible with python 3.7 and above
# and could accomodate the fontmanager function

test_requirements = ["pytest"]

setup(
    name=about["__title__"],
    python_requires=">=3.9",
    version=about["__version__"],
    install_requires=requirements,
    tests_require=test_requirements,
    description=about["__description__"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=about["__url__"],
    author=about["__author__"],
    author_email=about["__author_email__"],
    license=about["__license__"],
    packages=["hmseekr"],
    zip_safe=False,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.9",
    ],
    ext_modules = ext_modules,
    entry_points={
        "console_scripts": [
            "hmseekr_kmers = hmseekr.console_scripts:console_hmseekr_kmers",
            "hmseekr_train = hmseekr.console_scripts:console_hmseekr_train",
            "hmseekr_findhits = hmseekr.console_scripts:console_hmseekr_findhits",
            "hmseekr_gridsearch = hmseekr.console_scripts:console_hmseekr_gridsearch",
            "hmseekr = hmseekr.console_scripts:console_hmseekr_help"
        ]
    },
)
