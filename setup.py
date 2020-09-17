import os
import sys
import shutil
from subprocess import call
from setuptools import setup
from warnings import warn
import py_compile
from pathlib import Path


# Replace py_compile.compile with a function that calls it with doraise=True
# so stop when there is a syntax error
orig_py_compile = py_compile.compile


def doraise_py_compile(file, cfile=None, dfile=None, doraise=False):
    orig_py_compile(file, cfile=cfile, dfile=dfile, doraise=True)


py_compile.compile = doraise_py_compile

if sys.version_info.major != 3:
    raise RuntimeError("SEQC requires Python 3")
if sys.version_info.minor < 5:
    warn("Multiprocessing analysis methods may not function on Python versions < 3.5")

main_ns = {}

# get version
with open("src/seqc/version.py") as f:
    exec(f.read(), main_ns)

setup(
    name="seqc",
    version=main_ns["__version__"],
    description="Single Cell Sequencing Processing and QC Suite",
    author="Ambrose J. Carr",
    author_email="mail@ambrosejcarr.com",
    package_dir={"": "src"},
    package_data={"": ["*.r", "*.R"]},
    packages=[
        "seqc",
        "seqc.sequence",
        "seqc.alignment",
        "seqc.core",
        "seqc.stats",
        "seqc.summary",
        "seqc.notebooks",
    ],
    install_requires=[
        dep.strip() for dep in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    scripts=["src/scripts/SEQC"],
    extras_require={"GSEA_XML": ["html5lib", "lxml", "BeautifulSoup4"]},
    include_package_data=True,
)

# look for star
if not shutil.which("STAR"):
    warn("SEQC: STAR is not installed. SEQC will not be able to align files.")

# get location of setup.py
setup_dir = os.path.dirname(os.path.realpath(__file__))
seqc_dir = os.path.expanduser("~/.seqc/seqc")

print("setup_dir: {}".format(setup_dir))
print("seqc_dir: {}".format(seqc_dir))

if os.path.isdir(seqc_dir):
    shutil.rmtree(seqc_dir)


def ignore_test_and_tools(dir_, files):
    """Filter files to be moved by shutil.copytree. Ignore any hidden file and the
    test and tools directories, which are not needed by the remote instance.
    :param dir_: dummy variable, must be present to be passed to shutil.copytree()
    :param files: output of os.listdir(), files to be subjected to filtering
    :return list: list of files that should be filtered, and not copied.
    """
    return [f for f in files if (f == "test" or f.startswith("."))]


# install tools and a local copy of seqc.
# copy seqc repository
shutil.copytree(setup_dir, seqc_dir, ignore=ignore_test_and_tools)
shutil.make_archive(base_name=seqc_dir, format="gztar", root_dir=seqc_dir)
