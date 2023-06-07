import setuptools
from glob import glob

version = [l.strip() for l in open("dengue_ngs/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(
	name="dengue_wgs",
	version=version,
	packages=["dengue_wgs"],
	license="GPLv3",
	long_description="Dengue wgs command line tool",
	scripts= glob("scripts/*"),
)