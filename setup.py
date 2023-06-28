import setuptools
from glob import glob

version = [l.strip() for l in open("dengue_ngs/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(
	name="dengue_ngs",
	version=version,
	packages=["dengue_ngs"],
	license="GPLv3",
	long_description="Dengue wgs command line tool",
	scripts= glob("scripts/*"),
    data_files=[('share/dengue-ngs', ['data/sample_exclusion.txt'])],
)