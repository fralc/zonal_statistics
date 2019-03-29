import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="zonal-statistics-flc",
    version="0.0.1",
    author="Francesco Lo Conti",
    author_email="francesco.loconti@leitha.eu",
    description="Zonal statistics tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://gitlab-leitha.servizi.gr-u.it/LEI00025/zonal_statistics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
