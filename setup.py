from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["pip>=19.0.3", "wheel>=0.33.1", "numpy>=1.12.0", "nibabel>=2.1.0", "SciPy>=1.4.0"]]

setup(
    name="cluster_roi",
    version="1.1.0",
    author="Craddock, R C and James, G A and Holtzheimer, P E and Hu, X P and Mayberg, H S",
    author_email="cameron.craddock@gmail.com",
    description="A whole brain fMRI atlas generated via spatially constrained spectral clustering",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/ccraddock/cluster_roi/",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"],
)
