import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FatropOCPSpecification",
    version="0.0.2",
    author="Lander Vanroye",
    author_email="lander.vanroye@kuleuven.be",
    description="package for fatrop problem specification",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)