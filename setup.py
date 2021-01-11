import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gtfSymbiosis",
    version="0.0.1",
    author="Samuel Rosean",
    author_email="samrosean@gmail.com",
    packages=['primaryfunctions'],
    description="construction of transcriptome from GTF files, and comparison between GTf files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/samrosean/GTF_Symbiosis",
    ##packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
