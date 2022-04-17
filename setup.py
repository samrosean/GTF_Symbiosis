from setuptools import setup, find_packages

setup(name = 'gtfSymbiosis', version = '0.0.2', author = 'Samuel Rosean', author_email = 'samrosean@gmail.com', description = 'construction of transcriptome from GTF files, and comparison between GTF files', url = 'https://github.com/samrosean/GTF_Symbiosis', packages=find_packages(','),install_requires=['pandas', 'numpy', 'matplotlib', 'matplotlib_venn', 'prettytable', 'prettytable', 'networkx', 'venn', 'upsetplot', 'seaborn', 'pyranges'] )


