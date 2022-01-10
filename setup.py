from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'GEDSpy'
LONG_DESCRIPTION = 'GEDSpy is the python library for gene list enrichment with genes ontology, pathways and potential drugs'

# Setting up
setup(
        name="GEDSpy", 
        version=VERSION,
        author="Jakub Kubis",
        author_email="jbiosystem@gmail.com",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages("GEDSpy"),
        install_requires=['requests', 'bioservices', 'pandas', 'tqdm', 'seaborn', 'matplotlib', 'scipy', 'networkx', 'pyvis'],       
        keywords=['python', 'GO', 'pathways', 'drug', 'gene ontology'],
        license = 'MIT',
        classifiers = [
            "Development Status :: 3 - Alpha",
            "Programming Language :: Python :: 3",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: Linux :: Ubuntu",
        ],
        python_requires='>=3.6',
)


