import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="besca", # Replace with your own username
    version="2.1",
    description="collection of BEDA internal python functions for analysing single cell RNAseq data",
    long_description= 'please view https://pages.github.roche.com/BEDA/besca/',
    author='BEDA community',
    license='GPLv3',
    author_email='alice.julien-laferriere@roche.com',
    url='https://github.com/bedapub/besca',
    packages=setuptools.find_packages(),
    package_data={'besca.datasets.data':['*.h5ad'], 'besca.st':['*.css'],\
        'besca.datasets.mito_files': ['*.tsv'], 'besca.export':['reformat'],\
        'besca.datasets.genesets':['HumanCD45p_scseqCMs6.gmt','Immune.txt']},
    classifiers=[
        'Development Status :: Release 2.1 (Scanpy 1.4.1)',
        'Programming Language :: Python :: 3.7.1'
        'License :: GPLv3'
    ]
    )
