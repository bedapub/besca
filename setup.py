import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='besca',
        version='2.1',
        description='collection of BEDA internal python functions for analysing single cell RNAseq data',
        long_description= 'please view https://pages.github.com/bedapub/besca/',
        classifiers=[
          'Development Status :: Release 2.1 (Scanpy 1.4.1)',
          'Programming Language :: Python :: 3.7.1',
          'License :: GPLv3'],
        url='https://github.com/bedapub/besca',
        license='GPLv3',
        author='BEDA community',
        author_email='alice.julien-laferriere@roche.com',
        packages=find_packages(),
        zip_safe=False,
        package_data={'besca.datasets.data':['*.h5ad'], 'besca.st':['*.css'],\
            'besca.datasets.mito_files': ['*.tsv'], 'besca.export':['reformat'],\
            'besca.datasets.genesets':['HumanCD45p_scseqCMs6.gmt','Immune.txt']},
        install_requires=requires)      
