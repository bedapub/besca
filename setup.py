""" setup procedure """

if __name__ == '__main__':

    from setuptools import setup, find_packages
    import versioneer

    reqs = open('requirements.txt', encoding='utf-8').readlines()
    requires = [req.strip() for req in reqs]

    setup(name='besca',
          version=versioneer.get_version(),
          cmdclass=versioneer.get_cmdclass(),
          description='collection of BEDA internal python functions for'
                      'analysing single cell RNAseq data',
          long_description='View https://pages.github.com/bedapub/besca/',
          classifiers=[
              'Development Status :: Release 2.5 (Scanpy 1.8)',
              'Programming Language :: Python :: 3.7.1',
              'License :: GPLv3'],
          url='https://github.com/bedapub/besca',
          license='GPLv3',
          author='BEDA community',
          author_email='alice.julien-laferriere@roche.com',
          packages=find_packages(exclude=["devtools", "tests"]),
          zip_safe=False,
          package_data={'besca.datasets.data': ['*.h5ad'],
                        'besca.st': ['*.css'],
                        'besca.datasets.nomenclature': ['*.tsv'],
                        'besca.datasets.mito_files': ['*.tsv'],
                        'besca.datasets.genesets': ['*.gmt', '*.tsv']},
          install_requires=requires)
