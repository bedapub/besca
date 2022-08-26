""" setup procedure """

if __name__ == '__main__':

    from setuptools import setup, find_packages
    import versioneer

    reqs = open('requirements.txt', encoding='utf-8').readlines()
    requires = [req.strip() for req in reqs]

    with open("README.md", "r") as fh:
      long_description = fh.read()

    setup(name='besca',
          version=versioneer.get_version(),
          cmdclass=versioneer.get_cmdclass(),
          description='Collection of BEDA internal python functions for analysing single cell RNAseq data.',
          long_description=long_description,
          long_description_content_type="text/markdown",
          classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Developers',
            'Topic :: Software Development :: Build Tools',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.8'],
          url='https://github.com/bedapub/besca',
          license='GPLv3',
          author='BEDA community',
          author_email='manuel.kohler@roche.com',
          packages=find_packages(exclude=["devtools", "tests"]),
          zip_safe=False,
          package_data={'besca.datasets.data': ['*.h5ad'],
                        'besca.st': ['*.css'],
                        'besca.datasets.nomenclature': ['*.tsv'],
                        'besca.datasets.mito_files': ['*.tsv'],
                        'besca.datasets.genesets': ['*.gmt', '*.tsv']},
          install_requires=requires)
