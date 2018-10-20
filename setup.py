import pathotyping

from setuptools import setup

VERSION = pathotyping.__version__

with open('README.md') as fh:
    README = fh.read()

setup(
    name='patho_typing',
    version='{}'.format(VERSION),
    packages=['pathotyping',
              'pathotyping.modules'],
    package_dir={'pathotyping': 'pathotyping'},
    package_data={'pathotyping': ['../.git/*', '../.git/*/*', '../.git/*/*/*',
                                  'modules/seq_conf/*/*']},
    include_package_data=True,
    data_files=[('', ['LICENSE'])],
    install_requires=[
        'ReMatCh'
    ],
    description='In silico pathogenic typing directly from raw Illumina reads',
    long_description=README,
    long_description_content_type='text/markdown',
    keywords=['reference mapping', 'pathotyping', 'sequence presence/absence', 'Escherichia coli',
              'Haemophilus influenzae', 'Yersinia enterocolitica'],
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Environment :: Console',
        'Operating System :: Unix',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    url='https://github.com/B-UMMI/patho_typing',
    author='Miguel P. Machado',
    author_email='mpmachado@medicina.ulisboa.pt',
    license='GPL3',
    # To use entry_points with .py the first folder cannot have the same name of the script
    entry_points={
        'console_scripts': [
            'patho_typing.py = pathotyping.patho_typing:main'
        ]
    },
    python_requires='>=3.4'
)
