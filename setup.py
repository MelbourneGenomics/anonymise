#!/usr/bin/env python

from distutils.core import setup

setup(
    name='anonymise',
    version='0.0.1',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['anonymise'],
    package_dir={'anonymise': 'anonymise'},
    package_data={'anonymise': ['data/application_json_schema.txt']},
    entry_points={
        'console_scripts': ['anonymise = anonymise.anon:main']
    },
    url='https://github.com/bjpop/anonymise',
    license='LICENSE.txt',
    description=( 'Anonymise data'),
    long_description=('Anonymise data'),
    install_requires=[
        "jsonschema == 2.5.1",
        "functools32"
    ],
)
