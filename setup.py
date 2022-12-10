from setuptools import setup, find_packages
import os

with open('README.md', 'r', encoding='utf-8') as fh:
	long_description = fh.read()

setup(
	name = 'hairpin',
	version='2.0.0',
	author='ampersandor',
	author_email='ddong3525@gmail.com',
	description='hairpin calculator',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/ampersandor/hairpin',
	classifiers=[
		'Programming Language :: Python :: 3.8',
		'Operating System :: OS Independent',
	],
	package_dir={'': '.'},
	packages=find_packages(where='.'),
	python_requires='>=3.8'
)

