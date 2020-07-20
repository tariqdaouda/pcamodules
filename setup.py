from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

setup(
    name='pcamodules',

    version='0.1',

    description='',
    long_description="""""",
    url='',

    author='Tariq Daouda',
    author_email='',

    license='MIT',
    
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    python_requires='>=3.5',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],

    install_requires=["numpy", "sklearn"],

    keywords='',

    packages=find_packages(),

    entry_points={}
)
