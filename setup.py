from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['kmeansmodule.c'])
setup(name='mykmeanssp',
     version='1.0',
     description='Python wrapper for K-Means++ algorithm',
     ext_modules=[module])