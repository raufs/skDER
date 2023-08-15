from setuptools import setup
import os

setup(name='skDER',
      version='1.0.4',
      description='Program to select distinct representatives from an input set of microbial genomes.',
      url='http://github.com/raufs/skDER/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['skDER'],
      scripts=['bin/skder'],
      zip_safe=False)

try:
    os.system("g++ -o skDERcore skDERcore.cpp")
    os.system("g++ -o skDERsum skDERsum.cpp")
except:
    print('C++ compilation failed')
    pass
