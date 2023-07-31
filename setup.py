import os
import sys
from setuptools import setup

setup(name='skDER',
      version='1.0.1',
      description='Dynamic dereplication of microbial genomes.',
      url='http://github.com/raufs/skDER/',
      author='Rauf Salamzade',
      author_email='salamzader@gmail.com',
      license='BSD-3',
      packages=['skDER'],
      scripts=['skDER.py'],
      zip_safe=False)

try:
    os.system("g++ -o skDERcore skDERcore.cpp")
except:
    print('C++ compilation failed')
    pass
