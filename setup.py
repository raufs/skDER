#!/usr/bin/env python
from setuptools import setup
import os

if __name__ == "__main__":
    setup()
    try:
        os.system("g++ -o skDERsum skDERsum.cpp")
        os.system("g++ -o skDERcore skDERcore.cpp")
        os.system("mv skDERsum $CONDA_PREFIX/bin/")
        os.system("mv skDERcore $CONDA_PREFIX/bin/")
    except:
        print('C++ compilation or setup of executables failed')
