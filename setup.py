#!/usr/bin/env python
from setuptools import setup
import os

if __name__ == "__main__":
    setup()
    try:
        os.system("g++ -o skDERsum skDERsum.cpp")
        os.system("g++ -o skDERcore skDERcore.cpp")
    except:
        print('C++ compilation failed')
