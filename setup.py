import os
import shutil

from setuptools import setup

if shutil.which('g++') is not None:
    os.system('g++ -std=c++11 -O3 scib/knn_graph/knn_graph.cpp -o scib/knn_graph/knn_graph.o')
setup()
