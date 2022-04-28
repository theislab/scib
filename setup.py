import os

from setuptools import setup

os.system('g++ -std=c++11 -O3 scib/knn_graph/knn_graph.cpp -o scib/knn_graph/knn_graph.o')
setup()
