import subprocess
import sys

from setuptools import setup

try:
    cmd = (
        "g++ -std=c++11 -O3 scib/knn_graph/knn_graph.cpp -o scib/knn_graph/knn_graph.o"
    )
    sys.stdout.write("Compile knn_graph C++ code for LISI metric...\n")
    sys.stdout.flush()
    subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, text=True)
except subprocess.CalledProcessError as exc:
    sys.stdout.write(
        f"Failed to compile knn_graph for LISI - skipping...\n{exc.returncode}\n{exc.output}"
    )
    sys.stdout.flush()

setup()
