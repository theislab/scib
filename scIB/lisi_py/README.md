# Compile `compute_simpson_index_graph_lisi.pyx` 

We use a cythonized version of `graph_lisi` for performance reasons. Here, we describe briefly how to create to `.so` file, 
which can be imported as module in `Python 3` (Details can be found in the
[Cython NumPy Tutorial](https://cython.readthedocs.io/en/latest/src/userguide/numpy_tutorial.html)).

First, make sure that you have installed `cython`. Then, run:
```
$ cython -3 compute_simpson_index_graph_cy.pyx
```
The option `-3` enables to use `Python3` language. 

Compile the resulting `.c` file (**important**: we did not use `cimport numpy` in the header, if you want to, specify where 
the numpy c header files are located).
```
$ gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python3.7m -o compute_simpson_index_graph_cy.so compute_simpson_index_graph_cy.c
```
The `.so` file can be imported as module in `Python 3`:
```
from compute_simpson_index_graph_lisi_cy import compute_simpson_index_graph_lisi_cy
```
