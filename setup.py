from setuptools import setup

setup(name='scIB',
      version='0.1',
      description='Benchmark tools for single cell data integration',
      url='https://github.com/theislab/Benchmarking_data_integration',
      author='Malte Luecken, Maren Buettner, Daniel Strobl, Michaela Mueller',
      author_email='malte.luecken@helmholtz-muenchen.de',
      packages=['scIB'],
      package_data={'scIB': ['resources/*.txt']},
      zip_safe=False)

