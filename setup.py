from setuptools import setup

setup(name='scIB',
      version='0.1',
      description='Benchmark tools for single cell data integration',
      url='https://github.com/theislab/Benchmarking_data_integration',
      author='Daniel Strobl, Michaela Müller',
      author_email='',
      packages=['scIB', 'scanpy'],
      package_data={'scIB': ['resources/*.txt']},
      zip_safe=False)

