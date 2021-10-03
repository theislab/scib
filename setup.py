from setuptools import setup

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()
    requirements = [x for x in requirements if not x.startswith("#") and x != ""]
    
with open("requirements_extra.txt", "r") as f:
    requirements_extra = f.read().splitlines()
    requirements_extra = [x for x in requirements_extra if not x.startswith("#") and x != ""]


setup(name='scib',
      version='0.1.1',
      description='Benchmark tools for single cell data integration',
      author='Malte Luecken, Maren Buettner, Daniel Strobl, Michaela Mueller',
      author_email='malte.luecken@helmholtz-muenchen.de',
      packages=['scib', 'scib.metrics'],
      package_data={'scib': ['resources/*.txt', 'knn_graph/*']},
      zip_safe=False,
      license='MIT',
      url='https://github.com/theislab/scib',
      keywords = ['benchmark', 'single cell', 'data integration'],
      install_requires=requirements,
      extras_require={'integration': requirements_extra},
      classifiers=[
         'Development Status :: 3 - Alpha',      
         'Intended Audience :: Developers',      
         'Topic :: Software Development :: Build Tools',
         'License :: OSI Approved :: MIT License',   
         'Programming Language :: Python :: 3',      
         'Programming Language :: Python :: 3.7',
      ])
