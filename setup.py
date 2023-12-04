from setuptools import setup
#from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

#with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    #long_description = f.read()

setup(
  name = 'TARDIS: Targets Discoverer',
  packages = ['TARDIS'],
  version = '1.0',
  description = 'A pipeline that uses genome-scale metabolic network modeling to find new targets for antimicrobial drug-design.',
  #long_description = long_description,
  license='CC-BY-4.0',
  author = 'Camila Rodrigues Chaves, Pedro Torres',
  author_email = 'camilachaves@biof.ufrj.br',
  #url = 'https://github.com/monteirotorres/ProtCHOIR',
  #include_package_data = True,
  #package_data={'ProtCHOIR':['Contents/*.html', 'Contents/*.svg']},
  keywords = ['GEM',
              'targets',
              'antimicrobial',
              'drug-design',
              'network model',
              'metabolic network'],
  classifiers = [],
  py_modules=[
      'TARDIS.__main__',
      'TARDIS.FindTargets',
      'TARDIS.Homology-search',
      'TARDIS.Initialise'],
  entry_points={
      'console_scripts': [
          'TARDIS = TARDIS.__main__:main']}
)

  #install_requires=['cobrapy','contrabass','subprocess','biopython','pandas','carveme','ncbi-blast'])
