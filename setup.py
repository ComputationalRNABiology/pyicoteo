from __future__ import with_statement
from distutils.core import setup
from pyicoslib.parser import VERSION


setup(name='Pyicos',
      version=VERSION,
      description='Mapped reads analysis tool and library',
      author=u'Juan Gonzalez_Vallinas',
      author_email='juanramon.gonzalezvallinas@upf.edu',
      url='http://regulatorygenomics.upf.edu/pyicos',
      packages = ['pyicoslib.lib', 'pyicoslib.chromlen'],
      package_data={
          'pyicoslib.chromlen': [
              '*'
          ],
      },
      scripts = ['pyicos'],
      py_modules = ['pyicoslib.core', 'pyicoslib.turbomix','pyicoslib.utils', 'pyicoslib.parser', 'pyicoslib.defaults', 'pyicoslib.bam', 'pyicoslib.enrichment']
     )

