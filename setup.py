from distutils.core import setup

setup(name='Pyicos',
      version='0.8.5',
      description='Mapped reads analysis tool and library',
      author=u'Juan Gonzalez_Vallinas',
      author_email='juanramon.gonzalezvallinas@upf.edu',
      url='http://regulatorygenomics.upf.edu/pyicos',
      packages = ['pyicoslib.lib', 'chrdesc'],
      package_data = {'chrdesc' : ['mm8', 'mm9', 'hg18', 'hg19'] },
      scripts = ['pyicos'],
      py_modules = ['pyicoslib.core', 'pyicoslib.operations', 'pyicoslib.parser', 'pyicoslib.defaults']
     )

