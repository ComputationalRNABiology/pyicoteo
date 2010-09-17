from distutils.core import setup

setup(name='Pyicos',
      version='0.8.5',
      description='Mapped reads analysis tool and library',
      author=u'Juan Gonzalez_Vallinas',
      author_email='juanramon.gonzalezvallinas@upf.edu',
      url='http://regulatorygenomics.upf.edu/pyicos',
      package_data = {'chrdesc' : ['chrdesc/mm8', 'chrdesc/mm9', 'chrdesc/hg18', 'chrdesc/hg19'] },
      packages = ['pyicoslib.lib'],
      scripts = ['pyicos'],
      py_modules = ['pyicoslib.core', 'pyicoslib.operations', 'pyicoslib.parser', 'pyicoslib.defaults']
     )

