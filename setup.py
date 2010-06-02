from distutils.core import setup

setup(name='Pyicos',
      version='0.8.4',
      description='Mapped reads analysis tool and library',
      author=u'Juan Gonzalez_Vallinas',
      author_email='juanramon.gonzalezvallinas@upf.edu',
      url='http://regulatorygenomics.upf.edu/pyicos',
      packages = ['pyicoslib.lib'],
      scripts = ['pyicos'],
      py_modules = ['pyicoslib.core', 'pyicoslib.operations', 'pyicoslib.parser']
     )

