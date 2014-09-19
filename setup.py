from distutils.core import setup

setup(name='plldesigner',
      version='0.1',
      py_modules=['plldesigner', 'plldesigner.sdmod','plldesigner.pll',
                  'plldesigner.pnoise'],
      description='Python toolkit for PLLs design and phase noise analyisis',
      author='Juan F. Osorio',
      author_email='jfosorio@gmail.com',
      url='https://github.com/jfosorio/plldesigner'
     )
