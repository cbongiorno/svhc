from distutils.core import setup

setup(	name='svhc',
      	version='0.1',
      	packages=['svhc'],
	install_requires=[	'numpy',
				'pandas',
				'fastcluster',
				'matplotlib',
				'scipy',
				'collections',
				'seaborn',
				'multiprocessing',
				'functools'],
	scripts=['bin/svhc','bin/svhc_benchmark'],
	author_email='Christian Bongiorno',
	email='pvofeta@gmail.com',
	license='MIT',
	zip_safe=False

      )

