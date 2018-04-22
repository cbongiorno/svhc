from setuptools import setup

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
	scripts=['bin/svhc','bin/svhc_benchmark','bin/svhc_plot'],
	author='Christian Bongiorno',
	author_email='pvofeta@gmail.com',
	license='GPL',
	zip_safe=False

      )

