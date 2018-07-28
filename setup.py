from setuptools import setup

setup(	name='svhc',
      	version='0.29',
      	packages=['svhc'],
	install_requires=[	'numpy>=1.14.2',
				'pandas>=0.22.0',
				'fastcluster>=1.1.24',
				'matplotlib>=2.0.2',
				'scipy>=1.0.1',
				'seaborn>=0.8.1',
				'multiprocessing>=2.6.2.1',
				'joblib>=0.12.1'],
	scripts=['bin/svhc','bin/svhc_benchmark','bin/svhc_plot'],
	author='Christian Bongiorno',
	author_email='pvofeta@gmail.com',
	license='GPL',
	zip_safe=False

      )

