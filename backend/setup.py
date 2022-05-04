from setuptools import setup

setup(
    name='reaction_cime',
    packages=['reaction_cime'],
    include_package_data=True,
    install_requires=[
        'flask==2.0.1', # 2.0.1; set to 2.0.2 later --> 2.0.2 does not show error traceback.... (or 2.0.3 as soon as it is released)
        'flask_sqlalchemy==2.5.1',
        'sqlalchemy==1.4.30', # 1.3.20 # 1.4 includes caching
        'projection-space-explorer @ git+https://github.com/jku-vds-lab/projection-space-explorer.git@develop-cime#egg=projection-space-explorer&subdirectory=backend',
        'flask-cors==3.0.10',
        'rdkit-pypi==2021.9.2',
        'opentsne==0.6.1',
        'umap-learn==0.5.2',
        'scikit-learn',
        'gower==0.0.5',
        'pandas',
        'Pillow',
        'numba==0.55.1'
    ],
)
