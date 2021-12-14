from setuptools import setup

setup(
    name='reaction_cime',
    packages=['reaction_cime'],
    include_package_data=True,
    install_requires=[
        'flask',
        'flask_sqlalchemy',
        'projection-space-explorer @ git+https://github.com/jku-vds-lab/projection-space-explorer.git@develop-cime#egg=projection-space-explorer&subdirectory=backend',
        'flask-cors',
        'rdkit-pypi',
        'opentsne',
        'umap-learn',
        'scikit-learn',
        'gower'
    ],
)
