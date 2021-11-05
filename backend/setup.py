from setuptools import setup

setup(
    name='reaction-cime',
    packages=['reaction-cime'],
    include_package_data=True,
    install_requires=[
        'flask',
        'projection-space-explorer @ git+https://github.com/jku-vds-lab/projection-space-explorer.git@develop#egg=projection-space-explorer&subdirectory=backend',
        'flask-cors'
    ],
)
