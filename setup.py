from json import load
from pathlib import Path

from setuptools import find_packages, setup

here = Path(__file__).parent
pkg = load((here / "package.json").open())


def read_it(name):
    fn = here / name
    return fn.read_text() if fn.exists() else ""


def requirements(file):
    return [r.strip() for r in read_it(file).strip().split("\n")]


setup(
    name=pkg["name"].lower(),
    version=pkg["version"].replace("-SNAPSHOT", ".dev0"),
    url=pkg["homepage"],
    description=pkg["description"],
    long_description=read_it("README.md"),
    long_description_content_type="text/markdown",
    keywords=pkg.get("keywords", ""),
    author=pkg["author"]["name"],
    author_email=pkg["author"]["email"],
    license=pkg["license"],
    zip_safe=False,
    entry_points={
        "visyn.plugin": ["{0} = {0}:VisynPlugin".format(pkg["name"].replace("-", "_"))],
    },
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        # Pick your license as you wish (should match "license" above)
        "License :: OSI Approved :: " + ("BSD License" if pkg["license"] == "BSD-3-Clause" else pkg["license"]),
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.10",
    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=["docs", "tests*"]),
    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=requirements("requirements.txt"),
    extras_require={
        "develop": requirements("requirements_dev.txt"),
    },
    # Include all files from the MANIFEST.in file.
    include_package_data=True,
    package_data={},
    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files # noqa
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    data_files=[],  # [('my_data', ['data/data_file'])],
)
