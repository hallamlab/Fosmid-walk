import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

PKG_NAME = "fosmid-walk"
PKG_DIR = PKG_NAME.replace('-', '_')

pks = [PKG_DIR]
setuptools.setup(
    name=PKG_NAME,
    version="1.0",
    author="Avery, Kat, Tony",
    author_email="contacttonyliu@gmail.com",
    description="estimator for fosmid pool population size",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url=f"https://github.com/Tony-xy-Liu/{PKG_NAME}",
    project_urls={
        "Bug Tracker": f"https://github.com/Tony-xy-Liu/{PKG_NAME}/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    # packages=setuptools.find_packages(where="src"),
    packages=pks,
    # package_data={
    #     # "":["*.txt"],
    #     # "package-name": ["*.txt"],
    #     # "test_package": ["res/*.txt"],
    # },
    entry_points={
        'console_scripts': [
            f'foswalk = {PKG_DIR}:main',
        ]
    },
    python_requires=">=3.10",
)
