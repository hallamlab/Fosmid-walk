import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

pks = ['fosmid_walk']

setuptools.setup(
    name="fosmid-walk",
    version="0.1",
    author="Avery Noonan & Tony Liu",
    author_email="contacttonyliu@gmail.com",
    description="estimates population size of fosmid pool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Tony-xy-Liu/fosmid-walk",
    project_urls={
        "Bug Tracker": "https://github.com/Tony-xy-Liu/fosmid-walk/issues",
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
            'foswalk = fosmid_walk:main',
        ]
    },
    python_requires=">=3.10",
)