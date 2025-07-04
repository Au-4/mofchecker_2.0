from setuptools import setup, find_packages

setup(
    name="mofchecker",
    version="2.0",
    author="Xin Jin",
    author_email="xin.jin@epfl.ch",
    description="MOFChecker is an algorithm designed for MOF structure curation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Au-4/mofchecker_2.0",
    license = "ANCSA",
    license_files=["LICENSE"],
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "click",
        "networkx>=2.5",
        "backports.cached-property",
        "ase",
        "pyyaml",
        "structuregraph_helpers",
        "element_coder",
        "typing_extensions",
        "libconeangle",
        "psutil"
    ],
    extras_require={
        "pyeqeq": [
            "pyeqeq"
        ]
    },
    entry_points={
        "console_scripts": [
            "mofchecker=mofchecker.cli:run",
        ],
    },    
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)

