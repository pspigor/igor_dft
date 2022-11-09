import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="igor_dft",
    version="0.0.1",
    author="Igor Pereira", #<<<
    author_email="ipereira@peq.coppe.ufrj.br", #<<<
    description="A small cdft package", #<<<
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pspigor/igor_dft",  #<<<
    packages=setuptools.find_packages(
        where='.',
        include=['igor_dft*'],  # alternatively: `exclude=['additional*']`
        ),
    classifiers=[
        "Programming Language :: Python :: 3", 
        "Operating System :: OS Independent",
    ],
)
