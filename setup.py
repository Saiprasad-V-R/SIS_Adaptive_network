from setuptools import setup, find_packages

setup(
    name="sis_adaptive",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["networkx","numpy","scipy","matplotlib"],
    python_requires=">=3.8",
)

