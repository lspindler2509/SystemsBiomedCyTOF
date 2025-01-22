from setuptools import setup, find_packages

setup(
    name='acdc',
    version='0.1.0',
    description='Automated Cytometry-based Discovery of Cell Types (ACDC)',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Mahima Arunkumar',
    author_email='mahima.arunkumar@tum.de',
    packages=find_packages(),
    install_requires=[
        'numpy>=1.22.2',
        'pandas>=1.5.2',
        'scikit-learn>=1.2.0',
        'scipy>=1.10.1',
        'seaborn>=0.12.2',
        'matplotlib>=3.7.1',
        'phenograph>=1.5.7', 
        'fcsy>=0.10.0',
        #'flowio>=1.3.0',
        "umap-learn>=0.5.7",
        "plotly>=5.24.1",
        "packaging>=23.1",
        "tenacity>=6.2.0",
        "scanpy==1.10.3"
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
