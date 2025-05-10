import setuptools
import glob
# with open("README.md", "r") as fh:
#     long_description = fh.read()


setuptools.setup(
    name="obshelper", # Replace with your own username
    version="1.0.0	",
    author="Xu Chen",
    author_email="xuchen@nao.cas.cn",
    description="A package to assist with planning HI spectral line observations for FAST",
#     long_description_content_type="text/markdown",
#     url="",
    packages=['obshelper', 'obshelper.utils'],
    
    install_requires=['numpy>=1.12', 'matplotlib', 'scipy', 'pandas',
             'astropy', 'astroquery', 'ephem', 'hiviewer', 'ipyaladin' 
],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    python_requires='>=3.6',
    zip_safe=False,
)
