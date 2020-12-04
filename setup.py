from setuptools import setup

# load the __version__ variable without importing the package already
exec(open('mkpy3/version.py').read())

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='mkpy3',
      version=__version__,
      description='Python3 tools for the NASA TESS/K2/KEPLER astrophysical missions',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='http://github.com/KenMighell/mkpy3',
      author='Kenneth Mighell',
      author_email='kmighell@seti.org',
      license='MIT',
      packages=['mkpy3'],
      zip_safe=False,
      classifiers=[
          "Development Status :: 4 - Beta",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Astronomy",
          ],
      python_requires='>=3.6',
)

# EOF
