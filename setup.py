import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
      name='decouphage',
      version='0.0.1',
      description='decouphage - A tool to annotate phage genomes.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      author='Ricardo Pascal',
      author_email='voorloop@gmail.com',
      url='https://github.com/voorloopnul/decouphage',
      packages=setuptools.find_packages(),
      scripts=['decouphage.py'],
      python_requires='>=3.4',
      install_requires=[
            "biopython==1.79",
            "phanotate==1.5.0",
            "pandas==1.3.2",
            "click==8.0.1",
            "pytest==7.1.2"
      ]
)