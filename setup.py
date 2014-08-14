from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='PySeidon',
      version='v0.1',
      description='Suite of tools for FVCOM model',
      long_description=readme(),
      url='https://github.com/GrumpyNounours/PySeidon',
      author='Thomas Roc, Wesley Bowman, Jon Smith',
      author_email='thomas.roc@acadiau.ca,wesley.bowman23@gmail.com,lavieenroux20@gmail.com',
      maintainer='Thomas Roc',
      license='GNU Affero GPL v3.0',
      packages=['pyseidon'],
      package_dir={'pyseidon': ['PySeidon/fvcomClass','PySeidon/adcpClass', 'PySeidon/drifterClass','PySeidon/stationClass', 'PySeidon/tidegaugeClass','PySeidon/validationClass','PySeidon/utilities']},
      zip_safe=False)
