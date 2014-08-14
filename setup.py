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
      package_dir={'pyseidon': 'PySeidon/fvcomClass'},
      package_dir={'pyseidon': 'PySeidon/adcpClass'},
      package_dir={'pyseidon': 'PySeidon/drifterClass'},
      package_dir={'pyseidon': 'PySeidon/stationClass'},
      package_dir={'pyseidon': 'PySeidon/tidegaugeClass'},
      package_dir={'pyseidon': 'PySeidon/validationClass'},
      package_dir={'pyseidon': 'PySeidon/utilities'},
      zip_safe=False)
