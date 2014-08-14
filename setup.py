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
      packages=['PySeidon'],
      package_dir={'PySeidon': 'PySeidon/fvcomClass'},
      package_dir={'PySeidon': 'PySeidon/adcpClass'},
      package_dir={'PySeidon': 'PySeidon/drifterClass'},
      package_dir={'PySeidon': 'PySeidon/stationClass'},
      package_dir={'PySeidon': 'PySeidon/tidegaugeClass'},
      package_dir={'PySeidon': 'PySeidon/validationClass'},
      package_dir={'PySeidon': 'PySeidon/utilities'},
      zip_safe=False)
