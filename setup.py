from setuptools import setup, find_packages

with open('README.md') as readme_file:
    README = readme_file.read()

with open('HISTORY.md') as history_file:
    HISTORY = history_file.read()

setup_args = dict(
    name='molmag_ac_gui',
    version='0.2.2',
    description='A user interface and functions to work with magnetic relaxation',
    long_description_content_type="text/markdown",
    long_description=README + '\n\n' + HISTORY,
    license='MIT',
    packages=['molmag_ac_gui'],
    author='Emil A. Klahn',
    author_email='emil.klahn@gmail.com',
    keywords=['Magnetism', 'Molecular magnetism', 'Magnetic relaxation'],
    url='http://www.eklahn.com',
    download_url='https://pypi.org/project/molmag-ac-gui'
)

package_data = dict(
    data=["molmag_ac_gui/data/read_options.json"]
)

#Not best practice to do this, I've added the issue on Github
install_requires = [
    'lmfit',
    'matplotlib',
    'numpy',
    'pandas',
    'names'
    'PyQt5',
    'scipy',
]

if __name__ == '__main__':
    setup(**setup_args,
          install_requires=install_requires,
          package_data=package_data,
          include_package_data=True
          )