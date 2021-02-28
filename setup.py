from setuptools import setup, find_packages

with open('README.md') as readme_file:
    README = readme_file.read()

with open('HISTORY.md') as history_file:
    HISTORY = history_file.read()

setup_args = dict(
    name='molmag_ac_gui',
    version='0.2.0',
    description='A user interface and functions to work with magnetic relaxation',
    long_description_content_type="text/markdown",
    long_description=README + '\n\n' + HISTORY,
    license='MIT',
    packages=find_packages(),
    author='Emil A. Klahn',
    author_email='emil.klahn@gmail.com',
    keywords=['Magnetism', 'Molecular magnetism', 'Magnetic relaxation'],
    url='http://www.eklahn.com',
    download_url='https://pypi.org/project/molmag-ac-gui'
)

#Not best practice to do this, I've added the issue on Github
install_requires = [
    'asteval==0.9.21',
    'cycler==0.10.0',
    'future==0.18.2',
    'kiwisolver==1.3.1',
    'lmfit==1.0.1',
    'matplotlib==3.3.3',
    'numpy==1.19.4',
    'pandas==1.1.4',
    'Pillow==8.0.1',
    'pyparsing==2.4.7',
    'PyQt5==5.15.2',
    'PyQt5-sip==12.8.1',
    'python-dateutil==2.8.1',
    'pytz==2020.4',
    'scipy==1.5.4',
    'six==1.15.0',
    'uncertainties==3.1.5'
]

if __name__ == '__main__':
    setup(**setup_args, install_requires=install_requires)