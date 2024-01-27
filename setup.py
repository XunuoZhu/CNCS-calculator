from setuptools import setup,find_packages

requires_packages=['scipy','progressbar']

setup(
    name = 'CNCScalculator',
    version = '2.0',
    packages = find_packages(),
    install_requires = requires_packages,
    #url = 'https://github.com/jiujiezz/candris',
    author = 'Gu Xun & Zhou Zhan',
    author_email = 'xgu@iastate.edu & zhanzhou@zju.edu.cn',
    include_package_data=True
)

