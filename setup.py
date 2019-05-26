from setuptools import setup

setup(
    name='moltherm',
    version='0.0.1',
    packages=['moltherm', 'moltherm.screen', 'moltherm.screen.tests',
        'moltherm.compute', 'moltherm.compute.tests'],
    url='github.com/peiyuan-yu/Moltherm',
    license='MIT',
    author='Peiyuan Yu, Qi Wang, Evan Spotte-Smith',
    author_email='peiyuan@lbl.gov',
    description='High-throughput generation, optimization and calculation of molecules and chemical reactions.'
)
