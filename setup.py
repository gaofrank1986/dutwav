from setuptools import setup

setup(name='dutwav',version='0.1',
        description="pre-processing tool for DUTWAV",
        url='',
        author ='frank gao',
        author_email='gaosong2006@gmail.com',
        include_package_data = True,
        license='MIT',
        packages=['dutwav','dutwav.analytical'],#,'wave','mass']
        zip_safe=False)


