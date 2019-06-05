from setuptools import setup, find_packages

setup(name='xMCA',
      version='0.1',
      description='Maximum Covariance Analysis in xarray.',
      url='https://github.com/Yefee/xMCA',
      author='Chengfei He',
      author_email='che43@wisc.edu',
      include_package_data=True,
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        ],
      keywords='statistical analysis MCA SVD xarray',
      license='MIT',
      # packages=['xcesm','config'],
      packages=find_packages(),
      package_data={'xMCA': ['examples/data/*.nc']},
      install_requires=['xarray', 'numpy'],
      zip_safe=False)
print(find_packages())
