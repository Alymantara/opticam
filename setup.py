from setuptools import setup

setup(
  name = 'opticam',         # How you named your package folder (MyLib)
  packages = ['opticam'],   # Chose the same as "name"
  version = '1.0',      # Start with a small number and increase it with every change you make
  license='AGPLv3',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Exposure time calculator for OPTICam/SPM Observatory',   # Give a short description about your library
  author = 'Alexander Stone-Martinez',                   # Type in your name
  author_email = 'stonemaa@nmsu.edu',      # Type in your E-Mail
  url = 'https://github.com/Opticam',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/APOExposureTimeCalculator/APOExptime/archive/v_2.2.tar.gz',    # I explain this later on
  keywords = ['Exposure time', 'Astronomy', 'Opticam'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'astropy',
          'scipy',
          'numpy',
          'pyyaml',
          'wheel',
          'synphot'
      ],
  classifiers=[
    'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Programming Language :: Python :: 3.6',
  ],
  include_package_data=True,
  python_requires='>=3.6'
)
