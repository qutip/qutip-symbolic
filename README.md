# qutip-symbolic

[![build](https://github.com/qutip/qutip-symbolic/workflows/Tests/badge.svg)](https://github.com/qutip/qutip-symbolic/actions)
[![Documentation Status](https://readthedocs.org/projects/qutip-symbolic/badge/?version=latest)](https://qutip-symbolic.readthedocs.io/en/latest/)
[![PyPI version](https://badge.fury.io/py/qutip-symbolic.svg)](https://badge.fury.io/py/qutip-symbolic)

The qutip-symbolic package is a continuation of [SymPsi](https://github.com/sympsi/sympsi), which is itself a fork of [SymPy's](https://github.com/sympy/sympy/) ``sympy.physics.quantum`` package.

In contrast to SymPsi, qutip-symbolic aims to be an extension of ``sympy.physics.quantum``, rather than a fork. The intention is to slowly contribute changes back to SymPy while maintaining a compatibility layer so that users of qutip-symbolic do not need to continually update their ``from qutip_symbolic import ...`` statements.

Installation
------------

To install the package, download to source code and run

```
pip install qutip-symbolic
```

If you want to edit the source code, please download the source code and run the following command under the folder with `setup.cfg`.

```
pip install -e .
```

To build and test the documentation, additional packages need to be installed:

```
pip install matplotlib sphinx numpydoc sphinx_rtd_theme
```

Under the `docs` directory, use

```
make html
```

to build the documentation, or

```
make doctest
```

to test the code in the documentation.

Documentation
-------------

The latest documentation of `qutip-symbolic` hosted at [qutip-symbolic.readthedocs.io/](https://qutip-symbolic.readthedocs.io/en/latest/).

Testing
-------

To test the installation from a download of the source code, run from the `qutip-symbolic` directory

```
pytest tests
```

Support
-------

This package is supported and maintained by the same developers group as QuTiP.

[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-NumFOCUS-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://numfocus.org)
[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=flat)](https://unitary.fund)


QuTiP development is supported by [Nori's lab](http://dml.riken.jp/)
at RIKEN, by the University of Sherbrooke, by Chalmers University of Technology, by Macquarie University and by Aberystwyth University,
[among other supporting organizations](http://qutip.org/#supporting-organizations).

License
-------

[![license](https://img.shields.io/badge/license-New%20BSD-blue.svg)](http://en.wikipedia.org/wiki/BSD_licenses#3-clause_license_.28.22Revised_BSD_License.22.2C_.22New_BSD_License.22.2C_or_.22Modified_BSD_License.22.29)

You are free to use this software, with or without modification, provided that the conditions listed in the LICENSE.txt file are satisfied.
