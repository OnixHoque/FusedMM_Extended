# Shared Object for FusedMM

Shared Object of FusedMM that can be used with Python.

## System Requirements
Users need to have the following software/tools installed in their PC/server. The source code was compiled and run successfully in Linux (Ubuntu and Debian distributions).
```
GCC version >= 4.9
OpenMP version >= 4.5
```
## How to use?

- Clone the repository. Make sure you have GCC installed.
- Select the branch containing shared object code:
    `git checkout python_compatible_portable`
- Run
    `./configure`
- Run
    `./generate_shared_object.sh`
- Run
    `python fusedmm.py`


