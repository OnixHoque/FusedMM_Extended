# Shared Object for FusedMM

Shared Object of FusedMM that can be used with Python.

- Clone the repository. Make sure you have GCC installed.
- Select the branch containing shared object code:
    `git checkout python_compatible_portable`
- Run
    `./configure`
- Run
    `./generate_shared_object.sh`
- Run
    `python3 test.py`
