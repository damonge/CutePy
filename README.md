# CutePy

## Installation
1. Compile and install the C library. Should be a matter of `./configure; make; make install`. Since you probably don't have admin permissions, you'll need to do `./configure --prefix=/path/to/prefix/`, where `/path/to/prefix/` is wherever you can install stuff (e.g. `$HOME`). If you're installing on NERSC, you may also need to do `export CC=cc` (but maybe not).
2. If the installation folder is not in your path, you'll need to do:
```
export LD_LIBRARY_PATH=/path/to/prefix:$LD_LIBRARY_PATH
export LDFLAGS+=" -L/path/to/prefix/lib"
export CPPFLAGS+=" -I/path/to/prefix/include"
```
3. Install python wrapper. Simply do `python3 setup.py install --user` (I'm still assuming you don't have admin privileges).
