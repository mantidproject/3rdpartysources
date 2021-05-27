# Fortran

Fortran code for Indirect and Muon utilities.

To build you'll need Python, Numpy and Fortran and C compilers. 

Using conda is the recommended way of building. On windows the following works:

```
conda create --name fortran python=3.8
conda activate fortran
conda install numpy
conda install -c msys2 m2w64-gcc-libgfortran
mkdir build
cd build
cmake ..
cmake --build .
```

And on OSX and Linux

```
conda create --name fortran python=3.8
conda activate fortran
conda install numpy
conda install -c conda-forge fortran-compiler
mkdir build
cd build
cmake ..
cmake --build .
```

The python extension libraries will be placed in the ```build/bin``` folder.


