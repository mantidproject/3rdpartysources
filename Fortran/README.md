# Fortran

Fortran code for Indirect and Muon utilities.

To build you'll need Python, Numpy and Fortran and C compilers. 

Using conda is the recommended way of building. On windows the following works:

```
conda create --name fortran python=3.8
conda activate fortran
conda install numpy
conda install -c msys2 m2w64-gcc-fortran
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



Alternatively, you can build the libraries using the setup.py file in the Fortran/Indirect folder. To use this, create the minimal conda environment as shown above. Now do:


```
cd Fortran/Indirect
python setup.py bdist_wheel
```

to build a wheel. This wheel can be uploaded using twine:

```
twine upload ./dist/name_of_wheel
```

Linux wheels require a docker image. For more details see this blog https://uwekorn.com/2019/09/15/how-we-build-apache-arrows-manylinux-wheels.html
