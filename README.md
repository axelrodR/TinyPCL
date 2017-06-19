# TinyPCL
A slim point cloud library with minimal dependencies.
This library contains work contributed by developers and researchers in the Omek consortium
http://www.omek3d.org/about-omek

### **work in progress**

## To build:
TinyPCL uses the [premake](https://premake.github.io/) cross platform build system.
To build, first download premake5 from [here](https://premake.github.io/download.html) and put it in ./make
call premake5 with the target platform (e.g. premake5 vs2015).
For more information see: https://premake.github.io/

## Structure:

| dir               | content                    |
|-------------------|----------------------------|
| ./include/        | external library interface |
| ./src/common/     | common untilities          |
| ./src/*           | internal interface and implementations |


