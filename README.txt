Guidelines how to use the code: 

REQUIREMENTS:
 - Cmake version no less than 3.7.2
 - A suitable C++ compiler
 - The linear algebra library Eigen
 - The Boost libraries

 
Further the code uses the open source C++-package on cubic spline interpolation included in the class spline.h which is from:
https://kluge.in-chemnitz.de/opensource/spline/ 


After compilation using cmake and make, a target called "s2d" is produced 
which can be executed by: el3_s3d -f [path to]/cylinder=14.obj -c  [path to]/default2.cfg

By this the simulation is started on a 3d-cylinder domain O_h

The main code of the simulation is provided in the executable 
file "pde.cpp" and "main.cpp".

In "main.cpp" the relevant parameters are described 
in lines 12-20, where e.g. the different standard deviations in the fuzzy
representation of the geometric perturbations are set

The default configuration of the code is appropriate to
the 3d-simulation provided in the thesis Section 15.
