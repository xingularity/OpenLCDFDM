# OpenLCDFDM
## What is this?
This is an open source finite difference simulation code for liquid crystal device. This code is developed for multidimensional simulation, it is also developed to perform parallel simulations through OpenMP.
## Purpose of this code
Open source codes are important for creating collaborations, reducing duplicated works and distributing knowledges, they have become essential parts of improving human life and advancing human knowledges and civilizations. There are many open source simulation codes and libraries for various physics and engineering fields such as elastic mechanics, fluid dynamics, plasma physics, EM fields and semiconductor devices. However, there is no open source code to simulate liquid crystal devices even LCD display is now a major technology. This code is developed to fill the gap.
## Simulation methods
There are three major parts in the simulation of liquid crystal devices.

The first part is to simulate liquid crystal transition under external electric field, this transition is called Fréedericksz transition. This part is done by minimizing the Oseen-Frank free energy and solving Poisson's equation.

The second part is to simulate polarization change when light propagates through the liquid crystal device. There are two major methods to be applied here, they are extended Jones matrix method and Berreman 4X4 method.

The last part is colorimetry calculation. It calculates human color perception by the light emitted from LCD devices.
## Libraries and Compilers
1. Libraries
  1. [blitz++](http://sourceforge.net/projects/blitz/)
  2. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  3. Python3 (>=3.4)
    1. Matplotlib
    2. numpy
    3. scipy
    2. Cython
  4. OpenMP
2. Compiler
  1. g++(>=4.8) or clang++, support for C++11.

## License
This code is distributed under BSD license. please read the LICENSE file in the repository.

## Compile and installation
The complete setup.py hasn't been written yet. Users can manually run 
  python3 setup.py build
in the repository root and copy the compiled lcd1d\*.so file to the directory where he/she prepare to run lcd1d. There are examples in "/examples/lcd1d_example/" directory. One can copy the so file to this directory and try example scripts "lcd1d_example_\*.py".

## Example plots


## References
1. Optics of Liquid Crystal Display by Pochi Yeh and Claire Gu. ISBN: 0470181761
2. Fundamentals of Liquid Crystal Devices by Shin-Tson Wu and Deng-Ke Yang. ISBN: 978-0-470-03202-2
