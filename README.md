# md-scripts
Scripts for MD simulation

## top2psf.py  
A Python code generating a PSF file from a GROMACS top file including gromacs itp files.  
The command-line to excute the code will look like:  

    python3 top2psf.py -f topol.top -o output.psf

## bending_modulus
includes C++ codes to estimate bending modulues of lipid bilayers by calculating spectrum of longitudinal molecular
orientation fluctuations (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.109.028102).  
A binary file `bend_mod` to do spectrum calculation can be compiled as follows:   

    cd src/

and  

    make

The command-line to use the generated file will be:

    ./bend_mod -f input.txt

`input.txt` has some lines to specify analyzed files and lipids. See a sample file in `example` directory.  
After getting spectrum data with `bend_mod`, the bending modulus is estimated with the data around small wave numbers. 
A python code `coeff.py` can be used to see the spectrum data and to calculate the bending modulus.  
The example command-line working in `post` directory will be:  

    python3 coeff.py ret_long.dat 310 3

`ret_long.dat` is the spectrum data obtained using `bend_mod`. The second and third arguments are temperature 
and the number of the spectrum data used to compute the bending modulus.

## lj3d_education  
Training codes (for bachelor students in our group) to run MD simulation of simple lj fluid.  
