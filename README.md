Single degree of freedom elastic response
==========================================

Linear newmark method is used for the computation.

The script takes the following input :

T   : Spectral acceleration period for the computation  
z   : Damping ratio (z = 0.005 for 5% damping)  
ag  : The acceleration time history of input ground-motion  
Dtg : The discreet time interval at which the input acceleration is recorded  

The script returns the following:

S  
|-S.d : Spectral displacement  
|-S.v : pseudo spectral velocity [S.d * ( 2 * pi / T )]  
|-S.a : pseudo spectral acceleration [S.d * ( 2 * pi / T ) ^ 2]  

H  
|-H.d : Time history of the displacement response of the oscillator  
|-H.v : Time history of the velocity response of the oscillator  
|-H.a : Time history of the acceleration response of the oscillator  

Compilation
===========

Before using the script in Matlab, the C file needs to be compiled into a mex file.
To compile the C file type "mex linearNewmark_TH.c" in the matlab command line window and follow the onscreen instructions. After compilation you can call the getLinearSpectra_withTH.m as any standard Matlab function.

License
=======

The sofware is released under the MIT license. (http://opensource.org/licenses/MIT)