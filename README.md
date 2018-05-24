# Fisher-matrix-forecast-for-CMB-B-mode

Please cite arXiv:1605.03504 if you use this code or a modified version based on it.


This is a description of the iFisherCMB program. This code needs camb pre-installed, and needs the outputs from another code CMB-foreground-residuals.

Tensor mode need to be used

A)to run the program
  0. Running script fisherCMB is to run step 1-3 together. But they can be run 
     individually as following.
  1. inputgencamb.out is to generate inputs of camb.
  2. genspect is a script to generate powerspectrum with all inputs generated in step 1.
  3. fisherCMB.out is to generate the fisher matrix, which will be stored in fisher.txt
  
B)to compile (always compile with following)
  1. ifort -o inputgencamb.out inputgencamb.f90
  2. ifort -o fisherCMB.out fisherCMB.f90
  
C)to add more fiducial parameters
  1. need to modified inputgencamb.f90, and a little in fisherCMB.f90
  2. change the section A ,B and C (not here) in inputgencamb.f90 according to
      the descriptions. NOTICE: the order of parameters in the three sections need
	  to be consistent.
  3. Also add the executing commands in fisherCMB script too (or genspect script)
  4. After adding, compile with setup B)1., and run.
  5. At the end of fisherCMB.f90, add a format for the new parameter.
