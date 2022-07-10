# Makefile for compiling the Serial Spectral Code Version 2.1 
# This uses the Intel Fortran Compiler, with FFTW 
# ---------------------------------------------------------------------
compile = ifort -openmp -O3 -traceback
link = -lrfftw_threads -lfftw_threads -lrfftw -lfftw  -L/home/mohit.gupta/fftw_install/fftw2/lib
#compile = gfortran -g -fopenmp -fbacktrace
#link = -lrfftw_threads -lfftw_threads -lrfftw -lfftw  -L/home/jason/Downloads/fftw/fftw-2.1.5/lib
#compile = ifort -openmp -O3
#link = -lrfftw_threads -lfftw_threads -lrfftw -lfftw  -L/home/ssray/softwares/fftw-2.1.5/lib
#link = -lrfftw -lfftw  -L/home/rama/software/fftw2/lib -I/home/rama/software/fftw2/i
#---------------------------------------------------------------------
OBJ_rk =  2dflurk.o global2drk.o gen_force.o iniconf2d.o energy2d.o rk2.o \
           eval_jac.o apply_pbc.o getpsivel.o getrandom.o 
MOD_FILE= mod_2dflu.o
MOD_FILE1= mod_part_interp.o
test : 
	echo $(OBJ_FILE) 
# -------------------------------------------------------------
nsrk.exe:$(OBJ_rk) $(MOD_FILE) $(MOD_FILE1)
	$(compile) $(OBJ_rk) $(MOD_FILE) $(MOD_FILE1) -o nsrk.exe $(link)
#
2dflurk.o: 2dflurk.f90 $(MOD_FILE) $(MOD_FILE1) 
	$(compile) -c 2dflurk.f90 $(MOD_FILE) $(MOD_FILE1)
#         
global2drk.o: global2drk.f90 $(MOD_FILE) 
	$(compile) -c global2drk.f90 $(MOD_FILE) 
#
trumbbell.o: trumbbell.f90 $(MOD_FILE)
	$(compile) -c trumbbell.f90 $(MOD_FILE)
#
iniconf2d.o: iniconf2d.f90 $(MOD_FILE) 
	$(compile) -c iniconf2d.f90 $(MOD_FILE) 
#         
gen_force.o: gen_force.f90 $(MOD_FILE) 
	$(compile) -c gen_force.f90 $(MOD_FILE)
#                 
#forcen.o: forcen.f90 $(MOD_FILE) 
#	$(compile) -c forcen.f90 $(MOD_FILE)
#         
energy2d.o: energy2d.f90 $(MOD_FILE) 
	$(compile) -c energy2d.f90 $(MOD_FILE) 
#         
eval_jac.o: eval_jac.f90 $(MOD_FILE) $(MOD_FILE1) 
	$(compile)   -c   eval_jac.f90 $(MOD_FILE) $(MOD_FILE1) 
#
apply_pbc.o: apply_pbc.f90 $(MOD_FIlE) 
	$(compile) -c apply_pbc.f90 $(MOD_FILE)
#         
rk2.o: rk2.f90  $(MOD_FILE) 
	$(compile) -c rk2.f90 $(MOD_FILE) 
#       
getpsivel.o: getpsivel.f90 $(MOD_FILE)
	$(compile) -c getpsivel.f90 $(MOD_FILE)
#
getrandom.o: getrandom.f90 $(MOD_FILE)
	$(compile) -c getrandom.f90 $(MOD_FILE)
#
#plot_cont.o: plot_cont.f 
#		$(compile) -c plot_cont.f
#
mod_2dflu.o:  mod_2dflu.f90 
	$(compile) -c mod_2dflu.f90
mod_part_interp.o:  mod_part_interp.f90 
	$(compile) -c mod_part_interp.f90
# ------------------------------------------------------------------------ 
clean:
	rm -f *.o *.mod core *.cmp.f90 *.kmo work.* *.d
total_clean:
	rm -f *.o *.mod core *.exe
data_clean:
	rm -f *.dat *.out data_psi/* spectras/* 
## --------------------------------------------------------------------##  
