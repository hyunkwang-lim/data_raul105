# Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
# Licensed under the GRASP Open Source License V1.0 (see LICENSE file)

# Options for C compiler
CC = gcc
CLIBS= -L/usr/local/lib -I/usr/local/include -I/usr/local/grasp-deps/include -L../../libs/lib -L../settings/yamlsettings -Iyamlsettings
LIBPREFIX=/usr/local
CFLAGS=-std=c99 -g -Wall $(CLIBS)
FC = gfortran
f90 = $(FC)

EXE = #grasp_input_load_functions_template_rule #drivers/sdata/sdata_dump drivers/polder/test_polder_driver


INPUT_OBJECTS= grasp_input.o grasp_input_preloader.o grasp_input_load_functions.o grasp_driver_settings.o grasp_input_segment.o
DRIVERS_OBJECTS=$(shell ./load_drivers.sh)
TRANSFORMERS_OBJECTS=$(shell ./load_transformers.sh)
OBJECTS=$(DRIVERS_OBJECTS) $(TRANSFORMERS_OBJECTS) $(INPUT_OBJECTS) 

ifeq ($(notdir $(f90)), gfortran)
FORTRAN_LIBS = -L$(LIBGFORTRAN) -lgfortran
else
ifeq ($(notdir $(f90)), ifort)
# the -shared flag is necessary with gcc for avoiding an error "hidden symbol ... referenced by DSO"
# when linking with the ifort libraries
# ref. http://www.linuxquestions.org/questions/linux-software-2/gcc-linker-error-hidden-symbol-%60stat%27-in-usr-lib-libc_nonshared-a-stat-os-755962/
FORTRAN_LIBS = -L$(LIBIFORT) -lifcore -limf -lintlc -lifport -ldl -lsvml
endif
endif


all : $(EXE) grasp_input.o grasp_input_preloader.o libgrasp_input.a 


libgrasp_input.a : $(OBJECTS) grasp_input_load_function_template_rule
	ar cru libgrasp_input.a $(OBJECTS) 

grasp_input.o: grasp_input.c ../retrieval/src/mod_par_inv.inc ../retrieval/src/mod_par_OS.inc grasp_input_load_function_template_rule compile_functions
	$(CC) $(CFLAGS) -c grasp_input.c ../retrieval/src/mod_par_inv.inc ../retrieval/src/mod_par_OS.inc 

grasp_input_load_functions.o: grasp_input_load_functions.c ../retrieval/src/mod_par_inv.inc ../retrieval/src/mod_par_OS.inc grasp_input_load_function_template_rule
	$(CC) $(CFLAGS) -c grasp_input_load_functions.c ../retrieval/src/mod_par_inv.inc ../retrieval/src/mod_par_OS.inc
	
grasp_input_load_function_template_rule:
	./grasp_input_load_functions_template.sh > grasp_input_load_functions.c
	./grasp_input_functions_settings_template.sh > grasp_input_functions_settings.h
	
compile_functions:
	for driver in `ls drivers` ; do \
	    cd drivers/$$driver ; \
	    make CFLAGS="$(CFLAGS)" LDFLAGS="$(LDFLAGS)" f90=$(f90); \
	    cd - ; \
	done	
	for transformer in `ls transformers` ; do \
	    cd transformers/$$transformer ; \
	    make CFLAGS+="$(CFLAGS)" LDFLAGS+="$(LDFLAGS)" f90=$(f90); \
	    cd - ; \
	done		

grasp_input_preloader.o: grasp_input_preloader.c grasp_input_load_function_template_rule
	$(CC) $(CFLAGS) -c grasp_input_preloader.c -lm
	
clean:
	-rm -f $(EXE) *.o *.mod *~ *.a *.s */*.o */*/*.o *.dSYM
	-cd drivers/polder && make clean
	-cd drivers/sdata && make clean
	-rm -f grasp_input_load_functions.c grasp_input_functions_settings.h
	for driver in `ls drivers` ; do \
	    cd drivers/$$driver ; \
	    make clean; \
	    cd - ; \
	done
	for transformer in `ls transformers` ; do \
	    cd transformers/$$transformer ; \
	    make clean; \
	    cd - ; \
	done	

