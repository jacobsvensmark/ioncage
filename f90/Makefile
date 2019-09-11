############################################################################################
# MAKEFILE - By addition of new sourcecode in seperate files, make sure to add it here.    #
# Note that any change in source code requires a total recompilation using this make file. #
# To compile into an executable, type                                                      #
# > make                                                                                   #
# and the specified .out executable file is produced.                                      #
############################################################################################
FC              = gfortran-mp-8
CFLAGS          = -O2 -mmacosx-version-min=10.14 -o 
SOURCES         = precision_type.f90 math_module.f90 calc_module.f90 io_general_module.f90 \
                  io_ioncage_module.f90 dNdt_module.f90 main.f90
EXECUTABLE_NAME = ioncage.out
all: main.f90
	$(FC) $(CFLAGS) $(EXECUTABLE_NAME) $(SOURCES)
