F90FLAGS = -fbounds-check -fbacktrace -g
SRC = src
BIN = .bin

# where to locate the object code
#
OBJ= .obj

# Compile flags
FC=gfortran 
LD=gfortran -g
OPT= -c -O3 -Wall -pedantic
DBG= -c -g 

SRC= src
EXECUTABLE= test_GMRES.x

#
# object files
#
_OBJECTS = gmres_module.o conjugate_gradient.o test_GMRES.o
OBJECTS = $(patsubst %,$(OBJ)/%,$(_OBJECTS))

#
# executable 
#
$(EXECUTABLE) : $(OBJECTS)
	$(LD) -o $(EXECUTABLE) $(OBJECTS)
$(OBJ)/%.o : $(SRC)/%.f90 | $(OBJ)
	$(FC) $(F90FLAGS) -c -o $@ $< -J$(OBJ)/
$(OBJ) :
	mkdir -p $(OBJ)

.PHONY : clean

# clean up
clean :
	rm -f $(EXECUTABLE) 
	rm -rf $(OBJ)