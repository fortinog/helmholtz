# Hidden directory to keep all compiled files
OBJ= .obj

# Compile flags
FC=gfortran 
LD=gfortran -g
DBG= -c -g 
F90FLAGS = -fbounds-check -fbacktrace -g

SRC= src
EXECUTABLE= helmholtz_wave_fd.x

# Define object files
# _OBJECTS = gmres_module.o conjugate_gradient.o test_GMRES.o helmholtz_parameters.o helmholtz_wave_fd.o
_OBJECTS = helmholtz_parameters.o helmholtz_wave_fd.o
OBJECTS = $(patsubst %,$(OBJ)/%,$(_OBJECTS))

# Rules to build executable 
$(EXECUTABLE) : $(OBJECTS)
	$(LD) -o $(EXECUTABLE) $(OBJECTS)
$(OBJ)/%.o : $(SRC)/%.f90 | $(OBJ)
	$(LD) $(F90FLAGS) -c -o $@ $< -J$(OBJ)/
$(OBJ) :
	mkdir -p $(OBJ)

.PHONY : clean

# clean up
clean :
	rm -f *.x 
	rm -rf $(OBJ)