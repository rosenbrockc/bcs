LOG		= compile.log
#FFLAGS		= -c
.SILENT:

# sources
WSRCF90		=  \
		f90wrap_num_types.f90 \
		f90wrap_matrix_sets.f90 \
		f90wrap_laplace.f90 \
		f90wrap_bcs.f90 \

SRCF90		=  \
		../bcs/num_types.f90 \
		../bcs/matrix_sets.f90 \
		../bcs/laplace.f90 \
		../bcs/bcs.f90

OBJSF90		= $(SRCF90:.f90=.o)

all: $(OBJSF90)
	-rm _bcs.so
	echo "Creating shared library _bcs.so"
	-f2py-f90wrap -c -m _bcs $(WSRCF90) -I. >> /dev/null
	echo "done."

$(OBJSF90): %.o: %.f90
	echo "Compiling: $^"
	-$(F90) -c $^
	-$(F90) -c $(FFLAGS) $^ -o $@ >> $(LOG) 2>> $(LOG)
	echo "done."

clean:
	-rm *.o *.mod *.so
