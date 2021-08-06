target = ../main
object = \
  variables.o \
  myutil.o \
  myPBiCGSTAB_MPI.o \
  NewtonRaphson.o \
  sorsolver.o \
  bdrc.o \
  psifunc.o \
  gs_solver.o \
  Field.o \
  RgridMod.o \
  BgridMod.o \
  diagnosis.o \
  equset.o \
  writingMod.o \
  initialMod.o \
  main.o

# trisolver.o \
# spline_interp.o \
# NR_extremum.o \

# FC = gfortran

FFLAGS  =
LDFLAGS =

ifeq (${shell uname},Linux)
  ifeq (${shell hostname},fes1)
    FC = frtpx
    FFLAGS  += -Kfast
    LDFLAGS += -Kopenmp
  endif
  ifeq (${shell hostname},kolmo01)
    FC = mpiifort
    LDFLAGS += -qopenmp
  endif
  ifeq (${shell hostname},kolmogorov)
    FC = mpiifort
    FFLAGS  += -O3
    LDFLAGS += -qopenmp
  endif
  ifeq (${shell hostname},gauss)
    FC = mpif90
    LDFLAGS += -O3 -fopenmp
  endif
  ifeq (${shell hostname},evans)
    FC = mpifort
    LDFLAGS += -O3 -openmp
  endif
  ifeq (${shell hostname},FrenchQuarter)
    FC = mpifort
    LDFLAGS += -openmp
  endif
endif
ifeq (${shell uname},Darwin)
  FC = mpifort
  LDFLAGS += -qopenmp
  FFLAGS  += -Wl,-stack_size,0x40000000
endif
# ---  Debug  Options  --- #
# FFLAGS = -fbacktrace -ffpe-trap=invalid,zero,overflow
# FFLAGS = -g -traceback -CB
# FFLAGS = -g -traceback -check-pointers=rw
# FFLAGS = -g -traceback -check all -warn all -fpe0
# -----------------   Option End   ------------------ #



.SUFFIXES: .o .f90

%.o: %.f90
	$(FC) -c $< $(FFLAGS) $(LDFLAGS)

%.mod: %.f90 %.o
	@:

$(target): $(object)
	$(FC) -o $@ $(object) $(FFLAGS) $(LDFLAGS)

.PHONY:all clean
clean:
	$(RM) $(target) $(object) *.mod
all: clean $(target)
