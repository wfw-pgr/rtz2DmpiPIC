# ------------ #
# --  File  -- #
# ------------ #
target = ../main
object = \
  variables.o  \
  utilityMod.o \
  initialMod.o \
  mPBiCGSTAB.o \
  cylPoisMod.o \
  xyzPoisMod.o \
  BiCG4GSMod.o \
  EgridRZMod.o \
  RgridRZMod.o \
  BgridRZMod.o \
  jxBchckMod.o \
  lInterpMod.o \
  writingMod.o \
  main.o
# calculation or calc_type2
# -----------------    File  End   ------------------ #
  # BgridXZMod.o \


# --------------- #
# --  Options  -- #
# --------------- #
FC      = mpif90
FFLAGS  =
LDFLAGS =
# ===  Linux  === #
ifeq (${shell uname},Linux)
# --- kolmogorov --- #
  ifeq (${shell hostname},kolmogorov)
    FC       = mpiifort
    # FFLAGS  += -g -traceback -O1 -mcmodel=medium -shared-intel
    FFLAGS  += -O1 -mcmodel=medium -shared-intel
    LDFLAGS += -qopenmp
  endif
# ---   evans   --- #
  ifeq (${shell hostname},evans)
    FC       = mpiifort
    FFLAGS  += -fast
    LDFLAGS += -openmp
  endif
# ---   france  --- #
  ifeq (${shell hostname},FrenchQuarter)
    FC       = mpifort
    FFLAGS  += -O3
    LDFLAGS += -openmp
  endif
# ---   gauss   --- #
  ifeq (${shell hostname},gauss)
    FC       = mpif90
    FFLAGS  += -O3
    LDFLAGS += -fopenmp
  endif
endif
# ===  Mac OS  === #
ifeq (${shell uname},Darwin)
  FC = mpifort
  FFLAGS += -xHOST -O3 -no-prec-div -ipo
  FFLAGS  += -Wl,-stack_size,0x80000000
  LDFLAGS += -qopenmp
endif
# ---  Debug  Options  --- #
# FFLAGS = -fbacktrace -ffpe-trap=invalid,zero,overflow
# FFLAGS = -g -traceback -CB
# FFLAGS = -g -traceback -check-pointers=rw
# FFLAGS = -g -traceback -check all -warn all -fpe0
# -----------------   Option End   ------------------ #



# -------------- #
# -- MakeFile -- #
# -------------- #

# -- Rules -- #
.SUFFIXES: .o .f90
%.o: %.f90
	$(FC) -c $< $(FFLAGS) $(LDFLAGS)

%.mod: %.f90 %.o
	@:

$(target): $(object)
	$(FC) -o $@ $(object) $(FFLAGS) $(LDFLAGS)

# -- Orders -- #
.PHONY:clean all
all: clean $(target)
clean:
	$(RM) $(target) $(object) *.mod
# -----------------   Makefile End   ------------------ #
