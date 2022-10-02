# ------------ #
# --  File  -- #
# ------------ #
target = ../main
object = \
  variables.o  \
  shapeFnMod.o \
  mt19937Mod.o \
  randGenMod.o \
  sortingMod.o \
  lInterpMod.o \
  myBiCGSTAB.o \
  NwtnRphMod.o \
  fMPIComMod.o \
  rMPIComMod.o \
  pMPIComMod.o \
  jMPIComMod.o \
  displayMod.o \
  ptCheckMod.o \
  depBdrcMod.o \
  mxwBdrcMod.o \
  pcwBdrcMod.o \
  jcrBdrcMod.o \
  maxwellMod.o \
  smoothfMod.o \
  p2ndMomMod.o \
  momentsMod.o \
  ptcBdrcMod.o \
  collideMod.o \
  qchargeMod.o \
  currentMod.o \
  makeVDFMod.o \
  ptExchgMod.o \
  ptCoalsMod.o \
  ptSplitMod.o \
  ptOrderMod.o \
  kdtree2Mod.o \
  ppcCtrlMod.o \
  ptAccelMod.o \
  debugerMod.o \
  writingMod.o \
  setupGSMod.o \
  allocatMod.o \
  com4DDDMod.o \
  exchDDDMod.o \
  diffDDDMod.o \
  utilityMod.o \
  saveRunMod.o \
  initialMod.o \
  diagnosMod.o \
  recRateMod.o \
  getXOptMod.o \
  sFluxIpMod.o \
  sEnergyMod.o \
  main.o
# -----------------    File  End   ------------------ #
# dynmDomMod.o \
  

# --------------- #
# --  Options  -- #
# --------------- #
# common
FFLAGS  = 
LDFLAGS =
HOST:=${shell hostname}
RB:=${words ${filter reedbush-u%,${HOST}} }
PS:=${words ${filter fes%,${HOST}} }
# ---    ReedBush-U    --- #
ifeq (${RB},1)
  FC       = mpiifort
  FFLAGS  += -xHOST -O3 -no-prec-div -ipo -qopenmp
  LDFLAGS +=
endif
# --- Plasma Simulator --- #
ifeq (${PS},1)
  FC       = mpifrtpx
  FFLAGS  += -Kfast -Kopenmp
  LDFLAGS +=
endif
# ---      gauss       --- #
ifeq (${shell hostname},gauss)
  FC       = mpif90
  FFLAGS  += -O3 -fopenmp
  LDFLAGS +=
endif
# ---      evans       --- #
ifeq (${shell hostname},evans)
  FC       = mpifort
  FFLAGS  += -xHOST -O3 -no-prec-div -ipo -openmp
  LDFLAGS +=
endif
# ---  FrenchQuarter   --- #
ifeq (${shell hostname},FrenchQuarter)
  FC       = mpifort
  FFLAGS  += -xHOST -O3 -no-prec-div -ipo -openmp
  LDFLAGS +=
endif
# ---    kolmogorov    --- #
ifeq (${shell hostname},kolmogorov)
  FC       = mpiifort
  FFLAGS  += -xHOST -O3 -no-prec-div -ipo -qopenmp
  LDFLAGS += 
endif
# ---      euler       --- #
ifeq (${shell hostname},euler.local)
  FC=mpifort
  FFLAGS+=-xHOST -O3 -no-prec-div -ipo
  FFLAGS+=-Wl,-stack_size,0x80000000
  LDFLAGS+=-qopenmp
endif
ifeq (${shell hostname},maxwell.local)
  FC       =mpif90
#  FFLAGS  +=-O3 -fopenmp -fallow-argument-mismatch -fno-range-check
  FFLAGS  +=-O3 -fopenmp -fallow-argument-mismatch
  LDFLAGS +=
endif
ifeq (${shell hostname},maxwell.lan)
  FC       =mpif90
#  FFLAGS  +=-O3 -fopenmp -fallow-argument-mismatch -fno-range-check
  FFLAGS  +=-O3 -fopenmp -fallow-argument-mismatch
  LDFLAGS +=
endif

# ---  Debug  Options  --- #
# FFLAGS += -g
# FFLAGS   += -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow
# FFLAGS += -fbacktrace -ffpe-trap=invalid,zero,overflow
# FFLAGS += -g -traceback
# FFLAGS += -g -traceback -check-pointers=rw
# FFLAGS += -g -traceback -check all -warn all -fpe0
# -----------------   Option End   ------------------ #


# -------------- #
# -- MakeFile -- #
# -------------- #

# -- files -- #
FILES =src main
JOB   =job/$(shell echo `grep -E "::\s+job\s+=" variables.f90 | cut -d "'" -f2`)/

# -- Rules -- #
.SUFFIXES: .o .f90
%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(LDFLAGS)

%.mod: %.f90 %.o
	@:

$(target): $(object)
	$(FC) $(object) -o $@ $(FFLAGS) $(LDFLAGS)

# -- Orders -- #
.PHONY:all clean copy run
all: clean $(target)
clean:
	$(RM) $(target) $(object) *.mod ../run.sh
copy:
	cd ../; cp -r $(FILES) $(JOB)

# -----------------   Makefile End   ------------------ #
