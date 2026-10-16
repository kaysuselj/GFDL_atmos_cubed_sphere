# Single-face TL/AD test build (CPU, r8, no FMS/MPI). The forward files are compiled unchanged
# from ../ against the stubs in stubs/.
#   make            forward library + test driver
#   make taf        regenerate the TAF TL/AD code in taf/ (needs the FastOpt server, see taf/README)
FC      = gfortran
FFLAGS  = -fdefault-real-8 -fdefault-double-8 -O2 -g -fbacktrace -cpp -ffree-line-length-none
DEBUG   = -O0 -fcheck=bounds,pointer -ffpe-trap=invalid,zero,overflow
B       = build

FWD     = ../tp_core.F90 ../fv_fill.F90 ../fv_mapz.F90 ../fv_tracer2d.F90
FWD_O   = $(B)/stubs.o $(B)/tp_core.o $(B)/fv_fill.o $(B)/fv_mapz.o $(B)/fv_tracer2d.o

DRV_O   = $(B)/adj_face_mod.o $(B)/adj_tops.o
TAF_SRC = $(wildcard taf/*_tl.f90 taf/*_ad.f90)
TAF_O   = $(patsubst taf/%.f90,$(B)/%.o,$(TAF_SRC))

all: $(B)/test_adv

$(B):
	mkdir -p $(B)

$(B)/stubs.o: stubs/stubs.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/%.o: ../%.F90 $(B)/stubs.o
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/tp_core.o: $(B)/stubs.o
$(B)/fv_fill.o: $(B)/stubs.o
$(B)/fv_mapz.o: $(B)/fv_fill.o
$(B)/fv_tracer2d.o: $(B)/tp_core.o $(B)/fv_mapz.o $(B)/fv_fill.o

$(B)/%.o: driver/%.F90 $(FWD_O)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/adj_tops.o: $(B)/adj_face_mod.o
$(B)/%.o: taf/%.f90 $(B)/adj_tops.o
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/test_adv.o: $(TAF_O)

$(B)/test_adv: $(FWD_O) $(DRV_O) $(TAF_O) $(B)/test_adv.o
	$(FC) $(FFLAGS) -o $@ $^

debug: FFLAGS += $(DEBUG)
debug: all

clean:
	rm -rf $(B)

.PHONY: all debug clean

$(B)/test_fwd: $(FWD_O) $(DRV_O) $(B)/test_fwd.o
	$(FC) $(FFLAGS) -o $@ $^
