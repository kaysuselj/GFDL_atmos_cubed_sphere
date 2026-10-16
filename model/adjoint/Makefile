# Single-face TL/AD test build (CPU, r8, no FMS/MPI). The forward files are compiled unchanged
# from ../ against the stubs in stubs/.
#   make            forward library + test driver
#   taf/run_taf.sh  regenerate the TAF TL/AD code in taf/gen (needs the FastOpt server)
FC      = gfortran
FFLAGS  = -fdefault-real-8 -fdefault-double-8 -O2 -g -fbacktrace -cpp -ffree-line-length-none
DEBUG   = -O0 -fcheck=bounds,pointer -ffpe-trap=invalid,zero,overflow
B       = build

FWD     = ../tp_core.F90 ../fv_fill.F90 ../fv_mapz.F90 ../fv_tracer2d.F90
FWD_O   = $(B)/stubs.o $(B)/tp_core.o $(B)/fv_fill.o $(B)/fv_mapz.o $(B)/fv_tracer2d.o

DRV_O   = $(B)/adj_face_mod.o $(B)/adj_tops.o

# TAF output (taf/gen): modules with the forward, TL and AD routines; separate build directory,
# because the generated modules have the same names as the forward ones
BA      = build_ad
GEN     = stubs fv_fill tp_core fv_mapz fv_tracer2d adj_face_mod adj_tops
GEN_O   = $(patsubst %,$(BA)/%_ad.o,$(GEN))

all: $(B)/test_fwd $(BA)/test_adv

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
$(BA):
	mkdir -p $(BA)
$(BA)/%_ad.o: taf/gen/%_ad.f90 | $(BA)
	$(FC) $(FFLAGS) -J$(BA) -c $< -o $@
$(BA)/fv_fill_ad.o: $(BA)/stubs_ad.o
$(BA)/tp_core_ad.o: $(BA)/stubs_ad.o
$(BA)/fv_mapz_ad.o: $(BA)/fv_fill_ad.o
$(BA)/fv_tracer2d_ad.o: $(BA)/tp_core_ad.o $(BA)/fv_mapz_ad.o
$(BA)/adj_face_mod_ad.o: $(BA)/stubs_ad.o
$(BA)/adj_tops_ad.o: $(BA)/fv_tracer2d_ad.o $(BA)/adj_face_mod_ad.o
$(BA)/test_adv.o: driver/test_adv.F90 $(GEN_O)
	$(FC) $(FFLAGS) -J$(BA) -c $< -o $@
$(BA)/test_adv: $(GEN_O) $(BA)/test_adv.o
	$(FC) $(FFLAGS) -o $@ $^

debug: FFLAGS += $(DEBUG)
debug: all

clean:
	rm -rf $(B) $(BA)

.PHONY: all debug clean

$(B)/test_fwd: $(FWD_O) $(DRV_O) $(B)/test_fwd.o
	$(FC) $(FFLAGS) -o $@ $^
