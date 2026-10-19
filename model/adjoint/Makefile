# TL/AD test builds of the FV3 tracer advection (CPU, r8). The forward files are compiled unchanged
# from ../ against the stubs in stubs/.
#   make             build/test_fwd    original forward, one face (no MPI)
#                    build_ad/test_adv TAF TL/AD, one face (no MPI)
#   make mpi         build_mpi/test_adv, build_mpi/test_halo  TAF TL/AD on six faces (FMS, MPI)
#   make gchp        build_gchp/test_gchp_adapter: the GCHP adapter (gchp/adv_adj_discrete.F90) in a
#                    GCHP-like build (real fv_arrays.F90 and fv_mp_mod.F90, and the generated code as
#                    GCHP compiles it: taf/gen_gchp, all modules renamed, stb_* and adj_*)
#   taf/run_taf.sh   regenerate the TAF TL/AD code in taf/gen (needs the FastOpt server)
# FMS: fms/build_fms.py builds the needed part of geoschem/FMS 259759d1 into $(FMSB).
FC      = gfortran
MPIFC   = mpif90
# CPP definitions of GCHP's fvdycore build (fvdycore/CMakeLists.txt; hydrostatic: no MOIST_CAPPA/USE_COND)
CPPDEFS = -DMAPL_MODE -DSPMD -DTIMING
FFLAGS  = -fdefault-real-8 -fdefault-double-8 -O2 -g -fbacktrace -cpp -ffree-line-length-none $(CPPDEFS)
DEBUG   = -O0 -fcheck=bounds,pointer -ffpe-trap=invalid,zero,overflow
FMSB    ?= ../../../build_fms
FMSSRC  ?= ../../../FMS
NETCDF  = $(shell nf-config --flibs) -L$(shell brew --prefix netcdf)/lib
# fv_mp_mod.F90 calls MPI through mpif.h with mixed argument types (fine with Intel, an error in gfortran >= 10)
COMMFLAGS = -fallow-argument-mismatch

B       = build
FWD_O   = $(B)/stubs.o $(B)/adj_comm_single.o $(B)/tp_core.o $(B)/fv_fill.o $(B)/fv_mapz.o $(B)/fv_tracer2d.o
DRV_O   = $(B)/adj_face_mod.o $(B)/adj_load.o $(B)/adj_tops.o

# TAF output (taf/gen): modules with the forward, TL and AD routines; separate build directories,
# because the generated modules have the same names as the forward ones
GEN     = stubs fv_fill tp_core fv_mapz fv_tracer2d adj_face_mod adj_tops
BA      = build_ad
GEN_O   = $(patsubst %,$(BA)/%_ad.o,$(GEN))

# six faces: the stub modules are renamed (stb_*), so that the real FMS and fv_mp_mod can be linked
BM      = build_mpi
STUBMODS = platform_mod constants_mod field_manager_mod tracer_manager_mod fv_timing_mod mpp_mod \
           mpp_domains_mod fv_mp_mod fv_arrays_mod boundary_mod fv_grid_utils_mod fv_cmp_mod
RENAME  = $(foreach m,$(STUBMODS),-D$(m)=stb_$(m))
GENM_O  = $(patsubst %,$(BM)/taf/%_ad.o,$(GEN))
COMM_O  = $(BM)/comm/fv_arrays_fms.o $(BM)/comm/fv_mp_mod.o $(BM)/comm/adj_comm_fms.o

# GCHP-like build: the generated forward modules get other names too, because GCHP has the real ones
BG      = build_gchp
FWDMODS = tp_core_mod fv_tracer2d_mod fv_mapz_mod fv_fill_mod
RENAME_ALL = $(RENAME) $(foreach m,$(FWDMODS),-D$(m)=adj_$(m))
GENG_O  = $(patsubst %,$(BG)/taf/%_ad.o,$(GEN))
COMMG_O = $(BG)/comm/fv_arrays.o $(BG)/comm/fv_mp_mod.o $(BG)/comm/adj_comm_fms.o

all: $(B)/test_fwd $(BA)/test_adv
mpi: $(BM)/test_adv $(BM)/test_halo $(BM)/test_chain

# ---------------------------------------------------------------- build: original forward, one face
$(B):
	mkdir -p $(B)
$(B)/stubs.o: stubs/stubs.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/adj_comm_single.o: comm/adj_comm_single.F90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/%.o: ../%.F90 $(B)/stubs.o
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/fv_mapz.o: $(B)/fv_fill.o
$(B)/fv_tracer2d.o: $(B)/tp_core.o $(B)/fv_mapz.o $(B)/fv_fill.o
$(B)/%.o: driver/%.F90 $(FWD_O)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@
$(B)/adj_load.o: $(B)/adj_face_mod.o
$(B)/adj_tops.o: $(B)/adj_face_mod.o
$(B)/test_fwd.o: $(DRV_O)
$(B)/test_fwd: $(FWD_O) $(DRV_O) $(B)/test_fwd.o
	$(FC) $(FFLAGS) -o $@ $^

# ---------------------------------------------------------------- build_ad: TAF code, one face
$(BA):
	mkdir -p $(BA)
$(BA)/%_ad.o: taf/gen/%_ad.f90 | $(BA)
	$(FC) $(FFLAGS) -J$(BA) -c $< -o $@
$(BA)/adj_comm_single.o: comm/adj_comm_single.F90 | $(BA)
	$(FC) $(FFLAGS) -J$(BA) -c $< -o $@
$(BA)/fv_fill_ad.o $(BA)/tp_core_ad.o $(BA)/adj_face_mod_ad.o: $(BA)/stubs_ad.o
$(BA)/fv_mapz_ad.o: $(BA)/fv_fill_ad.o
$(BA)/fv_tracer2d_ad.o: $(BA)/tp_core_ad.o $(BA)/fv_mapz_ad.o
$(BA)/adj_tops_ad.o: $(BA)/fv_tracer2d_ad.o $(BA)/adj_face_mod_ad.o
$(BA)/adj_load.o: driver/adj_load.F90 $(GEN_O) $(BA)/adj_comm_single.o
	$(FC) $(FFLAGS) -J$(BA) -c $< -o $@
$(BA)/test_adv.o: driver/test_adv.F90 $(BA)/adj_load.o
	$(FC) $(FFLAGS) -J$(BA) -c $< -o $@
$(BA)/test_adv: $(GEN_O) $(BA)/adj_comm_single.o $(BA)/adj_load.o $(BA)/test_adv.o
	$(FC) $(FFLAGS) -o $@ $^

# ---------------------------------------------------------------- build_mpi: TAF code, six faces
$(BM)/comm $(BM)/taf:
	mkdir -p $@
$(BM)/comm/fv_arrays_fms.o: comm/fv_arrays_fms.F90 | $(BM)/comm
	$(MPIFC) $(FFLAGS) -I$(FMSB) -J$(BM)/comm -c $< -o $@
$(BM)/comm/fv_mp_mod.o: ../../tools/fv_mp_mod.F90 $(BM)/comm/fv_arrays_fms.o
	$(MPIFC) $(FFLAGS) $(COMMFLAGS) -I$(FMSB) -J$(BM)/comm -c $< -o $@
$(BM)/comm/adj_comm_fms.o: comm/adj_comm_fms.F90 $(BM)/comm/fv_mp_mod.o
	$(MPIFC) $(FFLAGS) $(COMMFLAGS) -I$(FMSB) -J$(BM)/comm -c $< -o $@
$(BM)/taf/%_ad.o: taf/gen/%_ad.f90 | $(BM)/taf
	$(MPIFC) $(FFLAGS) $(RENAME) -J$(BM)/taf -c $< -o $@
$(BM)/taf/fv_fill_ad.o $(BM)/taf/tp_core_ad.o $(BM)/taf/adj_face_mod_ad.o: $(BM)/taf/stubs_ad.o
$(BM)/taf/fv_mapz_ad.o: $(BM)/taf/fv_fill_ad.o
$(BM)/taf/fv_tracer2d_ad.o: $(BM)/taf/tp_core_ad.o $(BM)/taf/fv_mapz_ad.o
$(BM)/taf/adj_tops_ad.o: $(BM)/taf/fv_tracer2d_ad.o $(BM)/taf/adj_face_mod_ad.o
$(BM)/taf/adj_load.o: driver/adj_load.F90 $(GENM_O) $(COMM_O)
	$(MPIFC) $(FFLAGS) $(RENAME) -I$(BM)/comm -J$(BM)/taf -c $< -o $@
$(BM)/taf/test_adv.o: driver/test_adv.F90 $(BM)/taf/adj_load.o
	$(MPIFC) $(FFLAGS) $(RENAME) -I$(BM)/comm -J$(BM)/taf -c $< -o $@
$(BM)/test_adv: $(GENM_O) $(BM)/taf/adj_load.o $(BM)/taf/test_adv.o $(COMM_O)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(FMSB)/libfms.a $(NETCDF)
$(BM)/taf/test_chain.o: driver/test_chain.F90 $(BM)/taf/adj_load.o
	$(MPIFC) $(FFLAGS) $(RENAME) -I$(BM)/comm -J$(BM)/taf -c $< -o $@
$(BM)/test_chain: $(GENM_O) $(BM)/taf/adj_load.o $(BM)/taf/test_chain.o $(COMM_O)
	$(MPIFC) $(FFLAGS) -o $@ $^ $(FMSB)/libfms.a $(NETCDF)
$(BM)/comm/test_halo.o: driver/test_halo.F90 $(COMM_O)
	$(MPIFC) $(FFLAGS) -I$(BM)/comm -J$(BM)/comm -c $< -o $@
$(BM)/test_halo: $(COMM_O) $(BM)/comm/test_halo.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(FMSB)/libfms.a $(NETCDF)

# ---------------------------------------------------------------- build_gchp: adapter, GCHP-like
$(BG)/comm $(BG)/taf:
	mkdir -p $@
$(BG)/comm/fv_arrays.o: ../fv_arrays.F90 | $(BG)/comm
	$(MPIFC) $(FFLAGS) -I$(FMSSRC)/include -I$(FMSB) -J$(BG)/comm -c $< -o $@
$(BG)/comm/fv_mp_mod.o: ../../tools/fv_mp_mod.F90 $(BG)/comm/fv_arrays.o
	$(MPIFC) $(FFLAGS) $(COMMFLAGS) -I$(FMSB) -J$(BG)/comm -c $< -o $@
$(BG)/comm/adj_comm_fms.o: comm/adj_comm_fms.F90 $(BG)/comm/fv_mp_mod.o
	$(MPIFC) $(FFLAGS) $(COMMFLAGS) -I$(FMSB) -J$(BG)/comm -c $< -o $@
$(BG)/taf/%_ad.o: taf/gen_gchp/%_ad_gchp.F90 | $(BG)/taf
	$(MPIFC) $(FFLAGS) -J$(BG)/taf -c $< -o $@
$(BG)/taf/fv_fill_ad.o $(BG)/taf/tp_core_ad.o $(BG)/taf/adj_face_mod_ad.o: $(BG)/taf/stubs_ad.o
$(BG)/taf/fv_mapz_ad.o: $(BG)/taf/fv_fill_ad.o
$(BG)/taf/fv_tracer2d_ad.o: $(BG)/taf/tp_core_ad.o $(BG)/taf/fv_mapz_ad.o
$(BG)/taf/adj_tops_ad.o: $(BG)/taf/fv_tracer2d_ad.o $(BG)/taf/adj_face_mod_ad.o
$(BG)/taf/adj_load.o: driver/adj_load.F90 $(GENG_O) $(COMMG_O)
	$(MPIFC) $(FFLAGS) $(RENAME_ALL) -I$(BG)/comm -J$(BG)/taf -c $< -o $@
$(BG)/comm/adv_adj_discrete.o: gchp/adv_adj_discrete.F90 $(GENG_O) $(COMMG_O)
	$(MPIFC) $(FFLAGS) -I$(FMSB) -I$(BG)/taf -J$(BG)/comm -c $< -o $@
$(BG)/comm/test_gchp_adapter.o: driver/test_gchp_adapter.F90 $(BG)/comm/adv_adj_discrete.o $(BG)/taf/adj_load.o
	$(MPIFC) $(FFLAGS) -I$(FMSB) -I$(BG)/taf -J$(BG)/comm -c $< -o $@
$(BG)/test_gchp_adapter: $(GENG_O) $(BG)/taf/adj_load.o $(COMMG_O) $(BG)/comm/adv_adj_discrete.o $(BG)/comm/test_gchp_adapter.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(FMSB)/libfms.a $(NETCDF)
gchp: $(BG)/test_gchp_adapter

debug: FFLAGS += $(DEBUG)
debug: all

clean:
	rm -rf $(B) $(BA) $(BM) $(BG)

.PHONY: all mpi gchp debug clean
