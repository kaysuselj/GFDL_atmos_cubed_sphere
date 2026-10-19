# Discrete TL/AD of the offline tracer advection (model/adjoint), compiled into the fvdycore
# library when ADJOINT is on (fvdycore/CMakeLists.txt: include this file after the srcs list, see
# model/adjoint/gchp/patches/fvdycore_CMakeLists.patch). The sources get fvdycore's compiler
# options (R8 default reals, MAPL_MODE SPMD TIMING). The generated modules have their own names
# (stb_*, adj_*; taf/gen_gchp), so they do not clash with the forward tp_core_mod, fv_tracer2d_mod,
# fv_mapz_mod, fv_fill_mod, fv_mp_mod, mpp_*.
set (fv_adjoint_srcs
  model/adjoint/taf/gen_gchp/stubs_ad_gchp.F90
  model/adjoint/taf/gen_gchp/fv_fill_ad_gchp.F90
  model/adjoint/taf/gen_gchp/tp_core_ad_gchp.F90
  model/adjoint/taf/gen_gchp/fv_mapz_ad_gchp.F90
  model/adjoint/taf/gen_gchp/fv_tracer2d_ad_gchp.F90
  model/adjoint/taf/gen_gchp/adj_face_mod_ad_gchp.F90
  model/adjoint/taf/gen_gchp/adj_tops_ad_gchp.F90
  model/adjoint/comm/adj_comm_fms.F90
  model/adjoint/gchp/adv_adj_discrete.F90
  )
list (APPEND srcs ${fv_adjoint_srcs})
