!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!>@brief The module 'adj_face_mod' holds the passive inputs of one advection step on this rank's
!! subdomain (grid metrics, C-grid winds with halo, pressure edges, ak/bk). It is filled by
!! adj_load_mod (adj_load.F90, not seen by TAF) from the files written by prep/make_face_inputs.py.
!>@details The subdomain bounds are run-time values (one face on one process, or a part of a face
!! on 6 or 24 ranks); indices are face-global as in FV3 (1..npx-1 on each face).
module adj_face_mod

  use fv_arrays_mod,   only: fv_grid_type, fv_flags_type, fv_grid_bounds_type
  use fv_mp_mod,       only: ng, group_halo_update_type
  use mpp_domains_mod, only: domain2d

  implicit none
  public

  integer, parameter :: npx = 25, npy = 25, npz = 72, nq = 1
  real, parameter :: dt = 1800.

  integer, save :: is, ie, js, je, isd, ied, jsd, jed, tile

  type(fv_grid_bounds_type), save :: bd
  type(fv_grid_type), target, save :: gridstruct
  type(fv_flags_type), save :: flagstruct
  type(domain2d), save :: domain
  type(group_halo_update_type), save :: q_pack

  !> C-grid winds: x-faces (is:ie+1, jsd:jed), y-faces (isd:ied, js:je+1), with halo
  real, allocatable, save :: cx_e(:,:,:), cy_e(:,:,:), mfx_e(:,:,:), mfy_e(:,:,:)
  real, allocatable, save :: ple0(:,:,:), ple1(:,:,:), pleadv_gchp(:,:,:)
  !> q before the step with halo (single face: halo from the neighbouring faces in the file;
  !! FMS: filled by the halo exchange), and GCHP's q after the step (check value)
  real, allocatable, save :: q_in(:,:,:), q_out_gchp(:,:,:)
  real, save :: ak(npz+1), bk(npz+1), ptop

  !> diagnostics written by the top-level routines (adj_tops.F90)
  real, allocatable, save :: pleadv_out(:,:,:), dpa_out(:,:,:)

  !> mass fixer in top_step: the frozen factor adj_s0 (if adj_freeze_s) or s(q); s(q) is written to
  !! adj_s_last. adj_s0 is only read in the differentiated code, so it stays passive (a TL/AD copy
  !! of an active adj_s0 would carry over between calls).
  logical, save :: adj_freeze_s = .false.
  real, save :: adj_s0 = 1., adj_s_last = 1.

  !> level used by the single-level kernel test (fv_tp_2d)
  integer, save :: k_test = 40

end module adj_face_mod
