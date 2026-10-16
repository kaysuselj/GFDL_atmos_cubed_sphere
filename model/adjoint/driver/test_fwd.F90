!> Forward check of the single-face build: one step of offline_tracer_advection with the
!! operational scheme (hord_tr 12, kord_tr 8) against the GCHP trajectory (q after the step, PLEadv).
!! usage: test_fwd <face data dir> [hord kord]
program test_fwd
  use adj_face_mod
  use mpp_mod,      only: adj_sum_mode, adj_sum_saved
  use fv_mp_mod,    only: adj_qhalo_fill, adj_qhalo
  implicit none
  character(len=512) :: dir, arg
  real :: q(is:ie, js:je, npz, nq), d(is:ie, js:je, npz), q2(is:ie, js:je, npz, nq)
  real :: pleadv_ref(is:ie, js:je, npz+1)
  integer :: k

  call get_command_argument(1, dir)
  call load_face(trim(dir))
  flagstruct%hord_tr = 12; flagstruct%kord_tr = 8
  if (command_argument_count() >= 3) then
     call get_command_argument(2, arg); read(arg, *) flagstruct%hord_tr
     call get_command_argument(3, arg); read(arg, *) flagstruct%kord_tr
  endif
  q(:, :, :, 1) = q_in(is:ie, js:je, :)
  allocate(adj_qhalo(ied-isd+1, jed-jsd+1, npz, nq))
  adj_qhalo(:, :, :, 1) = q_in
  adj_qhalo_fill = .true.
  adj_sum_mode = 1
  q2 = q
  call top_offline(q)
  pleadv_ref = pleadv_out
  call top_step(q2)
  write(*,'(a,2es10.2,l3)') 'top_step vs offline_tracer_advection: max|dq|, max|dpleadv|, bit-identical:', &
       maxval(abs(q2-q)), maxval(abs(pleadv_out-pleadv_ref)), all(q2 == q) .and. all(pleadv_out == pleadv_ref)
  write(*,'(a,2i4)') 'hord_tr, kord_tr =', flagstruct%hord_tr, flagstruct%kord_tr
  write(*,'(a,es24.16)') 'mass-fixer scaling factor (this face) =', adj_sum_saved(1)/adj_sum_saved(2)
  ! GCHP's mass fixer uses global sums (s ~ 1), this one face only: compare q / s_face.
  ! k = npz is left out: GCHP's output there includes the surface fluxes.
  q = q / (adj_sum_saved(1)/adj_sum_saved(2))
  d = q(:, :, :, 1) - q_out_gchp
  write(*,'(a,es12.4)') 'k < npz, s divided out: max |q - q_gchp| / max q =', &
       maxval(abs(d(:,:,1:npz-1))) / maxval(q_out_gchp)
  write(*,'(a,es12.4)') '                          max |q_gchp - q_in| / max q =', &
       maxval(abs(q_out_gchp(:,:,1:npz-1) - q_in(is:ie,js:je,1:npz-1))) / maxval(q_out_gchp)
  write(*,'(a,es12.4)') 'max |q - q_gchp| / max q          =', maxval(abs(d)) / maxval(q_out_gchp)
  write(*,'(a,es12.4)') 'max |q - q_in| / max q (change)   =', maxval(abs(q(:,:,:,1) - q_in(is:ie,js:je,:))) / maxval(q_out_gchp)
  write(*,'(a,es12.4)') 'max |q_gchp - q_in| / max q       =', maxval(abs(q_out_gchp - q_in(is:ie,js:je,:))) / maxval(q_out_gchp)
  write(*,'(a,es12.4,a)') 'max |pleadv - PLEadv_gchp|        =', maxval(abs(pleadv_out - pleadv_gchp)), ' Pa'
  do k = 1, npz, 10
     write(*,'(a,i3,2es12.4)') '  k, rms(q-q_gchp), rms(q_gchp-q_in):', k, &
          sqrt(sum(d(:,:,k)**2)/(n*n)), sqrt(sum((q_out_gchp(:,:,k)-q_in(is:ie,js:je,k))**2)/(n*n))
  enddo
  ! the step with zero q halo (as in test_adv), for the comparison with the TAF-generated forward
  adj_qhalo_fill = .false.
  q2(:, :, :, 1) = q_in(is:ie, js:je, :)
  call top_step(q2)
  call dump_step(trim(dir), 'orig', q2)
end program test_fwd
