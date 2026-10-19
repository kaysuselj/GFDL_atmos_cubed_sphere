!> Dot-product test of the q halo update of the FMS back-end alone: H = GCHP's group halo update
!! (fv_mp_mod start/complete_group_halo_update) on the data-domain arrays of all ranks, H^T = FMS
!! mpp_update_domains_ad. <H x, y> = <x, H^T y>, sums over the data domains of all ranks.
!! usage: mpirun -n 6*layout**2 test_halo [layout]
program test_halo
  use adj_backend_mod, only: backend_init, backend_finalize, backend_rank, backend_nranks, &
                             backend_bounds, backend_allsum, halo_q, halo_q_ad
  implicit none
  integer, parameter :: ng = 3, nk = 72
  integer :: layout, is, ie, js, je, tile, i, j, k
  character(len=32) :: arg
  real, allocatable :: x(:,:,:), y(:,:,:), hx(:,:,:), hty(:,:,:)
  real :: lhs, rhs, a

  layout = 1
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg); read(arg, *) layout
  endif
  call backend_init(layout)
  call backend_bounds(is, ie, js, je, tile)
  allocate(x(is-ng:ie+ng, js-ng:je+ng, nk), y(is-ng:ie+ng, js-ng:je+ng, nk))
  do k = 1, nk
     do j = js-ng, je+ng
        do i = is-ng, ie+ng
           a = sin(12.9898*i + 78.233*j + 37.719*k + 4.581*tile)*43758.5453
           x(i,j,k) = a - floor(a) - 0.5
           a = sin(17.31*i + 51.77*j + 11.13*k + 9.27*tile)*24634.6345
           y(i,j,k) = a - floor(a) - 0.5
        enddo
     enddo
  enddo
  if (backend_rank() == 0) write(*,*) 'test_halo: forward halo update'
  flush(6)
  hx = x
  call halo_q(hx)
  if (backend_rank() == 0) write(*,*) 'test_halo: adjoint halo update'
  flush(6)
  hty = y
  call halo_q_ad(hty)
  if (backend_rank() == 0) write(*,*) 'test_halo: sums'
  flush(6)
  lhs = backend_allsum(sum(hx*y))
  rhs = backend_allsum(sum(x*hty))
  if (backend_rank() == 0) then
     write(*,'(a,i3,a)') 'q halo update, ', backend_nranks(), ' ranks'
     write(*,'(a,2es24.16,es10.2)') '  <H x, y>, <x, H^T y>, rel. diff:', lhs, rhs, abs(lhs-rhs)/abs(lhs)
     write(*,'(a,es10.2)') '  max |H x - x| on the compute domain (must be 0):', maxval(abs(hx(is:ie,js:je,:) - x(is:ie,js:je,:)))
  endif
  call backend_finalize()
end program test_halo
