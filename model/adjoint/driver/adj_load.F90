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

!>@brief The module 'adj_load_mod' fills adj_face_mod for this rank from the files of
!! prep/make_face_inputs.py (not seen by TAF).
!>@details Single face (adj_backend_mod%multiface = .false.): <dir> holds one face; the subdomain
!! is the whole face and all halos (q, winds) come from the file.
!! Six faces: <dir>/f1 ... <dir>/f6; the rank's subdomain is cut out of its face. Grid metrics
!! come from the file (halo included); winds and q get their halos from the FMS exchange, as in
!! offline_tracer_advection, and are compared with the halos in the file.
module adj_load_mod

  use adj_face_mod
  use adj_backend_mod, only: multiface, backend_bounds, backend_cgrid_halo, backend_allmax, backend_rank, &
                             backend_nranks
  use mpp_domains_mod, only: adj_set_cgrid_halo

  implicit none
  private
  public :: load_step, read_r8, dump_field

  integer, parameter :: nf = npx - 1           !< cells per face edge
  integer, parameter :: fsd = 1 - ng, fed = nf + ng
  logical, save :: verbose = .true.

contains

  subroutine read_r8(path, a, nelem)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nelem
    real, intent(out) :: a(nelem)
    integer :: u, ios, fsize
    open(newunit=u, file=trim(path), access='stream', form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
       write(*,*) 'cannot open ', trim(path)
       error stop
    endif
    inquire(unit=u, size=fsize)
    if (fsize /= 8*nelem) then
       write(*,*) trim(path), ': file has ', fsize/8, ' values, expected ', nelem
       error stop
    endif
    read(u) a
    close(u)
  end subroutine read_r8

  !> load one step; a second call replaces the step (the halo checks are printed unless quiet)
  subroutine load_step(dir, quiet)
    character(len=*), intent(in) :: dir
    logical, intent(in), optional :: quiet
    character(len=512) :: fd
    real, allocatable :: g2(:,:), g3(:,:,:)
    real :: wx(1:nf+1, fsd:fed, npz), wy(fsd:fed, 1:nf+1, npz)
    real :: d

    verbose = .true.
    if (present(quiet)) verbose = .not. quiet
    if (allocated(gridstruct%area)) call free_step()
    if (multiface) then
       call backend_bounds(is, ie, js, je, tile)
       write(fd, '(a,a,i0)') trim(dir), '/f', tile
    else
       is = 1; ie = nf; js = 1; je = nf; tile = 1
       fd = dir
    endif
    isd = is - ng; ied = ie + ng; jsd = js - ng; jed = je + ng
    bd%is = is; bd%ie = ie; bd%js = js; bd%je = je
    bd%isd = isd; bd%ied = ied; bd%jsd = jsd; bd%jed = jed
    gridstruct%sw_corner = is == 1  .and. js == 1
    gridstruct%se_corner = ie == nf .and. js == 1
    gridstruct%ne_corner = ie == nf .and. je == nf
    gridstruct%nw_corner = is == 1  .and. je == nf

    allocate(gridstruct%area(isd:ied,jsd:jed), gridstruct%rarea(isd:ied,jsd:jed))
    allocate(gridstruct%dxa(isd:ied,jsd:jed), gridstruct%dya(isd:ied,jsd:jed))
    allocate(gridstruct%dx(isd:ied,jsd:jed+1), gridstruct%dy(isd:ied+1,jsd:jed))
    allocate(gridstruct%rdxc(isd:ied+1,jsd:jed), gridstruct%rdyc(isd:ied,jsd:jed+1))
    allocate(gridstruct%del6_u(isd:ied,jsd:jed+1), gridstruct%del6_v(isd:ied+1,jsd:jed))
    allocate(gridstruct%cosa_s(isd:ied,jsd:jed), gridstruct%rsin2(isd:ied,jsd:jed))
    allocate(gridstruct%sin_sg(isd:ied,jsd:jed,9), gridstruct%area_64(isd:ied,jsd:jed))
    allocate(g2(fsd:fed, fsd:fed))
    call read_r8(trim(fd)//'/area.r8', g2, size(g2));  gridstruct%area = g2(isd:ied, jsd:jed)
    call read_r8(trim(fd)//'/dxa.r8',  g2, size(g2));  gridstruct%dxa  = g2(isd:ied, jsd:jed)
    call read_r8(trim(fd)//'/dya.r8',  g2, size(g2));  gridstruct%dya  = g2(isd:ied, jsd:jed)
    deallocate(g2)
    allocate(g2(fsd:fed, fsd:fed+1))
    call read_r8(trim(fd)//'/dx.r8', g2, size(g2));    gridstruct%dx = g2(isd:ied, jsd:jed+1)
    deallocate(g2)
    allocate(g2(fsd:fed+1, fsd:fed))
    call read_r8(trim(fd)//'/dy.r8', g2, size(g2));    gridstruct%dy = g2(isd:ied+1, jsd:jed)
    deallocate(g2)
    allocate(g3(fsd:fed, fsd:fed, 9))
    call read_r8(trim(fd)//'/sin_sg.r8', g3, size(g3)); gridstruct%sin_sg = g3(isd:ied, jsd:jed, :)
    deallocate(g3)
    gridstruct%rarea   = 1. / gridstruct%area
    gridstruct%area_64 = gridstruct%area
    gridstruct%rdxc = 0.; gridstruct%rdyc = 0.; gridstruct%del6_u = 0.; gridstruct%del6_v = 0.
    gridstruct%cosa_s = 0.; gridstruct%rsin2 = 1.
    gridstruct%da_min = minval(gridstruct%area(is:ie,js:je))

    allocate(cx_e(is:ie+1,jsd:jed,npz), mfx_e(is:ie+1,jsd:jed,npz))
    allocate(cy_e(isd:ied,js:je+1,npz), mfy_e(isd:ied,js:je+1,npz))
    call read_r8(trim(fd)//'/cx.r8',  wx, size(wx)); cx_e  = wx(is:ie+1, jsd:jed, :)
    call read_r8(trim(fd)//'/mfx.r8', wx, size(wx)); mfx_e = wx(is:ie+1, jsd:jed, :)
    call read_r8(trim(fd)//'/cy.r8',  wy, size(wy)); cy_e  = wy(isd:ied, js:je+1, :)
    call read_r8(trim(fd)//'/mfy.r8', wy, size(wy)); mfy_e = wy(isd:ied, js:je+1, :)
    if (multiface) then
       ! halos from the FMS exchange (only the interior west/south values are used as input)
       call cgrid_from_fms(cx_e, cy_e, 'CX/CY')
       call cgrid_from_fms(mfx_e, mfy_e, 'MFX/MFY')
    endif
    call adj_set_cgrid_halo(is, ie, js, je, ng, npz, cx_e, cy_e, mfx_e, mfy_e)

    allocate(g3(1:nf, 1:nf, npz+1))
    allocate(ple0(is:ie,js:je,npz+1), ple1(is:ie,js:je,npz+1), pleadv_gchp(is:ie,js:je,npz+1))
    call read_r8(trim(fd)//'/ple0.r8',   g3, size(g3)); ple0        = g3(is:ie, js:je, :)
    call read_r8(trim(fd)//'/ple1.r8',   g3, size(g3)); ple1        = g3(is:ie, js:je, :)
    call read_r8(trim(fd)//'/pleadv.r8', g3, size(g3)); pleadv_gchp = g3(is:ie, js:je, :)
    deallocate(g3)
    allocate(g3(fsd:fed, fsd:fed, npz))
    allocate(q_in(isd:ied,jsd:jed,npz), q_out_gchp(is:ie,js:je,npz))
    call read_r8(trim(fd)//'/q_in.r8', g3, size(g3)); q_in = g3(isd:ied, jsd:jed, :)
    if (multiface) then
       q_in(:, :, :) = merge(q_in, 0., spread(interior_mask(), 3, npz))
       call adj_halo_q(q_in, ied-isd+1, jed-jsd+1, npz)
       d = maxval(abs(q_in - g3(isd:ied, jsd:jed, :)), mask=spread(face_mask(), 3, npz))
       d = backend_allmax(d)
       if (backend_rank() == 0 .and. verbose) &
            write(*,'(a,es10.2)') 'q halo: FMS exchange vs neighbouring faces in the file, max |diff| =', d
    endif
    deallocate(g3)
    allocate(g3(1:nf, 1:nf, npz))
    call read_r8(trim(fd)//'/q_out_gchp.r8', g3, size(g3)); q_out_gchp = g3(is:ie, js:je, :)
    deallocate(g3)
    call read_r8(trim(fd)//'/ak.r8', ak, npz+1)
    call read_r8(trim(fd)//'/bk.r8', bk, npz+1)
    ptop = ak(1)
    allocate(pleadv_out(is:ie, js:je, npz+1), dpa_out(is:ie, js:je, npz))
  end subroutine load_step

  subroutine free_step()
    deallocate(gridstruct%area, gridstruct%rarea, gridstruct%dxa, gridstruct%dya, gridstruct%dx, gridstruct%dy)
    deallocate(gridstruct%rdxc, gridstruct%rdyc, gridstruct%del6_u, gridstruct%del6_v, gridstruct%cosa_s)
    deallocate(gridstruct%rsin2, gridstruct%sin_sg, gridstruct%area_64)
    deallocate(cx_e, cy_e, mfx_e, mfy_e, ple0, ple1, pleadv_gchp, q_in, q_out_gchp, pleadv_out, dpa_out)
  end subroutine free_step

  !> fill the halo of a C-grid pair from its interior by the FMS calls of offline_tracer_advection,
  !! and report the largest difference to the halo values in the file (prep's geometric rotation)
  subroutine cgrid_from_fms(x, y, name)
    real, intent(inout) :: x(is:ie+1, jsd:jed, npz), y(isd:ied, js:je+1, npz)
    character(len=*), intent(in) :: name
    real :: xf(is:ie+1, jsd:jed, npz), yf(isd:ied, js:je+1, npz), d
    xf = x; yf = y
    call backend_cgrid_halo(x(is:ie, js:je, :), y(is:ie, js:je, :), x, y)
    d = max(maxval(abs(x - xf)), maxval(abs(y - yf))) / max(maxval(abs(xf)), maxval(abs(yf)), tiny(1.))
    d = backend_allmax(d)
    if (backend_rank() == 0 .and. verbose) &
         write(*,'(3a,es10.2)') 'C-grid halo ', name, ': FMS exchange vs file, max rel. diff =', d
  end subroutine cgrid_from_fms

  !> .true. on the compute domain
  function interior_mask() result(m)
    logical :: m(isd:ied, jsd:jed)
    m = .false.
    m(is:ie, js:je) = .true.
  end function interior_mask

  !> .true. on the cells of the data domain that exist (not the cube-corner ghosts)
  function face_mask() result(m)
    logical :: m(isd:ied, jsd:jed)
    integer :: i, j
    do j = jsd, jed
       do i = isd, ied
          m(i,j) = .not. ((i < 1 .or. i > nf) .and. (j < 1 .or. j > nf))
       enddo
    enddo
  end function face_mask

  !> write a rank's piece of a field to <dir>/<tag>_n<ranks>_t<tile>_i<is>_j<js>.r8
  subroutine dump_field(dir, tag, a, nelem)
    character(len=*), intent(in) :: dir, tag
    integer, intent(in) :: nelem
    real, intent(in) :: a(nelem)
    character(len=600) :: path
    integer :: u
    write(path, '(4a,i0,a,i0,a,i0,a,i0,a)') trim(dir), '/', trim(tag), '_n', backend_nranks(), '_t', tile, &
         '_i', is, '_j', js, '.r8'
    open(newunit=u, file=trim(path), access='stream', form='unformatted', status='replace')
    write(u) a
    close(u)
  end subroutine dump_field

end module adj_load_mod
