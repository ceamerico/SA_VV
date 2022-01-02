!>
!! \brief "mod_solvers" module provides 5 and 9 diagonal linear systems' solvers.
!!        It also provides procedures for calculation of norms.
!!
module mod_solvers

   use mod_class_ifile
   use mod_class_solver_abstract
   use mod_class_solver_msi
   use mod_class_solver_tdma

   implicit none

contains

   !> \brief Initialize solvers
   subroutine solvers_init(nx, ny, ifile, solver5d, solver9d)
      implicit none
      integer,                               intent(in)  :: nx       !< Number of volumes in x
      integer,                               intent(in)  :: ny       !< Number of volumes in y
      class(class_ifile),                    intent(in)  :: ifile    !< Input file
      class(class_solver_abstract), pointer, intent(out) :: solver5d !< Solver for 5 diag. systems
      class(class_solver_abstract), pointer, intent(out) :: solver9d !< Solver for 9 diag. systems

      ! Inner variables
      character(len=100) :: ksolver ! Kind of solver
      integer            :: nitm_u  ! Maximum number of iterations for solving the linear systems for u, v and T
      integer            :: nitm_p  ! Maximum number of iterations for solving the linear system for p
      real(8)            :: tol_u   ! Tolerance in the MSI for solving the linear systems for u, v and T
      real(8)            :: tol_p   ! Tolerance in the MSI for solving the linear system for p

      ! Selecting the solver
      call ifile%get_value(ksolver, "solver")
      call ifile%get_value(nitm_u,  "nitm_u")
      call ifile%get_value(nitm_p,  "nitm_p")
      call ifile%get_value(tol_u,    "tol_u")
      call ifile%get_value(tol_p,    "tol_p")

      ! Allocating the solvers
      if ( trim(ksolver) == "MSI") then

         allocate( class_solver_msi5d::solver5d )
         allocate( class_solver_msi9d::solver9d )

      else if ( trim(ksolver) == "TDMA") then

         allocate( class_solver_tdma5d::solver5d )
         allocate( class_solver_tdma9d::solver9d )

      else

         write(*,*) "solvers_init:"
         write(*,*) "Unkown option: ", trim(ksolver)
         write(*,*) "Stopping..."
         stop

      end if


      ! Initializing the solver5d
      select type (solver5d)

         type is ( class_solver_msi5d )

            call solver5d%init(nx, ny, nitm_p, tol_p)

         type is ( class_solver_tdma5d )

            call solver5d%init(nx, ny, nitm_p)

      end select


      ! Initializing the solver9d
      select type (solver9d)

         type is ( class_solver_msi9d )

            call solver9d%init(nx, ny, nitm_u, tol_u)

         type is ( class_solver_tdma9d )

            call solver9d%init(nx, ny, nitm_u)

      end select

   end subroutine



   subroutine norm_l1_5d( nx, ny, var, b, a, norm)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: var ! Unknown variable
      real(8), dimension (nx*ny),   intent(in) :: b   ! Source vector of the linear system
      real(8), dimension (nx*ny,5), intent(in) :: a   ! Coefficients of the linear system
      real(8), intent(inout) :: norm

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe

      ! Norm is calculated taking into account only real volumes

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1

            norm = norm + abs(        &
               + a(np,1) * var(nps) &
               + a(np,2) * var(npw) &
               + a(np,3) * var(np ) &
               + a(np,4) * var(npe) &
               + a(np,5) * var(npn) &
               - b(np) )

         end do
      end do

   end subroutine norm_l1_5d


   subroutine norm_l1_9d( nx, ny, var, b, a, norm)
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: var ! Unknown variable
      real(8), dimension (nx*ny),   intent(in) :: b   ! Source vector of the linear system
      real(8), dimension (nx*ny,9), intent(in) :: a   ! Coefficients of the linear system
      real(8), intent(inout) :: norm

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne

      ! Norm is calculated taking into account only real volumes

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npsw = nps - 1
            npse = nps + 1
            npnw = npn - 1
            npne = npn + 1

            norm = norm + abs(        &
               + a(np,1) * var(npsw) &
               + a(np,2) * var(nps ) &
               + a(np,3) * var(npse) &
               + a(np,4) * var(npw ) &
               + a(np,5) * var(np  ) &
               + a(np,6) * var(npe ) &
               + a(np,7) * var(npnw) &
               + a(np,8) * var(npn ) &
               + a(np,9) * var(npne) &
               - b(np) )

         end do
      end do

   end subroutine norm_l1_9d


   !> \brief Calculates the norm L1 of a 9-diagonal linear system A . x = b
   !! If |b| > 10 * epsilon, the relative norm is calculated, otherwise the
   !! absolute norm is applied.
   subroutine norm_l1_9d_relative( nx, ny, var, b, a, norm) ! InOutput: last one
      implicit none
      integer, intent(in) :: nx   ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in) :: ny   ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: var ! Unknown variable
      real(8), dimension (nx*ny),   intent(in) :: b   ! Source vector of the linear system
      real(8), dimension (nx*ny,9), intent(in) :: a   ! Coefficients of the linear system
      real(8), intent(inout) :: norm ! Norm L1

      ! Parameters

      real(8), parameter :: eps = 10.d0 * epsilon(1.d0)

      ! Auxiliary variables
      integer :: i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne

      real(8) :: norma ! Norm L1 of the linear system A . x = b
      real(8) :: normb ! Norm L1 of the source b

      ! Norm is calculated taking into account only real volumes

      norma = 0.d0

      normb = 0.d0

      do j = 2, ny-1

         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npsw = nps - 1
            npse = nps + 1
            npnw = npn - 1
            npne = npn + 1

            norma = norma + abs(        &
               + a(np,1) * var(npsw) &
               + a(np,2) * var(nps ) &
               + a(np,3) * var(npse) &
               + a(np,4) * var(npw ) &
               + a(np,5) * var(np  ) &
               + a(np,6) * var(npe ) &
               + a(np,7) * var(npnw) &
               + a(np,8) * var(npn ) &
               + a(np,9) * var(npne) &
               - b(np) )

            normb = normb + abs( b(np) )

         end do

      end do

      ! Summing the input norm with the calculated norm

      if ( normb > eps ) then ! Uses relative norm

         norm = norm + norma / normb

      else

         norm = norm + norma ! Uses absolute norm

      end if

   end subroutine norm_l1_9d_relative



   !> \brief Calculates the norm l1 of the source B of linear system Ax=B, taking into account only real volumes.
   real(8) function norm_l1_b( nx, ny, b)
      implicit none
      integer, intent(in) :: nx   !< Number of volumes in the csi direction (real+fictitious)
      integer, intent(in) :: ny   !< Number of volumes in the eta direction (real+fictitious)
      real(8), dimension (nx*ny),   intent(in) :: b   !< Source of the linear system

      ! Auxiliary variables
      integer :: i, j, np
      real(8) :: norm

      ! Norm is calculated taking into account only real volumes

      norm = 0.d0

      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i

            norm = norm + abs( b(np) )

         end do
      end do

      norm_l1_b = norm

   end function

end module
