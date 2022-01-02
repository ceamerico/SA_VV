!>
!! \brief mod_class_solver_tdma provides the implementation of TDMA solver
!!        for linear system of 5 and 9 diagonals.
!!
module mod_class_solver_tdma

   use mod_class_solver_abstract
   use mtdma2d5
   use mtdma2d9

   implicit none

   ! Makes everything private, except otherwise stated
   private

    !> \brief Public class that solves a linear system of 5D using TDMA
   type, public, extends(class_solver_abstract) :: class_solver_tdma5d

      integer :: nx   !< Number of volumes in x
      integer :: ny   !< Number of volumes in y
      integer :: nitm !< Maximum number of iterations

   contains

      procedure, public, pass :: init  => init5d  !< Initializer
      procedure, public, pass :: solve => solve5d !< Solves the linear system

   end type

    !> \brief Public class that solves a linear system of 9D using TDMA
   type, public, extends(class_solver_abstract) :: class_solver_tdma9d

      integer :: nx   !< Number of volumes in x
      integer :: ny   !< Number of volumes in y
      integer :: nitm !< Maximum number of iterations

   contains

      procedure, public, pass :: init  => init9d  !< Initializer
      procedure, public, pass :: solve => solve9d !< Solves the linear system

   end type


contains

   !> \brief Initializes the object
   subroutine init5d(this, nx, ny, nitm)
      implicit none
      class(class_solver_tdma5d) :: this !< Reference to this object
      integer,        intent(in) :: nx   !< Number of volumes in x
      integer,        intent(in) :: ny   !< Number of volumes in y
      integer,        intent(in) :: nitm !< Maximum number of iterations

      this%nx   = nx
      this%ny   = ny
      this%nitm = nitm

   end subroutine


   !> \brief Computes the solution of a linear system a . sol = b
   subroutine solve5d(this, a, b, sol)
      implicit none
      class(class_solver_tdma5d)             :: this !< A reference to this object
      real(8), dimension(:,:), intent(in)    :: a    !< Matrix of coefficients
      real(8), dimension(:),   intent(in)    :: b    !< Source vector
      real(8), dimension(:),   intent(inout) :: sol  !< Solution vector

      call tdma2d5(this%nx, this%ny, this%nitm, a, b, sol)

   end subroutine


   !> \brief Initializes the object
   subroutine init9d(this, nx, ny, nitm)
      implicit none
      class(class_solver_tdma9d) :: this !< Reference to this object
      integer,        intent(in) :: nx   !< Number of volumes in x
      integer,        intent(in) :: ny   !< Number of volumes in y
      integer,        intent(in) :: nitm !< Maximum number of iterations

      this%nx   = nx
      this%ny   = ny
      this%nitm = nitm

   end subroutine


   !> \brief Computes the solution of a linear system a . sol = b
   subroutine solve9d(this, a, b, sol)
      implicit none
      class(class_solver_tdma9d)             :: this !< A reference to this object
      real(8), dimension(:,:), intent(in)    :: a    !< Matrix of coefficients
      real(8), dimension(:),   intent(in)    :: b    !< Source vector
      real(8), dimension(:),   intent(inout) :: sol  !< Solution vector

      call tdma2d9(this%nx, this%ny, this%nitm, a, b, sol)

   end subroutine

end module
