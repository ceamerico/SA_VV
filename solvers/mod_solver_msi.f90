!>
!! \brief mod_class_solver_msi provides the implementation of MSI solver
!!        for linear system of 5 and 9 diagonals.
!!
module mod_class_solver_msi

   use mod_class_solver_abstract
   use msi2d5
   use msi2d9

   implicit none

   ! Makes everything private, except otherwise stated
   private

    !> \brief Public class that solves a linear system of 5D using MSI
   type, public, extends(class_solver_abstract) :: class_solver_msi5d

      integer :: nx   !< Number of volumes in x
      integer :: ny   !< Number of volumes in y
      integer :: nxy  !< nxy = nx * ny
      integer :: nitm !< Maximum number of iterations
      real(8) :: tol  !< Tolerance for the residual of the linear system
      real(8), allocatable, dimension(:,:) :: dl !< Lower matrix of the MSI method
      real(8), allocatable, dimension(:,:) :: du !< Upper matrix of the MSI method

   contains

      procedure, public, pass :: init  => init5d  !< Initializer
      procedure, public, pass :: solve => solve5d !< Solves the linear system

   end type

    !> \brief Public class that solves a linear system of 9D using MSI
   type, public, extends(class_solver_abstract) :: class_solver_msi9d

      integer :: nx   !< Number of volumes in x
      integer :: ny   !< Number of volumes in y
      integer :: nxy  !< nxy = nx * ny
      integer :: nitm !< Maximum number of iterations
      real(8) :: tol  !< Tolerance for the residual of the linear system
      real(8), allocatable, dimension(:,:) :: dl !< Lower matrix of the MSI method
      real(8), allocatable, dimension(:,:) :: du !< Upper matrix of the MSI method

   contains

      procedure, public, pass :: init  => init9d  !< Initializer
      procedure, public, pass :: solve => solve9d !< Solves the linear system

   end type


contains

   !> \brief Initializes the object
   subroutine init5d(this, nx, ny, nitm, tol)
      implicit none
      class(class_solver_msi5d) :: this !< Reference to this object
      integer,       intent(in) :: nx   !< Number of volumes in x
      integer,       intent(in) :: ny   !< Number of volumes in y
      integer,       intent(in) :: nitm !< Maximum number of iterations
      real(8),       intent(in) :: tol  !< Tolerance for the residual

      this%nx   = nx
      this%ny   = ny
      this%nxy  = nx*ny
      this%nitm = nitm
      this%tol  = tol

      allocate( this%dl(this%nxy,4), this%du(this%nxy,3) )

      this%dl = 0.d0
      this%du = 0.d0

   end subroutine


   !> \brief Computes the solution of a linear system a . sol = b
   subroutine solve5d(this, a, b, sol)
      implicit none
      class(class_solver_msi5d)              :: this !< A reference to this object
      real(8), dimension(:,:), intent(in)    :: a    !< Matrix of coefficients
      real(8), dimension(:),   intent(in)    :: b    !< Source vector
      real(8), dimension(:),   intent(inout) :: sol  !< Solution vector

      ! LU decomposition
      call lu2d5( a, this%dl, this%du, this%nxy, this%nx, this%ny)

      ! Solution of the linear system using MSI
      call fb2d5( a, this%dl, this%du, this%ny, this%nx, this%nxy, sol &
         , this%tol, this%nitm, b)

   end subroutine

   !> \brief Initializes the object
   subroutine init9d(this, nx, ny, nitm, tol)
      implicit none
      class(class_solver_msi9d) :: this !< Reference to this object
      integer,       intent(in) :: nx   !< Number of volumes in x
      integer,       intent(in) :: ny   !< Number of volumes in y
      integer,       intent(in) :: nitm !< Maximum number of iterations
      real(8),       intent(in) :: tol  !< Tolerance for the residual

      this%nx   = nx
      this%ny   = ny
      this%nxy  = nx*ny
      this%nitm = nitm
      this%tol  = tol

      allocate( this%dl(this%nxy,5), this%du(this%nxy,4) )

      this%dl = 0.d0
      this%du = 0.d0

   end subroutine


   !> \brief Computes the solution of a linear system a . sol = b
   subroutine solve9d(this, a, b, sol)
      implicit none
      class(class_solver_msi9d)              :: this !< A reference to this object
      real(8), dimension(:,:), intent(in)    :: a    !< Matrix of coefficients
      real(8), dimension(:),   intent(in)    :: b    !< Source vector
      real(8), dimension(:),   intent(inout) :: sol  !< Solution vector

      ! LU decomposition
      call lu2d9( a, this%dl, this%du, this%nxy, this%nx, this%ny)

      ! Solution of the linear system using MSI
      call fb2d9( a, this%dl, this%du, this%ny, this%nx, this%nxy, sol &
         , this%tol, this%nitm, b)

   end subroutine

end module
