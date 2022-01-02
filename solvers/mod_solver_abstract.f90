!>
!! \brief Defines an abstract class for creating linear system's solvers
!!
module mod_class_solver_abstract

   implicit none

   ! Makes everything private, except otherwise stated
   private

   !> \brief Public class defining an interface for linear system's solvers
   type, abstract, public :: class_solver_abstract

   contains
      procedure(solve_interface), pass, deferred, public :: solve !< Solves the LS
   end type


   abstract interface

      !> \brief Computes the solution of a linear system a . sol = b
      subroutine solve_interface(this, a, b, sol)
         import class_solver_abstract
         implicit none
         class(class_solver_abstract)           :: this !< A reference to this object
         real(8), dimension(:,:), intent(in)    :: a    !< Matrix of coefficients
         real(8), dimension(:),   intent(in)    :: b    !< Source vector
         real(8), dimension(:),   intent(inout) :: sol  !< Solution vector
      end subroutine

   end interface

end module
