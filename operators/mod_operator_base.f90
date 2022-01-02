!>
!!  \brief mod_operator_base provides the basic '+', '-', '*' and '/'
!!  operations for handling discrete matrices and source
!!  terms from the discretization of partial differential
!!  equations.
!!
!!  The result of the discretization is a linear system A . x = b.
!!  Where A is the matrix of coefficients, x is the unknown vector
!!  and b is the source term.
!!
!!  One important point is that the discrete matrix and source term
!!  are grouped in a single matrix.
!!
!!  Suppose the discretization generates a 5-diagonal matrix A
!!  with coefficientes As, Aw, Ap, Ae, An and source b. Then the
!!  numerical matrix will have 6 diagonals: diagonals 1-5 are reserved
!!  to the coefficients As ... An and the last one is reserved for
!!  the source term b.
!!  Similarly, if the discretization generates a 9-diagonal matrix A
!!  the numerical matrix will have 10 diagonals: diagonals 1-9 are reserved for
!!  the coefficients Asw ... Ane and the last one is reserved for the
!!  the source term b.
!!
!!  Operations of this module where designed to operate matrices A with
!!  1, 5 and 9 diagonals
!!
!!  The .plus. operator can sum matrices of different number of diagonals.
!!  This is performed by summing corresponding coefficients of each
!!  matrix, e.g., Ae=a(:,6) for a 9-diagonal matrix and Ae=a(:,4) for a
!!  5-diagonal matrix.
!!
!!  The .plus. operator also can sum matrices and vectors. Vectors are added
!!  to the source term of the numerical matrix.
!!
!!  The .minus. operator is similar to the .plus. operator.
!!
!!  The * operator multiplies vector and matrices. Each element of the
!!  vector is multiplied by each element of A and b of the discrete equation.
!!
!!  The / operator divides matrices by vectors. Each element of the matrix
!!  is divided by the corresponding element of the vector.
!!

module mod_operator_base
  implicit none

  ! Making all functions private, except otherwise stated
  private

  ! Public methods
  public :: operator(.plus.), operator(.minus.), operator(*), operator(/)!, assignment(.equals.)

  ! Interfaces
  interface operator (.plus.)
    module procedure :: addOp, addOpVec, addVecOp
  end interface

  interface operator (.minus.)
    module procedure subtractOp, subtractOpVec, subtractVecOp
  end interface

  interface operator (*)
    module procedure multiplyOpVec, multiplyVecOp
  end interface

  interface operator (/)
    module procedure divideOpVec
  end interface

contains

  !> \brief Adds the operator op1 and the vector vec
  pure function addOpVec(op1, vec) result(res)
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    real(8), dimension(:),   intent(in) :: vec !< RHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1),size(op1,2)) :: res !< Sum operator

    integer :: ld

    ld = ubound(op1,2)

    res(:,1:ld-1) = op1(:,1:ld-1)

    res(:,ld) = op1(:,ld) + vec

  end function


  !> \brief Adds the vector vec and the operator op1
  pure function addVecOp(vec, op1) result(res)
    real(8), dimension(:),   intent(in) :: vec !< RHS discrete operator
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1),size(op1,2)) :: res !< Sum operator

    integer :: ld

    ld = ubound(op1,2)

    res(:,1:ld-1) = op1(:,1:ld-1)

    res(:,ld) = op1(:,ld) + vec

  end function


  !> \brief Subtracts the vector vec from the operator op1
  pure function subtractOpVec(op1, vec) result(res)
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    real(8), dimension(:),   intent(in) :: vec !< RHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1),size(op1,2)) :: res !< Sum operator

    integer :: ld

    ld = ubound(op1,2)

    res(:,1:ld-1) = op1(:,1:ld-1)

    res(:,ld) = op1(:,ld) - vec

  end function


  !> \brief Subtracts the operator op1 from the vector vec
  pure function subtractVecOp(vec, op1) result(res)
    real(8), dimension(:),   intent(in) :: vec !< RHS discrete operator
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1),size(op1,2)) :: res !< Sum operator

    integer :: ld

    ld = ubound(op1,2)

    res(:,1:ld-1) = -op1(:,1:ld-1)

    res(:,ld) =  vec - op1(:,ld)
  end function


  !> \brief Multiplies the operator op1 and the vector vec
  pure function multiplyOpVec(op1, vec) result(res)
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    real(8), dimension(:),   intent(in) :: vec !< RHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1),size(op1,2)) :: res !< Sum operator

    integer :: ld, k

    ld = ubound(op1,2)

    do k = 1, ld
      res(:,k) = op1(:,k) * vec
    end do

  end function


  !> \brief Multiplies the vector vecop1 and the operator
  pure function multiplyVecOp(vec, op1) result(res)
    real(8), dimension(:),   intent(in) :: vec !< RHS discrete operator
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1),size(op1,2)) :: res !< Sum operator

    integer :: ld, k

    ld = ubound(op1,2)

    do k = 1, ld
      res(:,k) = op1(:,k) * vec
    end do

  end function


  !> \brief Divides the operator op1 by the vector vec
  pure function divideOpVec(op1, vec) result(res)
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    real(8), dimension(:),   intent(in) :: vec !< RHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1),size(op1,2)) :: res !< Sum operator

    integer :: ld, k

    ld = ubound(op1,2)

    do k = 1, ld
      res(:,k) = op1(:,k) / vec
    end do

  end function


  !> \brief Adds op1 and op2 operators
  pure function addOp(op1, op2) result(res)
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    real(8), dimension(:,:), intent(in) :: op2 !< RHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1), max(size(op1,2),size(op2,2))) :: res !< Sum operator

    integer :: ND1, ND2

    ! Getting the number of diagonals + 1 source term of each operator
    ND1 = size(op1,2)
    ND2 = size(op2,2)

    ! If Op1 >= Op2
    if ( ND1 >= ND2 ) then

      select case (ND1)

          ! Op1 9-diagonal + 1 source term
        case (10)

          select case (ND2)

            ! Op2 9-diagonal + 1 source term
            case (10)
              res = op1 + op2

            ! Op2 5-diagonal + 1 source term
            case (6)
              res(:, 1) = op1(:, 1)
              res(:, 2) = op1(:, 2) + op2(:,1)
              res(:, 3) = op1(:, 3)
              res(:, 4) = op1(:, 4) + op2(:,2)
              res(:, 5) = op1(:, 5) + op2(:,3)
              res(:, 6) = op1(:, 6) + op2(:,4)
              res(:, 7) = op1(:, 7)
              res(:, 8) = op1(:, 8) + op2(:,5)
              res(:, 9) = op1(:, 9)
              res(:,10) = op1(:,10) + op2(:,6)

            ! Op2 1-diagonal + 1 source term
            case (2)
              res(:, 1) = op1(:, 1)
              res(:, 2) = op1(:, 2)
              res(:, 3) = op1(:, 3)
              res(:, 4) = op1(:, 4)
              res(:, 5) = op1(:, 5) + op2(:,1)
              res(:, 6) = op1(:, 6)
              res(:, 7) = op1(:, 7)
              res(:, 8) = op1(:, 8)
              res(:, 9) = op1(:, 9)
              res(:,10) = op1(:,10) + op2(:,2)

            ! Impossible to sum, returning a huge number as error indicator
            case default
              res = huge(1.d0)
          end select

        ! Op1 5-diagonal + 1 source term
        case (6)

          select case (ND2)

            ! Op2 5-diagonal + 1 source term
            case (6)
              res = op1 + op2

            ! Op2 1-diagonal + 1 source term
            case (2)
              res(:,1) = op1(:,1)
              res(:,2) = op1(:,2)
              res(:,3) = op1(:,3) + op2(:,1)
              res(:,4) = op1(:,4)
              res(:,5) = op1(:,5)
              res(:,6) = op1(:,6) + op2(:,2)

            ! Impossible to sum, returning a huge number as error indicator
            case default
              res = huge(1.d0)
          end select

        case (2)

          select case (ND2)
            ! Op2 1-diagonal + 1 source term
            case (2)
              res = op1 + op2

            ! Impossible to sum, returning a huge number as error indicator
            case default
              res = huge(1.d0)
          end select

        ! Impossible to sum, returning a huge number as error indicator
        case default
          res = huge(1.d0)
      end select

    ! If Op1 < Op2
    else

      select case (ND2)

        ! op2 9-diagonal + 1 source term
        case (10)

          select case (ND1)

            ! op1 9-diagonal + 1 source term
            case (10)
              res = op2 + op1

            ! op1 5-diagonal + 1 source term
            case (6)
              res(:, 1) = op2(:, 1)
              res(:, 2) = op2(:, 2) + op1(:,1)
              res(:, 3) = op2(:, 3)
              res(:, 4) = op2(:, 4) + op1(:,2)
              res(:, 5) = op2(:, 5) + op1(:,3)
              res(:, 6) = op2(:, 6) + op1(:,4)
              res(:, 7) = op2(:, 7)
              res(:, 8) = op2(:, 8) + op1(:,5)
              res(:, 9) = op2(:, 9)
              res(:,10) = op2(:,10) + op1(:,6)

            ! op1 1-diagonal + 1 source term
            case (2)
              res(:, 1) = op2(:, 1)
              res(:, 2) = op2(:, 2)
              res(:, 3) = op2(:, 3)
              res(:, 4) = op2(:, 4)
              res(:, 5) = op2(:, 5) + op1(:,1)
              res(:, 6) = op2(:, 6)
              res(:, 7) = op2(:, 7)
              res(:, 8) = op2(:, 8)
              res(:, 9) = op2(:, 9)
              res(:,10) = op2(:,10) + op1(:,2)

            ! Impossible to sum, returning a huge number as error indicator
            case default
              res = huge(1.d0)
          end select

        ! op2 5-diagonal + 1 source term
        case (6)

          select case (ND1)

            ! op1 5-diagonal + 1 source term
            case (6)
              res = op2 + op1

            ! op1 1-diagonal + 1 source term
            case (2)
              res(:,1) = op2(:,1)
              res(:,2) = op2(:,2)
              res(:,3) = op2(:,3) + op1(:,1)
              res(:,4) = op2(:,4)
              res(:,5) = op2(:,5)
              res(:,6) = op2(:,6) + op1(:,2)

            ! Impossible to sum, returning a huge number as error indicator
            case default
              res = huge(1.d0)
          end select

        case (2)

          select case (ND1)
            ! op1 1-diagonal + 1 source term
            case (2)
              res = op2 + op1

            ! Impossible to sum, returning a huge number as error indicator
            case default
              res = huge(1.d0)
          end select

        ! Impossible to sum, returning a huge number as error indicator
        case default
          res = huge(1.d0)
      end select

    end if

  end function


  !> \brief Subtracts op2 from op1
  pure function subtractOp(op1, op2) result(res)
    real(8), dimension(:,:), intent(in) :: op1 !< LHS discrete operator
    real(8), dimension(:,:), intent(in) :: op2 !< RHS discrete operator
    ! The sum operator must be of the same size of the greatest operator
    real(8), dimension(size(op1,1), max(size(op1,2),size(op2,2))) :: res !< Sum operator

    res = addOp(op1, -op2)

  end function

end module
