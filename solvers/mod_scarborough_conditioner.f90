
!> \brief  Gives a linear system conditioner that satisfies the Scarborough condition

module mod_scarborough_conditioner

  implicit none

contains

  !> \brief Modifies the coefficients and source terms of
  !! the linear system so that the Scarborough condition
  !! is achieved with maximum value sc.
  subroutine scarborough_conditioner(nx, ny, sc, phi, a, b)
    implicit none
    integer,                      intent(in)    :: nx !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny !< Number of real partitions + 2 fictitious in eta
    real(8),                      intent(in)    :: sc !< Scarborough condition
    real(8), dimension(nx*ny),    intent(in)    :: phi!< Estimate for the unknown phi at the center of the volumes
    real(8), dimension(nx*ny, 9), intent(inout) :: a  !< Coefficients of the linear system
    real(8), dimension(nx*ny),    intent(inout) :: b  !< Source of the linear system

    ! Inner variables
    integer :: i, j, np
    real(8) :: sc1
    real(8) :: lambda

    ! Centers of the volumes
    do i = 1, nx

      do j = 1, ny

        np  = nx * (j-1) + i

        ! Current Scarborough constant
        sc1 = sum(abs(a(np,:))) / abs(a(np,5)) - 1.d0

        if ( sc1 > sc ) then

          ! Calculating lambda to get the desired Scarborough constant
          lambda = sc1/sc - 1.d0

          b(np) = b(np) + lambda * a(np,5) * phi(np)

          a(np,5) = a(np,5) * ( 1.d0 + lambda )

        end if

      end do

    end do

  end subroutine

end module
