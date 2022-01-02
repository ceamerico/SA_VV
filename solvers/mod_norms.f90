
!> \brief Provides methods to calculate norms

module mod_norms

  implicit none

contains

  !> \brief Calculates the norm L2
  real(8) function normL2(nx, ny, phi)
    implicit none
    integer,                   intent(in) :: nx    !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in) :: ny    !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny), intent(in) :: phi   !< Estimate for the unknown phi at the center of the volumes

    ! Inner variables
    integer :: i, j, np

    normL2 = 0.d0

    ! Centers of the real volumes
    do i = 2, nx-1
      do j = 2, ny-1

        np  = nx * (j-1) + i

        normL2 = normL2 + phi(np)**2

      end do
    end do

    normL2 = sqrt(normL2/(dble(nx-2)*dble(ny-2)))

  end function

end module
