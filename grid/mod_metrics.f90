module mod_metrics
  implicit none

contains

  !> \brief Calculates the Jacobian of transformation from
  !! (x,y) to (xi,eta) at the center of each real volume
  !! using CDS scheme.
  pure function JpCDS(nx, ny, x, y)
    implicit none
    integer,                   intent(in) :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in) :: ny     !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny), intent(in) :: x      !< x at the ne corner of volume P
    real(8), dimension(nx*ny), intent(in) :: y      !< y at the ne corner of volume P
    real(8), dimension(nx*ny)             :: JpCDS  !< Jacobian at the center of each real volume

    ! Inner variables
    integer :: i, j, np, npw, nps, npsw
    real(8) :: dxk, dyk, dxe, dye

    JpCDS = 0.d0

    do i = 2, nx-1
      do j = 2, ny-1
        np   = nx * (j-1) + i
        nps  = np - nx
        npw  = np - 1
        npsw = nps - 1

        dxk = ( (x(np)+x(nps)) - (x(npw)+x(npsw)) ) / 2.d0
        dyk = ( (y(np)+y(nps)) - (y(npw)+y(npsw)) ) / 2.d0

        dxe = ( (x(np)+x(npw)) - (x(nps)+x(npsw)) ) / 2.d0
        dye = ( (y(np)+y(npw)) - (y(nps)+y(npsw)) ) / 2.d0

        JpCDS(np) = 1.d0 / (dxk * dye - dxe * dyk)

      end do
    end do

  end function


  !> \brief Calculates the alpha metric of transformation from
  !! (x,y) to (xi,eta) at the center of each real volume
  !! using CDS scheme.
  pure function alphapCDS(nx, ny, x, y)
    implicit none
    integer,                   intent(in) :: nx        !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in) :: ny        !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny), intent(in) :: x         !< x at the ne corner of volume P
    real(8), dimension(nx*ny), intent(in) :: y         !< y at the ne corner of volume P
    real(8), dimension(nx*ny)             :: alphapCDS !< alpha at the center of each real volume

    ! Inner variables
    integer :: i, j, np, npw, nps, npsw
    real(8) :: dxe, dye

    alphapCDS = 0.d0

    do i = 2, nx-1
      do j = 2, ny-1
        np   = nx * (j-1) + i
        nps  = np - nx
        npw  = np - 1
        npsw = nps - 1

        dxe = ( (x(np)+x(npw)) - (x(nps)+x(npsw)) ) / 2.d0
        dye = ( (y(np)+y(npw)) - (y(nps)+y(npsw)) ) / 2.d0

        alphapCDS(np) = dxe ** 2.d0 + dye ** 2.d0

      end do
    end do

  end function


  !> \brief Calculates the beta metric of transformation from
  !! (x,y) to (xi,eta) at the center of each real volume
  !! using CDS scheme.
  pure function betapCDS(nx, ny, x, y)
    implicit none
    integer,                   intent(in) :: nx       !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in) :: ny       !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny), intent(in) :: x        !< x at the ne corner of volume P
    real(8), dimension(nx*ny), intent(in) :: y        !< y at the ne corner of volume P
    real(8), dimension(nx*ny)             :: betapCDS !< beta at the center of each real volume

    ! Inner variables
    integer :: i, j, np, npw, nps, npsw
    real(8) :: dxk, dyk, dxe, dye

    betapCDS = 0.d0

    do i = 2, nx-1
      do j = 2, ny-1
        np   = nx * (j-1) + i
        nps  = np - nx
        npw  = np - 1
        npsw = nps - 1

        dxk = ( (x(np)+x(nps)) - (x(npw)+x(npsw)) ) / 2.d0
        dyk = ( (y(np)+y(nps)) - (y(npw)+y(npsw)) ) / 2.d0

        dxe = ( (x(np)+x(npw)) - (x(nps)+x(npsw)) ) / 2.d0
        dye = ( (y(np)+y(npw)) - (y(nps)+y(npsw)) ) / 2.d0

        betapCDS(np) = dxk * dxe + dyk * dye

      end do
    end do

  end function



  !> \brief Calculates the gamma metric of transformation from
  !! (x,y) to (xi,eta) at the center of each real volume
  !! using CDS scheme.
  pure function gammapCDS(nx, ny, x, y)
    implicit none
    integer,                   intent(in) :: nx        !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in) :: ny        !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny), intent(in) :: x         !< x at the ne corner of volume P
    real(8), dimension(nx*ny), intent(in) :: y         !< y at the ne corner of volume P
    real(8), dimension(nx*ny)             :: gammapCDS !< gamma at the center of each real volume

    ! Inner variables
    integer :: i, j, np, npw, nps, npsw
    real(8) :: dxk, dyk

    gammapCDS = 0.d0

    do i = 2, nx-1
      do j = 2, ny-1
        np   = nx * (j-1) + i
        nps  = np - nx
        npw  = np - 1
        npsw = nps - 1

        dxk = ( (x(np)+x(nps)) - (x(npw)+x(npsw)) ) / 2.d0
        dyk = ( (y(np)+y(nps)) - (y(npw)+y(npsw)) ) / 2.d0

        gammapCDS(np) = dxk * dxk + dyk * dyk

      end do
    end do

  end function

end module
