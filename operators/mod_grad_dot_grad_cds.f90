
!> \brief Provides the discretized grad(psi).grad(phi) operator
!! using the CDS-2 discretization scheme.

module mod_grad_dot_grad_cds
  implicit none

    contains

  !>
  !!  Discretizes the operator grad(psi) .dot. grad(phi) in
  !!  the numerical coordinate system (xi,eta) using the CDS-2 scheme,
  !!  where
  !!    psi    = an arbitrary function
  !!    phi    = the unknown
  !!
  !!  This subroutine requires:
  !!            Jp: the Jacobian of the transformation at the center of each volume.
  !! alphap, betap,
  !!        gammap: The metrics of the transformation at the center of each volume.
  !!    psie, psin: the function psi on the east and north faces of each volume
  !!
  !!  This subroutine returns:
  !!    b = the source term of the discrete operator.
  !!
  pure function grad_dot_grad_cds(nx, ny, Jp, alphap, betap, gammap, psie, psin, phie, phin) result(b)
    integer,                      intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny),    intent(in)    :: Jp     !< Jacobian of volume P
    real(8), dimension(nx*ny),    intent(in)    :: alphap !< Metric alpha of the center of the volume P
    real(8), dimension(nx*ny),    intent(in)    :: betap  !< Metric beta of the center of the volume P
    real(8), dimension(nx*ny),    intent(in)    :: gammap !< Metric gamma of the center of the volume P
    real(8), dimension(nx*ny),    intent(in)    :: psie   !< Function psi at the center of face east of volume P
    real(8), dimension(nx*ny),    intent(in)    :: psin   !< Function psi at the center of face north of volume P
    real(8), dimension(nx*ny),    intent(in)    :: phie   !< Function phi at the center of face east of volume P
    real(8), dimension(nx*ny),    intent(in)    :: phin   !< Function phi at the center of face north of volume P
    real(8), dimension(nx*ny)                   :: b      !< Source term

    ! Inner variables
    integer :: i, j, np, npw, npe, nps, npn
    real(8) :: dphidxi, dphideta, dpsidxi, dpsideta

    ! Calculating the coefficients and source terms for real volumes
    do i = 2, nx-1
      do j = 2, ny-1

        np   = nx * (j-1) + i
        nps  = np - nx
        npn  = np + nx
        npw  = np - 1
        npe  = np + 1

        dphidxi  = phie(np) - phie(npw)
        dphideta = phin(np) - phin(nps)
        dpsidxi  = psie(np) - psie(npw)
        dpsideta = psin(np) - psin(nps)

        b(np) = - Jp(np)**2.d0 * (&
          & + alphap(np) * dpsidxi  * dphidxi  &
          & + gammap(np) * dpsideta * dphideta &
          & - betap(np)  * ( dphideta * dpsidxi + dpsideta * dphidxi ) &
          & )

      end do
    end do

  end function

end module
