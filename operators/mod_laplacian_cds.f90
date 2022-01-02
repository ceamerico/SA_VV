
!> \brief Provides the discrete div(psi*grad(phi)) operator
!!        using the CDS-2 discretization scheme.

module mod_laplacian_cds
  implicit none

    contains

  !>
  !!  Discretizes the operator div(psi*grad(phi)) in
  !!  the numerical coordinate system (xi,eta) using the CDS-2 scheme,
  !!  where
  !!    psi    = an arbitrary function
  !!    phi    = the unknown
  !!
  !!  This subroutine requires:
  !!    re, rn, rp: the radius on the east and north faces of each volume
  !!                as well as the radius at the center of each volume
  !!                (cartesian coordinates if r=1 and cylindrical if r=y).
  !!    Je, Jn, Jp: the Jacobian of the transformation on the east and north faces of each volume
  !!                as well as the Jacobian at the center of each volume.
  !! alphae, betae: Functions alpha and beta, grid metrics dependents. They are calculated on east faces of each volume.
  !!  betan, gaman: Functions beta and gama, grid metrics dependents. They are calculated on north faces of each volume.
  !!    psie, psin: the function psi on the east and north faces of each volume
  !!
  !!  This function returns the matrix 'a', where
  !!    a(:,1:9) = coefficients of the matrix of the discrete operator.
  !!    a(:,10)  = source of the discrete operator.
  !!
  pure function laplacian_cds(nx, ny, re, rn, rp, Jp, Je, Jn, alphae, betan, betae, gaman, psie, psin) result(a)
    integer,                      intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny),    intent(in)    :: re     !< Radius of the center of east face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: rn     !< Radius of the center of north face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: rp     !< Radius of the center of volume P
    real(8), dimension(nx*ny),    intent(in)    :: Jp     !< Jacobian of volume P
    real(8), dimension(nx*ny),    intent(in)    :: Jn     !< Jacobian of the center of north face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: Je     !< Jacobian of the center of east face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: alphae !< Metric alpha of the center of east face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: betan  !< Metric beta of the center of north face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: betae  !< Metric beta of the center of east face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: gaman  !< Metric gama of the center of north face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: psie   !< Function psi at east face
    real(8), dimension(nx*ny),    intent(in)    :: psin   !< Function psi at north face
    real(8), dimension(nx*ny, 10)               :: a      !< Coefficients and source of the linear system

    ! Inner variables
    integer :: i, j, np, nps, npw

    ! f = psi * r * J * gama or f = psi * r * J * alpha or f = psi * r * J * beta on the center of the faces
    real(8) :: f1e, f2e, f1w, f2w, f1n, f2n, f1s, f2s

    ! Calculating the coefficients and source terms for real volumes
    do i = 2, nx-1
      do j = 2, ny-1

        np  = nx * (j-1) + i
        nps = np - nx
        npw = np - 1

        f1w = re(npw) * psie(npw) * Je(npw) * alphae(npw)
        f2w = re(npw) * psie(npw) * Je(npw) * betae(npw)
        f1e = re(np ) * psie(np ) * Je(np ) * alphae(np )
        f2e = re(np ) * psie(np ) * Je(np ) * betae(np )
        f1s = rn(nps) * psin(nps) * Jn(nps) * gaman(nps)
        f2s = rn(nps) * psin(nps) * Jn(nps) * betan(nps)
        f1n = rn(np ) * psin(np ) * Jn(np ) * gaman(np )
        f2n = rn(np ) * psin(np ) * Jn(np ) * betan(np )

        a(np,2) = ( - f2w / 4.d0 + f1s + f2e / 4.d0 ) * Jp(np) /  rp(np) ! S
        a(np,4) = ( - f2s / 4.d0 + f1w + f2n / 4.d0 ) * Jp(np) /  rp(np) ! W
        a(np,6) = ( + f2s / 4.d0 + f1e - f2n / 4.d0 ) * Jp(np) /  rp(np) ! E
        a(np,8) = ( + f2w / 4.d0 + f1n - f2e / 4.d0 ) * Jp(np) /  rp(np) ! N

        a(np,1) = ( -f2w - f2s ) * Jp(np) / ( 4.d0 * rp(np) ) ! SW
        a(np,3) = ( +f2e + f2s ) * Jp(np) / ( 4.d0 * rp(np) ) ! SE
        a(np,7) = ( +f2n + f2w ) * Jp(np) / ( 4.d0 * rp(np) ) ! NW
        a(np,9) = ( -f2n - f2e ) * Jp(np) / ( 4.d0 * rp(np) ) ! NE

        a(np,5) = - ( a(np,2) + a(np,4) + a(np,6) + a(np,8) )  ! this ap can also be computed as: a(np,5) = a(np,5) + signal *  (-f1w - f1s - f1n - f1e ) * Jp(np) / rp(np)

      end do

    end do

    ! Source term
    a(:,10) = 0.d0

  end function

end module
