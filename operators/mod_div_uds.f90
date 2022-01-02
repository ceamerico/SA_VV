
!> \brief Provides the discretized div(psi*phi*vec(U)) operator
!! using the UDS discretization scheme.

module mod_div_uds
  implicit none

contains

  !>
  !!  Discretizes the operator signal * div(psi*phi*vec(U)) in
  !!  the numerical coordinate system (xi,eta) using the pure UDS scheme,
  !!  where
  !!    signal = signal to add (1.0d0) or to subtract (-1.d0) this operator
  !!    psi    = an arbitrary function
  !!    phi    = the unknown
  !!    vec(U) = the velocity vector
  !!
  !!  This subroutine requires:
  !!    re, rn, rp: the radius on the east and north faces of each volume
  !!                as well as the radius at the center of each volume
  !!                (cartesian coordinates if r=1 and cylindrical if r=y)
  !!            Jp: the Jacobian of the transformation
  !!           phi: an estimate of phi
  !!    psie, psin: the function psi on the east and north faces of each volume
  !!      Uce, Vcn: the contravariant components of the velocity vector on the east and north faces
  !!
  !!  This function returns:
  !!    a(:,1:9) = coefficients of the discrete operator
  !!    a(:,10)  = source term of the discrete operator
  !!
   pure function div_uds(nx, ny, re, rn, rp, Jp, psie, psin, Uce, Vcn, phi) result(a)
    integer,                      intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny),    intent(in)    :: re     !< Radius of the center of east face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: rn     !< Radius of the center of north face of volume P
    real(8), dimension(nx*ny),    intent(in)    :: rp     !< Radius of the center of volume P
    real(8), dimension(nx*ny),    intent(in)    :: Jp     !< Jacobian of volume P
    real(8), dimension(nx*ny),    intent(in)    :: psie   !< Function psi at east face
    real(8), dimension(nx*ny),    intent(in)    :: psin   !< Function psi at north face
    real(8), dimension(nx*ny),    intent(in)    :: Uce    !< Contravariant velocity U at east face
    real(8), dimension(nx*ny),    intent(in)    :: Vcn    !< Contravariant velocity V at north face
    real(8), dimension(nx*ny),    intent(in)    :: phi    !< Estimate for the unknown phi at the center of the volumes
    real(8), dimension(nx*ny, 6)                :: a      !< Coefficients and source of the linear system

    ! Inner variables
    integer :: i, j, np, nps, npw

    ! Coefficients for UDS scheme
    real(8) :: aw, ae, as, an

    ! Coefficients of the matrix
    real(8) :: Avw, Ave, Avs, Avn

    ! f = psi * r * U or f = psi * r * V on the center of the faces
    real(8) :: fe, fw, fn, fs

    ! Calculating the coefficients and source terms for real volumes
    do i = 2, nx-1
      do j = 2, ny-1

        np  = nx * (j-1) + i
        nps = np - nx
        npw = np - 1

        fw = psie(npw) * re(npw) * Uce(npw)
        fe = psie(np ) * re(np ) * Uce(np )
        fs = psin(nps) * rn(nps) * Vcn(nps)
        fn = psin(np ) * rn(np ) * Vcn(np )

        as = dsign( 0.5d0, Vcn(nps) )
        an = dsign( 0.5d0, Vcn(np)  )
        aw = dsign( 0.5d0, Uce(npw) )
        ae = dsign( 0.5d0, Uce(np)  )

        Avs = - fs * (0.5d0+as) * Jp(np) / rp(np) ! S
        Avw = - fw * (0.5d0+aw) * Jp(np) / rp(np) ! W
        Ave = + fe * (0.5d0-ae) * Jp(np) / rp(np) ! E
        Avn = + fn * (0.5d0-an) * Jp(np) / rp(np) ! N

        a(np,1) = Avs ! S
        a(np,2) = Avw ! W
        a(np,4) = Ave ! E
        a(np,5) = Avn ! N

        a(np,3) = -( Avs + Avw + Ave + Avn )

        ! Source term
        a(np,6) = - phi(np) * (fe-fw+fn-fs) * Jp(np) / rp(np)
      end do

    end do

  end function

end module
