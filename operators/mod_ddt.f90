
!> \brief Provides the discrete d(psi*phi)/dt operator
module mod_ddt
  implicit none

    contains

  !>
  !!  Discretizes the operator signal* d(psi*phi)/dt in
  !!  the numerical coordinate system (xi,eta,t)
  !!  where:
  !!    signal = signal to add (1.0d0) or to subtract (-1.d0) this operator
  !!    psi    = an arbitrary function
  !!    phi    = the unknown
  !!
  !!  This subroutine requires:
  !!    nx   : Number of real partitions + 2 fictitious in xi
  !!    ny   : Number of real partitions + 2 fictitious in eta
  !!    psi  : the function psi on the center of each volume
  !!    psia : the function psi on the center of each volume at previous time step
  !!    phia : the solution phi on the center of each volume at previous time step
  !!    dt   : timestep
  !!
  !!  This function returns:
  !!    a(1) = coefficients of the discrete operator
  !!    a(2) = source term of the discrete operator

  pure function ddt(nx, ny, dt, psia, phia, psi) result(a)
    integer,                      intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta
    real(8),                      intent(in)    :: dt     !< Time step
    real(8), dimension(nx*ny),    intent(in)    :: psia   !< Function psi at the center of the volume in the previous temporal iteration
    real(8), dimension(nx*ny),    intent(in)    :: phia   !< Function phi at the center of the volume in the previous temporal iteration
    real(8), dimension(nx*ny),    intent(in)    :: psi    !< Function psi at the center of the volume
    real(8), dimension(nx*ny, 2)                :: a      !< Coefficients and source of the linear system

    ! As a reminder, the discretized equation is: [ (psi * phi)_P - (psia * phia)_P ] / dt,
    a(:,1) = psi / dt

    ! That's why [(psia_P * phia_P / dt ] is added to b, when written as: ap * phi_p + sum_ab * phi_ab = bp
    a(:,2) = phia * psia / dt

  end function

end module
