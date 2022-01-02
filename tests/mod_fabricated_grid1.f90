module mod_fabricated_grid1
  implicit none

  real(8), parameter :: pi = acos(-1.d0)

contains

  real(8) pure function xA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    xA = 1.2D+0*cos((5.0D-1*pi*(eta-1.0D+0))/(ny-2.0D+0)-7.5D-1*pi)&
      &*((1.0D+0*(xi-1.0D+0))/(nx-2.0D+0)+5.0D-1)
  end function


  real(8) pure function yA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    yA = -5.0D-1*sin((5.0D-1*pi*(eta-1.0D+0))/(ny-2.0D+0)-7.5D-1*pi&
      &)*((1.0D+0*(xi-1.0D+0))/(nx-2.0D+0)+5.0D-1)
  end function


  real(8) pure function xksi(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    xksi = (1.2D+0*cos((5.0D-1*pi*(eta-1.0D+0))/(ny-2.0D+0)-7.5D-1*&
      &pi))/(nx-2.0D+0)
  end function


  real(8) pure function xeta(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    xeta = -(6.0D-1*pi*sin((5.0D-1*pi*(eta-1.0D+0))/(ny-2.0D+0)-7.5&
      &D-1*pi)*((1.0D+0*(xi-1.0D+0))/(nx-2.0D+0)+5.0D-1))/(ny-2.0D+0)
  end function


  real(8) pure function yksi(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    yksi = -(5.0D-1*sin((5.0D-1*pi*(eta-1.0D+0))/(ny-2.0D+0)-7.5D-1*&
      &pi))/(nx-2.0D+0)
  end function


  real(8) pure function yeta(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    yeta = -(2.5D-1*pi*cos((5.0D-1*pi*(eta-1.0D+0))/(ny-2.0D+0)-7.5&
      &D-1*pi)*((1.0D+0*(xi-1.0D+0))/(nx-2.0D+0)+5.0D-1))/(ny-2.0D+0)
  end function


  real(8) pure function JA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    JA = 1.d0/(xksi(nx, ny, xi, eta)*yeta(nx, ny, xi, eta)-xeta(nx, ny, xi, eta)*yksi(nx, ny, xi, eta))
  end function


  real(8) pure function alphaA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    alphaA = xeta(nx, ny, xi, eta)**2.d0+yeta(nx, ny, xi, eta)**2.d0
  end function


  real(8) pure function betaA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    betaA = xksi(nx, ny, xi, eta)*xeta(nx, ny, xi, eta)+yksi(nx, ny, xi, eta)*yeta(nx, ny, xi, eta)
  end function


  real(8) pure function gammaA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    gammaA = xksi(nx, ny, xi, eta)**2.d0+yksi(nx, ny, xi, eta)**2.d0
  end function

end module
