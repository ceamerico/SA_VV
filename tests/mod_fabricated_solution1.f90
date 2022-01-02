module mod_fabricated_solution1

  ! mod_fabricated_functions_operators1 express the functions and operators
  ! in terms of (x,y), while mod_fabricated_grid1 gives x=x(xi,eta), y=y(xi,eta).
  ! This module combines the above mentioned modules to give functions and operators
  ! as functions of xi and eta.
  use mod_fabricated_grid1
  use mod_fabricated_functions_operators1

  implicit none

  private :: phi0, psi0, u0, v0, Div0, Lap0, GradDotGrad0, GradDotGrad_phiphi0, Lap_zeta0, omg0, dvx0, duy0, psi_p_minus_d0

contains

  real(8) pure function phiA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta

    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    phiA = phi0(x,y)

  end function

  real(8) pure function psiA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta

    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    psiA = psi0(x,y)

  end function


  real(8) pure function UcA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta

    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    UcA = u0(x,y) * yeta(nx, ny, xi, eta) - v0(x,y) * xeta(nx, ny, xi, eta)

  end function


  real(8) pure function VcA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta

    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    VcA = v0(x,y) * xksi(nx, ny, xi, eta) - u0(x,y) * yksi(nx, ny, xi, eta)

  end function


  real(8) pure function DivA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    DivA = Div0(x,y)

  end function

  real(8) pure function LapA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    LapA = Lap0(x,y)

  end function

  real(8) pure function GradDotGradA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    GradDotGradA = GradDotGrad0(x,y)

  end function
  
  real(8) pure function GradDotGrad_phiphiA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    GradDotGrad_phiphiA = GradDotGrad_phiphi0(x,y)

  end function
  
  real(8) pure function Lap_zetaA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    Lap_zetaA = Lap_zeta0(x,y)

  end function  
  
  real(8) pure function uA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta

    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    uA = u0(x,y)

  end function
  
  real(8) pure function vA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta

    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    vA = v0(x,y)

  end function
  
  
    real(8) pure function omgA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    omgA = omg0(x,y)

    end function  
    
    real(8) pure function dvxA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    dvxA = dvx0(x,y)

    end function  
    
    real(8) pure function duyA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    duyA = dvx0(x,y)

    end function     
    
    
    real(8) pure function psi_p_minus_dA(nx, ny, xi, eta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    real(8), intent(in) :: xi
    real(8), intent(in) :: eta
    real(8) :: x, y

    x = xA(nx, ny, xi, eta)
    y = yA(nx, ny, xi, eta)

    
    psi_p_minus_dA = psi_p_minus_d0(x,y)

  end function   
  
end module
