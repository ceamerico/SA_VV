module mod_fabricated_functions_operators1

  implicit none

  ! Coefficient to perturb mass conservation. mcc=0 => mass conservation.
  ! A small perturbation (e.g. 0.1) may be useful to check the robustness
  ! of the numerical scheme.
  real(8), parameter :: mcc = 0.1d0

contains


  real(8) pure function phi0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    phi0 = 3.0D-1*sin(x)*cos(y)+1.0D+0
  end function


  real(8) pure function psi0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    psi0 = 2.0D-1*sin(x*y)+1.0D+0
  end function


  real(8) pure function u0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    u0 = -((25.d0*mcc*cos(x)+25.d0)*sin(x)*y*sin(y)+((-50.d0*mcc*cos(x))-50.d0)*&
    &sin(x)*cos(y)-200.d0*mcc*cos(x)-200.d0)/((12.d0*sin(x)*cos(y)+40.d0)*sin(x*y)+&
    & 60.d0*sin(x)*cos(y)+200.d0)
  end function


  real(8) pure function v0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    v0 = -(1.25D-1*cos(x)*(mcc*sin(x)+1.0D+0)*y*cos(y))/((3.0D-1*sin(&
    &x)*cos(y)+1.0D+0)*(2.0D-1*sin(x*y)+1.0D+0))
  end function


  real(8) pure function Div0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    Div0 = ((-(y*(3.0D-1*sin(x)*cos(y)+1.0D+0)*((-25*mcc*sin(x)**2*y*&
      &sin(y))+cos(x)*(25*mcc*cos(x)+25)*y*sin(y)+50*mcc*sin(x)**2*cos(y&
      &)+cos(x)*((-50*mcc*cos(x))-50)*cos(y)+200*mcc*sin(x))*(2.0D-1*sin&
      &(x*y)+1.0D+0))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+2&
      &00))-(3.0D-1*cos(x)*y*cos(y)*((25*mcc*cos(x)+25)*sin(x)*y*sin(y)+&
      &((-50*mcc*cos(x))-50)*sin(x)*cos(y)-200*mcc*cos(x)-200)*(2.0D-1*s&
      &in(x*y)+1.0D+0))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)&
      &+200)-(2.0D-1*y**2*(3.0D-1*sin(x)*cos(y)+1.0D+0)*((25*mcc*cos(x)+&
      &25)*sin(x)*y*sin(y)+((-50*mcc*cos(x))-50)*sin(x)*cos(y)-200*mcc*c&
      &os(x)-200)*cos(x*y))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*co&
      &s(y)+200)+(y*(3.0D-1*sin(x)*cos(y)+1.0D+0)*((25*mcc*cos(x)+25)*si&
      &n(x)*y*sin(y)+((-50*mcc*cos(x))-50)*sin(x)*cos(y)-200*mcc*cos(x)-&
      &200)*(2.0D-1*sin(x*y)+1.0D+0)*(12*cos(x)*cos(y)*sin(x*y)+y*(12*si&
      &n(x)*cos(y)+40)*cos(x*y)+60*cos(x)*cos(y)))/((12*sin(x)*cos(y)+40&
      &)*sin(x*y)+60*sin(x)*cos(y)+200)**2+1.25D-1*cos(x)*(mcc*sin(x)+1.&
      &0D+0)*y**2*sin(y)-2.5D-1*cos(x)*(mcc*sin(x)+1.0D+0)*y*cos(y))/y
  end function


  real(8) pure function Lap0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    Lap0 = ((-3.0D-1*sin(x)*sin(y)*(2.0D-1*sin(x*y)+1.0D+0))-6.0D-1*s&
      &in(x)*y*cos(y)*(2.0D-1*sin(x*y)+1.0D+0)-6.0D-2*x*sin(x)*y*sin(y)*&
      &cos(x*y)+6.0D-2*cos(x)*y**2*cos(y)*cos(x*y))/y
  end function


  real(8) pure function GradDotGrad0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    GradDotGrad0 = 6.0D-2*cos(x)*y*cos(y)*cos(x*y)-6.0D-2*x*sin(x)*si&
      &n(y)*cos(x*y)
  end function
  
  real(8) pure function GradDotGrad_phiphi0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y  
    GradDotGrad_phiphi0 = 8.999999999999999d-2*sin(x)**2*sin(y)**2+8.&
        &999999999999999d-2*cos(x)**2*cos(y)**2
  end function
    
    real(8) pure function Lap_zeta0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    Lap_zeta0 = ((-3.d-1*sin(x)*sin(y)*(2.d-1*sin(x*y)+1.d+0)*(2.d-1*&
    &sin(x*y)+3.d-1*sin(x)*cos(y)+2.d+0))-6.d-1*sin(x)*y*cos(y)*(2.d-1&
    &*sin(x*y)+1.d+0)*(2.d-1*sin(x*y)+3.d-1*sin(x)*cos(y)+2.d+0)-6.d-2&
    &*x*sin(x)*y*sin(y)*cos(x*y)*(2.d-1*sin(x*y)+3.d-1*sin(x)*cos(y)+2&
    &.d+0)+6.d-2*cos(x)*y**2*cos(y)*cos(x*y)*(2.d-1*sin(x*y)+3.d-1*sin&
    &(x)*cos(y)+2.d+0)+3.d-1*cos(x)*y*cos(y)*(2.d-1*y*cos(x*y)+3.d-1*c&
    &os(x)*cos(y))*(2.d-1*sin(x*y)+1.d+0)-3.d-1*sin(x)*y*sin(y)*(2.d-1&
    &*x*cos(x*y)-3.d-1*sin(x)*sin(y))*(2.d-1*sin(x*y)+1.d+0))/y
    end function
    
    real(8) pure function omg0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    omg0 = abs(((25*mcc*cos(x)+25)*sin(x)*sin(y)-((-50*mcc*cos(x))-50&
    &)*sin(x)*sin(y)+(25*mcc*cos(x)+25)*sin(x)*y*cos(y))/((12*sin(x)*c&
    &os(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200)-(((25*mcc*cos(x)+25)*sin&
    &(x)*y*sin(y)+((-50*mcc*cos(x))-50)*sin(x)*cos(y)-200*mcc*cos(x)-2&
    &00)*((-12*sin(x)*sin(y)*sin(x*y))+x*(12*sin(x)*cos(y)+40)*cos(x*y&
    &)-60*sin(x)*sin(y)))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*co&
    &s(y)+200)**2+(1.25d-1*sin(x)*(mcc*sin(x)+1.d+0)*y*cos(y))/((3.d-1&
    &*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))-(1.25d-1*mcc*cos(x)&
    &**2*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0)&
    &)+(3.75d-2*cos(x)**2*(mcc*sin(x)+1.d+0)*y*cos(y)**2)/((3.d-1*sin(&
    &x)*cos(y)+1.d+0)**2*(2.d-1*sin(x*y)+1.d+0))+(2.5d-2*cos(x)*(mcc*s&
    &in(x)+1.d+0)*y**2*cos(y)*cos(x*y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(&
    &2.d-1*sin(x*y)+1.d+0)**2))
    end function
    
    real(8) pure function duy0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    duy0 =  (((25*mcc*cos(x)+25)*sin(x)*y*sin(y)+((-50*mcc*cos(x))-50&
    &)*sin(x)*cos(y)-200*mcc*cos(x)-200)*((-12*sin(x)*sin(y)*sin(x*y))&
    &+x*(12*sin(x)*cos(y)+40)*cos(x*y)-60*sin(x)*sin(y)))/((12*sin(x)*&
    &cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200)**2-((25*mcc*cos(x)+25)*&
    &sin(x)*sin(y)-((-50*mcc*cos(x))-50)*sin(x)*sin(y)+(25*mcc*cos(x)+&
    &25)*sin(x)*y*cos(y))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*co&
    &s(y)+200)

    end function
    
    real(8) pure function dvx0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y
    dvx0 =  (1.25d-1*sin(x)*(mcc*sin(x)+1.d+0)*y*cos(y))/((3.d-1*sin(&
    &x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))-(1.25d-1*mcc*cos(x)**2*y&
    &*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))+(3.&
    &75d-2*cos(x)**2*(mcc*sin(x)+1.d+0)*y*cos(y)**2)/((3.d-1*sin(x)*co&
    &s(y)+1.d+0)**2*(2.d-1*sin(x*y)+1.d+0))+(2.5d-2*cos(x)*(mcc*sin(x)&
    &+1.d+0)*y**2*cos(y)*cos(x*y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1&
    &*sin(x*y)+1.d+0)**2)

    end function

    real(8) pure function psi_p_minus_d0(x, y)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: y  
    psi_p_minus_d0 = (2.d-1*sin(x*y)+1.d+0)*(1.355d-1*(1.d+0-1.2d+0*exp(-(5.d-1*(3&
    &.d-1*sin(x)*cos(y)+1.d+0)**2)/(2.d-1*sin(x*y)+1.d+0)**2))*(3.d-1*&
    &sin(x)*cos(y)+1.d+0)*(abs(((25*mcc*cos(x)+25)*sin(x)*sin(y)-((-50&
    &*mcc*cos(x))-50)*sin(x)*sin(y)+(25*mcc*cos(x)+25)*sin(x)*y*cos(y)&
    &)/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200)-(((25*mcc&
    &*cos(x)+25)*sin(x)*y*sin(y)+((-50*mcc*cos(x))-50)*sin(x)*cos(y)-2&
    &00*mcc*cos(x)-200)*((-12*sin(x)*sin(y)*sin(x*y))+x*(12*sin(x)*cos&
    &(y)+40)*cos(x*y)-60*sin(x)*sin(y)))/((12*sin(x)*cos(y)+40)*sin(x*&
    &y)+60*sin(x)*cos(y)+200)**2+(1.25d-1*sin(x)*(mcc*sin(x)+1.d+0)*y*&
    &cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))-(1.2&
    &5d-1*mcc*cos(x)**2*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*&
    &sin(x*y)+1.d+0))+(3.75d-2*cos(x)**2*(mcc*sin(x)+1.d+0)*y*cos(y)**&
    &2)/((3.d-1*sin(x)*cos(y)+1.d+0)**2*(2.d-1*sin(x*y)+1.d+0))+(2.5d-&
    &2*cos(x)*(mcc*sin(x)+1.d+0)*y**2*cos(y)*cos(x*y))/((3.d-1*sin(x)*&
    &cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0)**2))+5.948839976204641d+0*(3&
    &.d-1*sin(x)*cos(y)+1.d+0)*(1.d+0-(3.d-1*sin(x)*cos(y)+1.d+0)/(((3&
    &.d-1*sin(x)*cos(y)+1.d+0)**4/(((3.d-1*sin(x)*cos(y)+1.d+0)**3/(2.&
    &d-1*sin(x*y)+1.d+0)**3+3.579109999999999d+2)*(2.d-1*sin(x*y)+1.d+&
    &0)**4)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))))-(3.d-1*sin(x)*cos(y)+1.d+&
    &0)**2*((6.494896984028202d+0*(0.3d0*((4.431938439877681d+4*(3.d-1*s&
    &in(x)*cos(y)+1.d+0)**6)/(abs(((25*mcc*cos(x)+25)*sin(x)*sin(y)-((&
    &-50*mcc*cos(x))-50)*sin(x)*sin(y)+(25*mcc*cos(x)+25)*sin(x)*y*cos&
    &(y))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200)-(((25*&
    &mcc*cos(x)+25)*sin(x)*y*sin(y)+((-50*mcc*cos(x))-50)*sin(x)*cos(y&
    &)-200*mcc*cos(x)-200)*((-12*sin(x)*sin(y)*sin(x*y))+x*(12*sin(x)*&
    &cos(y)+40)*cos(x*y)-60*sin(x)*sin(y)))/((12*sin(x)*cos(y)+40)*sin&
    &(x*y)+60*sin(x)*cos(y)+200)**2+(1.25d-1*sin(x)*(mcc*sin(x)+1.d+0)&
    &*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))-(&
    &1.25d-1*mcc*cos(x)**2*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d&
    &-1*sin(x*y)+1.d+0))+(3.75d-2*cos(x)**2*(mcc*sin(x)+1.d+0)*y*cos(y&
    &)**2)/((3.d-1*sin(x)*cos(y)+1.d+0)**2*(2.d-1*sin(x*y)+1.d+0))+(2.&
    &5d-2*cos(x)*(mcc*sin(x)+1.d+0)*y**2*cos(y)*cos(x*y))/((3.d-1*sin(&
    &x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0)**2))+5.948839976204641d+0&
    &*(3.d-1*sin(x)*cos(y)+1.d+0)*(1.d+0-(3.d-1*sin(x)*cos(y)+1.d+0)/(&
    &((3.d-1*sin(x)*cos(y)+1.d+0)**4/(((3.d-1*sin(x)*cos(y)+1.d+0)**3/&
    &(2.d-1*sin(x*y)+1.d+0)**3+3.579109999999999d+2)*(2.d-1*sin(x*y)+1&
    &.d+0)**4)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))))**6-(5.948839976204641d&
    &+0*(3.d-1*sin(x)*cos(y)+1.d+0))/(abs(((25*mcc*cos(x)+25)*sin(x)*s&
    &in(y)-((-50*mcc*cos(x))-50)*sin(x)*sin(y)+(25*mcc*cos(x)+25)*sin(&
    &x)*y*cos(y))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200&
    &)-(((25*mcc*cos(x)+25)*sin(x)*y*sin(y)+((-50*mcc*cos(x))-50)*sin(&
    &x)*cos(y)-200*mcc*cos(x)-200)*((-12*sin(x)*sin(y)*sin(x*y))+x*(12&
    &*sin(x)*cos(y)+40)*cos(x*y)-60*sin(x)*sin(y)))/((12*sin(x)*cos(y)&
    &+40)*sin(x*y)+60*sin(x)*cos(y)+200)**2+(1.25d-1*sin(x)*(mcc*sin(x&
    &)+1.d+0)*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1&
    &.d+0))-(1.25d-1*mcc*cos(x)**2*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d&
    &+0)*(2.d-1*sin(x*y)+1.d+0))+(3.75d-2*cos(x)**2*(mcc*sin(x)+1.d+0)&
    &*y*cos(y)**2)/((3.d-1*sin(x)*cos(y)+1.d+0)**2*(2.d-1*sin(x*y)+1.d&
    &+0))+(2.5d-2*cos(x)*(mcc*sin(x)+1.d+0)*y**2*cos(y)*cos(x*y))/((3.&
    &d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0)**2))+5.9488399762&
    &04641d+0*(3.d-1*sin(x)*cos(y)+1.d+0)*(1.d+0-(3.d-1*sin(x)*cos(y)+&
    &1.d+0)/(((3.d-1*sin(x)*cos(y)+1.d+0)**4/(((3.d-1*sin(x)*cos(y)+1.&
    &d+0)**3/(2.d-1*sin(x*y)+1.d+0)**3+3.579109999999999d+2)*(2.d-1*si&
    &n(x*y)+1.d+0)**4)+1.d+0)*(2.d-1*sin(x*y)+1.d+0)))))+(5.9488399762&
    &04641d+0*(3.d-1*sin(x)*cos(y)+1.d+0))/(abs(((25*mcc*cos(x)+25)*si&
    &n(x)*sin(y)-((-50*mcc*cos(x))-50)*sin(x)*sin(y)+(25*mcc*cos(x)+25&
    &)*sin(x)*y*cos(y))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(&
    &y)+200)-(((25*mcc*cos(x)+25)*sin(x)*y*sin(y)+((-50*mcc*cos(x))-50&
    &)*sin(x)*cos(y)-200*mcc*cos(x)-200)*((-12*sin(x)*sin(y)*sin(x*y))&
    &+x*(12*sin(x)*cos(y)+40)*cos(x*y)-60*sin(x)*sin(y)))/((12*sin(x)*&
    &cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200)**2+(1.25d-1*sin(x)*(mcc&
    &*sin(x)+1.d+0)*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(&
    &x*y)+1.d+0))-(1.25d-1*mcc*cos(x)**2*y*cos(y))/((3.d-1*sin(x)*cos(&
    &y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))+(3.75d-2*cos(x)**2*(mcc*sin(x)+&
    &1.d+0)*y*cos(y)**2)/((3.d-1*sin(x)*cos(y)+1.d+0)**2*(2.d-1*sin(x*&
    &y)+1.d+0))+(2.5d-2*cos(x)*(mcc*sin(x)+1.d+0)*y**2*cos(y)*cos(x*y)&
    &)/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0)**2))+5.9488&
    &39976204641d+0*(3.d-1*sin(x)*cos(y)+1.d+0)*(1.d+0-(3.d-1*sin(x)*c&
    &os(y)+1.d+0)/(((3.d-1*sin(x)*cos(y)+1.d+0)**4/(((3.d-1*sin(x)*cos&
    &(y)+1.d+0)**3/(2.d-1*sin(x*y)+1.d+0)**3+3.579109999999999d+2)*(2.&
    &d-1*sin(x*y)+1.d+0)**4)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))))))/((0.3d0*&
    &((4.431938439877681d+4*(3.d-1*sin(x)*cos(y)+1.d+0)**6)/(abs(((25*&
    &mcc*cos(x)+25)*sin(x)*sin(y)-((-50*mcc*cos(x))-50)*sin(x)*sin(y)+&
    &(25*mcc*cos(x)+25)*sin(x)*y*cos(y))/((12*sin(x)*cos(y)+40)*sin(x*&
    &y)+60*sin(x)*cos(y)+200)-(((25*mcc*cos(x)+25)*sin(x)*y*sin(y)+((-&
    &50*mcc*cos(x))-50)*sin(x)*cos(y)-200*mcc*cos(x)-200)*((-12*sin(x)&
    &*sin(y)*sin(x*y))+x*(12*sin(x)*cos(y)+40)*cos(x*y)-60*sin(x)*sin(&
    &y)))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200)**2+(1.&
    &25d-1*sin(x)*(mcc*sin(x)+1.d+0)*y*cos(y))/((3.d-1*sin(x)*cos(y)+1&
    &.d+0)*(2.d-1*sin(x*y)+1.d+0))-(1.25d-1*mcc*cos(x)**2*y*cos(y))/((&
    &3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))+(3.75d-2*cos(x&
    &)**2*(mcc*sin(x)+1.d+0)*y*cos(y)**2)/((3.d-1*sin(x)*cos(y)+1.d+0)&
    &**2*(2.d-1*sin(x*y)+1.d+0))+(2.5d-2*cos(x)*(mcc*sin(x)+1.d+0)*y**&
    &2*cos(y)*cos(x*y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1&
    &.d+0)**2))+5.948839976204641d+0*(3.d-1*sin(x)*cos(y)+1.d+0)*(1.d+&
    &0-(3.d-1*sin(x)*cos(y)+1.d+0)/(((3.d-1*sin(x)*cos(y)+1.d+0)**4/((&
    &(3.d-1*sin(x)*cos(y)+1.d+0)**3/(2.d-1*sin(x*y)+1.d+0)**3+3.579109&
    &999999999d+2)*(2.d-1*sin(x*y)+1.d+0)**4)+1.d+0)*(2.d-1*sin(x*y)+1&
    &.d+0))))**6-(5.948839976204641d+0*(3.d-1*sin(x)*cos(y)+1.d+0))/(a&
    &bs(((25*mcc*cos(x)+25)*sin(x)*sin(y)-((-50*mcc*cos(x))-50)*sin(x)&
    &*sin(y)+(25*mcc*cos(x)+25)*sin(x)*y*cos(y))/((12*sin(x)*cos(y)+40&
    &)*sin(x*y)+60*sin(x)*cos(y)+200)-(((25*mcc*cos(x)+25)*sin(x)*y*si&
    &n(y)+((-50*mcc*cos(x))-50)*sin(x)*cos(y)-200*mcc*cos(x)-200)*((-1&
    &2*sin(x)*sin(y)*sin(x*y))+x*(12*sin(x)*cos(y)+40)*cos(x*y)-60*sin&
    &(x)*sin(y)))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200&
    &)**2+(1.25d-1*sin(x)*(mcc*sin(x)+1.d+0)*y*cos(y))/((3.d-1*sin(x)*&
    &cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))-(1.25d-1*mcc*cos(x)**2*y*co&
    &s(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))+(3.75d&
    &-2*cos(x)**2*(mcc*sin(x)+1.d+0)*y*cos(y)**2)/((3.d-1*sin(x)*cos(y&
    &)+1.d+0)**2*(2.d-1*sin(x*y)+1.d+0))+(2.5d-2*cos(x)*(mcc*sin(x)+1.&
    &d+0)*y**2*cos(y)*cos(x*y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*si&
    &n(x*y)+1.d+0)**2))+5.948839976204641d+0*(3.d-1*sin(x)*cos(y)+1.d+&
    &0)*(1.d+0-(3.d-1*sin(x)*cos(y)+1.d+0)/(((3.d-1*sin(x)*cos(y)+1.d+&
    &0)**4/(((3.d-1*sin(x)*cos(y)+1.d+0)**3/(2.d-1*sin(x*y)+1.d+0)**3+&
    &3.579109999999999d+2)*(2.d-1*sin(x*y)+1.d+0)**4)+1.d+0)*(2.d-1*si&
    &n(x*y)+1.d+0)))))+(5.948839976204641d+0*(3.d-1*sin(x)*cos(y)+1.d+&
    &0))/(abs(((25*mcc*cos(x)+25)*sin(x)*sin(y)-((-50*mcc*cos(x))-50)*&
    &sin(x)*sin(y)+(25*mcc*cos(x)+25)*sin(x)*y*cos(y))/((12*sin(x)*cos&
    &(y)+40)*sin(x*y)+60*sin(x)*cos(y)+200)-(((25*mcc*cos(x)+25)*sin(x&
    &)*y*sin(y)+((-50*mcc*cos(x))-50)*sin(x)*cos(y)-200*mcc*cos(x)-200&
    &)*((-12*sin(x)*sin(y)*sin(x*y))+x*(12*sin(x)*cos(y)+40)*cos(x*y)-&
    &60*sin(x)*sin(y)))/((12*sin(x)*cos(y)+40)*sin(x*y)+60*sin(x)*cos(&
    &y)+200)**2+(1.25d-1*sin(x)*(mcc*sin(x)+1.d+0)*y*cos(y))/((3.d-1*s&
    &in(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))-(1.25d-1*mcc*cos(x)**&
    &2*y*cos(y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.d-1*sin(x*y)+1.d+0))+&
    &(3.75d-2*cos(x)**2*(mcc*sin(x)+1.d+0)*y*cos(y)**2)/((3.d-1*sin(x)&
    &*cos(y)+1.d+0)**2*(2.d-1*sin(x*y)+1.d+0))+(2.5d-2*cos(x)*(mcc*sin&
    &(x)+1.d+0)*y**2*cos(y)*cos(x*y))/((3.d-1*sin(x)*cos(y)+1.d+0)*(2.&
    &d-1*sin(x*y)+1.d+0)**2))+5.948839976204641d+0*(3.d-1*sin(x)*cos(y&
    &)+1.d+0)*(1.d+0-(3.d-1*sin(x)*cos(y)+1.d+0)/(((3.d-1*sin(x)*cos(y&
    &)+1.d+0)**4/(((3.d-1*sin(x)*cos(y)+1.d+0)**3/(2.d-1*sin(x*y)+1.d+&
    &0)**3+3.579109999999999d+2)*(2.d-1*sin(x*y)+1.d+0)**4)+1.d+0)*(2.&
    &d-1*sin(x*y)+1.d+0)))))**6+6.4d+1)**1.6666666666666666d-1-9.67281&
    &3801308746d-1*exp(-(5.d-1*(3.d-1*sin(x)*cos(y)+1.d+0)**2)/(2.d-1*&
    &sin(x*y)+1.d+0)**2)))
  end function
    
end module
