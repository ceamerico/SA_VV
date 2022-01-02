module mod_spalart_allmaras_verification

  use mod_fabricated_solution1
  use mod_spallart_allmaras
  use mod_operator_base
  use mod_ddt
  use mod_div_uds
  use mod_laplacian_cds
  use mod_grad_dot_grad_cds
  use mod_class_solver_msi
  use mod_scarborough_conditioner
  use mod_norms
  use mod_metrics

  implicit none

  ! Vectors used in this module.

  ! Caution:
  ! It is not a good practice to create local variables in modules,
  ! because this makes the module less reusable in the same program.

  real(8), allocatable, dimension(:)   :: x
  real(8), allocatable, dimension(:)   :: y
  real(8), allocatable, dimension(:)   :: re
  real(8), allocatable, dimension(:)   :: rn
  real(8), allocatable, dimension(:)   :: rp
  real(8), allocatable, dimension(:)   :: Jp
  real(8), allocatable, dimension(:)   :: Jn
  real(8), allocatable, dimension(:)   :: Je
  real(8), allocatable, dimension(:)   :: alphae
  real(8), allocatable, dimension(:)   :: alphap
  real(8), allocatable, dimension(:)   :: betap
  real(8), allocatable, dimension(:)   :: betan
  real(8), allocatable, dimension(:)   :: betae
  real(8), allocatable, dimension(:)   :: gammap
  real(8), allocatable, dimension(:)   :: gamman
  real(8), allocatable, dimension(:)   :: Uce
  real(8), allocatable, dimension(:)   :: Vcn
  real(8), allocatable, dimension(:)   :: psie
  real(8), allocatable, dimension(:)   :: psin
  real(8), allocatable, dimension(:)   :: psi
  real(8), allocatable, dimension(:)   :: phie
  real(8), allocatable, dimension(:)   :: phin
  real(8), allocatable, dimension(:)   :: phi
  real(8), allocatable, dimension(:)   :: phiAS
  real(8), allocatable, dimension(:)   :: Source
  
  real(8), allocatable, dimension(:)   :: xp
  real(8), allocatable, dimension(:)   :: yp
  real(8), allocatable, dimension(:)   :: vt
  real(8), allocatable, dimension(:)   :: d
  real(8), allocatable, dimension(:)   :: hat_nu
  real(8), allocatable, dimension(:)   :: hat_nua
  real(8), allocatable, dimension(:)   :: hat_nun
  real(8), allocatable, dimension(:)   :: hat_nue
  real(8), allocatable, dimension(:)   :: vlp
  real(8), allocatable, dimension(:)   :: ro
  real(8), allocatable, dimension(:)   :: roa
  real(8), allocatable, dimension(:)   :: ron
  real(8), allocatable, dimension(:)   :: roe
  real(8), allocatable, dimension(:)   :: omg
  real(8), allocatable, dimension(:)   :: real_psi
  real(8), allocatable, dimension(:)   :: real_psie
  real(8), allocatable, dimension(:)   :: real_psin
  real(8), allocatable, dimension(:)   :: u
  real(8), allocatable, dimension(:)   :: v
    
  contains
    
    subroutine verify_SA_equation_code(nx, ny, dt)
    implicit none
    integer, intent(in) :: nx !< Number of real partitions + 2 fictitious in xi
    integer, intent(in) :: ny !< Number of real partitions + 2 fictitious in eta
    real(8), intent(in) :: dt !< Time stepping

    ! Inner variables
    integer :: it
    integer :: nitm    = 10
    real(8) :: tol     = 1.d-4
    integer :: nitmax  = 3000
    real(8) :: normtol = 1.D-14
    real(8) :: norm0   = 1.d0
    real(8) :: norm1   = 1.d0
    real(8) :: sc      = 1.99d0 ! Desired Scarborough constant
        
    real(8), dimension(nx*ny,10) :: equation

    type(class_solver_msi9d) :: solver

    ! Initializing the solver
    call solver%init(nx, ny, nitm, tol)

    ! Creating the vectors
    call create_vectors(nx, ny)

    ! Initializing the vectors
    call initialize_vectors(nx, ny)

    Jp     = JpCDS(nx, ny, x, y)
    alphap = alphapCDS(nx, ny, x, y)
    betap  = betapCDS(nx, ny, x, y)
    gammap = gammapCDS(nx, ny, x, y)
 
    !The initialize subroutine calculates the first value for hat_nu at center and volume surfaces. Also calculates d and vt. 
    !For this verification solution calculation, two major simplifications were made:
    ! • The entire vector d are equal to 1.d0, since none analytical equation for d was obtained.
    ! • The function r do not evaluate the min between f(tilde_s, hat_nu, d, kappa) and r_lim.
    
    !call initialize_spalart_allmaras (x, y, xp, yp, nx, ny, psi, psi**2, phi, hat_nua, phie, phin, d, vt)
    d = 1.d0
    
    ! Extrapolating psi to fictitious based on CDS and Dirichlet boundary conditions
    call extrapolate_to_fictitious_cds(nx, ny, psiA, psi)

    ! Interpolating psi on faces using psi at nodes
    call interpolate_on_all_faces_cds(nx, ny, psi, psie, psin)
    
    ! Initializing the outter cycle
    do it = 1, nitmax
        
        ! Since u and v fields are already prescribed, the vorticity magnitute do not change over iterations.
        ! For physical relevant problems, Jp should not become negative. If analytical and numerical omg becomes significantly different, the code may nt converge.
        call get_omg( nx, ny, x, y, Jp, u, v, omg ) 

        ! To overwrite this problem, vorticity becomes the analytical vorticity.
        omg = S_source(nx,ny,omgA)
        
        Source = &
          & -                         S_source(nx,ny,DivA)                & ! This line calculates the source for Div. analytical term.
          & + ( 1.d0 / sigma )     *  S_source(nx,ny,Lap_zetaA)           & ! This line calculates the source for Lap. analytical term.
          & - ((psi+phi) / sigma ) *  S_source(nx,ny,GradDotGradA)        & ! This line calculates the source for grad (phi) • grad (psi) analytical term. 
          & + ( cb2 * psi/ sigma ) *  S_source(nx,ny,GradDotGrad_phiphiA) & ! This line calculates the source for grad (phi) • grad (phi) analytical term.
          & +                         S_source(nx,ny, psi_p_minus_dA)       ! This line calculates the source for psi *(P-D) analytical component.
    
        ! Interpolating phi on faces
        call interpolate_on_all_faces_cds(nx, ny, phi, phie, phin)
      
        call get_real_psi(nx, ny, psi, phi, psi**2, real_psi)                       ! These two lines are responsible for real_psi calculation, where real_psi is:
        call interpolate_on_all_faces_cds(nx, ny, real_psi, real_psie, real_psin)   ! real_psi = ro*(phi+vlp/ro). And also, responsible to evaluate real_psi at faces.
      
        equation =  &
        &                                  ddt(nx, ny, dt, psi, phi, psi)                                                                    & ! Numerical approach for the temporal component - Since the manufactured solution for psi do not change at any cicle, psi = psia = ro = roa.
        & .plus.                           div_uds(nx, ny, re, rn, rp, Jp, psie, psin, Uce, Vcn, phi)                                        & ! Numerical approach for the advective component.
        & .minus. ( 1.d0 / sigma )      *  laplacian_cds(nx, ny, re, rn, rp, Jp, Je, Jn, alphae, betan, betae, gamman, real_psie, real_psin) & ! Numerical approach for the laplacian component with real_psi at faces.
        & .plus.  (( psi+phi )/ sigma ) *  grad_dot_grad_cds(nx, ny, Jp, alphap, betap, gammap, psie, psin, phie, phin)                      & ! Numerical approach for the first scalar product component.
        & .minus. ( cb2 * psi / sigma ) *  grad_dot_grad_cds(nx, ny, Jp, alphap, betap, gammap, phie, phin, phie, phin)                      & ! Numerical approach for the second scalar product component.
        & .minus.                          psi_prod_minus_dest(nx, ny, d, omg, psi, psi**2, phi)                                             & ! Numerical approach for the second attempt for the psi*(P-D) simplified component.        
        & .minus. Source                                                                                                                       ! Source term addition.

      ! Calculating the boundary conditions
      call boundary_conditions(nx, ny, phi, equation(:,1:9), equation(:,10))

      ! Applying the Scarborough conditioner
      call scarborough_conditioner(nx, ny, sc, phi, equation(:,1:9), equation(:,10))

      ! Solving the linear system
      call solver%solve(equation(:,1:9), equation(:,10), phi)

      ! Calculating the norm of the error
      norm1 = normL2(nx, ny, phi-phiAS)

      write(*,*) norm1

      if ( abs(norm1-norm0) < normtol ) exit

      norm0 = norm1

    end do
    
    ! For mach2d code, remember:
    ! After SA equation convergence, calculate the kinematic eddy viscosity and turbulent thermal conductivity. 
    ! and update absolute viscosity and thermal conductivity with: call get_turbulent_viscosity(nx, ny, ro, hat_nu, vlp, vt)  and call get_turbulent_thermal_conductivity(nx, ny, ro, hat_nu, vlp, vt, cp, kp, kt)     
    
    end subroutine
    
  ! Instantiates the vectors
  subroutine create_vectors(nx, ny)
    implicit none
    integer, intent(in) :: nx !< Number of real partitions + 2 fictitious in xi
    integer, intent(in) :: ny !< Number of real partitions + 2 fictitious in eta

    allocate(x(nx*ny))
    allocate(y(nx*ny))
    allocate(re(nx*ny))
    allocate(rn(nx*ny))
    allocate(rp(nx*ny))
    allocate(Jp(nx*ny))
    allocate(Jn(nx*ny))
    allocate(Je(nx*ny))
    allocate(alphae(nx*ny))
    allocate(alphap(nx*ny))
    allocate(betap(nx*ny))
    allocate(betan(nx*ny))
    allocate(betae(nx*ny))
    allocate(gammap(nx*ny))
    allocate(gamman(nx*ny))
    allocate(Uce(nx*ny))
    allocate(Vcn(nx*ny))
    allocate(psie(nx*ny))
    allocate(psin(nx*ny))
    allocate(psi(nx*ny))
    allocate(phie(nx*ny))
    allocate(phin(nx*ny))
    allocate(phi(nx*ny))
    allocate(phiAS(nx*ny))
    allocate(Source(nx*ny))
    
    allocate(xp(nx*ny))
    allocate(yp(nx*ny))
    allocate(vt(nx*ny))
    allocate(d(nx*ny)) 
    allocate(hat_nu(nx*ny)) 
    allocate(hat_nua(nx*ny))   
    allocate(hat_nue(nx*ny))     
    allocate(hat_nun(nx*ny))     
    allocate(vlp(nx*ny)) 
    allocate(ro(nx*ny)) 
    allocate(roa(nx*ny)) 
    allocate(ron(nx*ny)) 
    allocate(roe(nx*ny)) 
    allocate(omg(nx*ny)) 
    allocate(real_psi(nx*ny)) 
    allocate(real_psie(nx*ny)) 
    allocate(real_psin(nx*ny))   
    allocate(u(nx*ny))
    allocate(v(nx*ny))

    
    x    = 0.d0
    y    = 0.d0
    re   = 0.d0
    rn   = 0.d0
    rp   = 0.d0
    Jp   = 0.d0
    Jn   = 0.d0
    Je   = 0.d0
    alphae = 0.d0
    alphap = 0.d0
    betap = 0.d0
    betae = 0.d0
    betan = 0.d0
    gammap = 0.d0
    gamman = 0.d0
    Uce  = 0.d0
    Vcn  = 0.d0
    psie = 0.d0
    psin = 0.d0
    psi  = 0.d0
    phie = 0.d0
    phin = 0.d0
    phi  = 0.d0
    phiAS = 0.d0
    Source = 0.d0
    
    xp = 0.d0
    yp = 0.d0
    vt = 0.d0
    d = 0.d0
    hat_nu = 0.d0
    hat_nua = 0.d0
    hat_nue = 0.d0
    hat_nun = 0.d0
    vlp = 0.d0
    ro = 0.d0
    roa = 0.d0
    ron = 0.d0
    roe = 0.d0
    omg = 0.d0
    real_psi = 0.d0
    real_psie = 0.d0
    real_psin = 0.d0   
    u = 0.d0
    v = 0.d0

  end subroutine


  ! Initializes the vectors
  subroutine initialize_vectors(nx, ny)
    implicit none
    integer, intent(in) :: nx !< Number of real partitions + 2 fictitious in xi
    integer, intent(in) :: ny !< Number of real partitions + 2 fictitious in eta

    ! Inner variables
    integer :: i, j, np
    real(8) :: xi, eta


    ! Calculating the vertices of the volumes
    do i = 1, nx
      do j = 1, ny

        np  = nx * (j-1) + i

        ! Center of the east face
        xi  = dble(i)
        eta = dble(j)

        x(np) = xA(nx, ny, xi, eta)
        y(np) = yA(nx, ny, xi, eta)

      end do

    end do
    
    ! Calculating centroids of the volumes
    do i = 1, nx
      do j = 1, ny

        np  = nx * (j-1) + i

        ! Center of the volume
        xi  = dble(i)-0.5d0
        eta = dble(j)-0.5d0
        
        xp(np) = xA(nx, ny, xi, eta)
        yp(np) = yA(nx, ny, xi, eta)

      end do
    end do

    ! East faces of all volumes
    do i = 1, nx-1

      do j = 2, ny-1

        np  = nx * (j-1) + i

        ! Center of the east face
        xi  = dble(i)
        eta = dble(j)-0.5d0

        re(np)    = yA(nx, ny, xi, eta)
        psie(np)  = psiA(nx, ny, xi,eta)
        Je(np)    = JA(nx, ny, xi,eta)
        alphae(np)= alphaA(nx, ny, xi, eta)
        betae(np) = betaA(nx, ny, xi, eta)
        Uce(np)   = UcA(nx, ny, xi, eta)
      end do

    end do

    ! North faces of all volumes
    do i = 2, nx-1

      do j = 1, ny-1

        np  = nx * (j-1) + i

        ! Center of the east face
        xi  = dble(i)-0.5d0
        eta = dble(j)

        rn(np)    = yA(nx, ny, xi, eta)
        psin(np)  = psiA(nx, ny, xi,eta)
        Jn(np)    = JA(nx, ny, xi,eta)
        betan(np) = betaA(nx, ny, xi, eta)
        gamman(np) = gammaA(nx, ny, xi, eta)
        Vcn(np)    = VcA(nx, ny, xi, eta)
      end do

    end do

    ! Centers of the real volumes
    do i = 2, nx-1

      do j = 2, ny-1

        np  = nx * (j-1) + i

        ! Center of the volume
        xi  = dble(i)-0.5d0
        eta = dble(j)-0.5d0

        rp(np)     = yA(nx, ny, xi,eta)
        Jp(np)     = JA(nx, ny, xi,eta)
        alphap(np) = alphaA(nx, ny, xi,eta)
        betap(np)  = betaA(nx, ny, xi,eta)
        gammap(np) = gammaA(nx, ny, xi,eta)
        psi(np)    = psiA(nx, ny, xi,eta)
        phiAS(np)  = phiA(nx, ny, xi,eta)
        u(np)      = uA(nx, ny, xi, eta)
        v(np)      = vA(nx, ny, xi, eta)
      end do

    end do

  end subroutine


  ! Creates the vector of source term for the problem under study
  pure function S_source(nx, ny, f) result(b)
    implicit none
    integer,                   intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny)                :: b      !< Source of the linear system

    interface function
      real(8) pure function f(nx, ny, xi, eta)
        implicit none
        integer, intent(in) :: nx  !< Number of real partitions + 2 fictitious in xi
        integer, intent(in) :: ny  !< Number of real partitions + 2 fictitious in eta
        real(8), intent(in) :: xi  !< Coordinate in the numerical domain
        real(8), intent(in) :: eta !< Coordinate in the numerical domain
      end function
    end interface

    ! Inner variables
    integer :: i, j, np
    real(8) :: xi, eta


    ! Centers of the real volumes
    do i = 2, nx-1

      do j = 2, ny-1

        np  = nx * (j-1) + i

        ! Center of the volume
        xi  = dble(i)-0.5d0
        eta = dble(j)-0.5d0

        b(np) = f(nx, ny, xi, eta)

      end do

    end do

  end function


  ! Sets the boundary conditions
  subroutine boundary_conditions(nx, ny, phi, a, b)
    implicit none
    integer,                      intent(in)    :: nx !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny),    intent(in)    :: phi!< Estimate for the unknown phi at the center of the volumes
    real(8), dimension(nx*ny, 9), intent(inout) :: a  !< Coefficients of the linear system
    real(8), dimension(nx*ny),    intent(inout) :: b  !< Source of the linear system

    ! Inner variables
    integer :: i, j, np, npw, npe, nps, npn
    real(8) :: xi, eta, phic

    ! West boundary

    i = 1

    do j = 2, ny-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)
      eta = dble(j)-0.5d0

      ! Dirichlet boundary condition for the CDS scheme
      a(np,:) = 0.d0
      a(np,5) = 1.d0
      a(np,6) = 1.d0
      b(np)   = 2.0d0 * phiA(nx, ny, xi, eta)

    end do


    ! East boundary

    i = nx

    do j = 2, ny-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)-1.d0
      eta = dble(j)-0.5d0

      ! Dirichlet boundary condition for the CDS scheme
      a(np,:) = 0.d0
      a(np,5) = 1.d0
      a(np,4) = 1.d0
      b(np)   = 2.0d0 * phiA(nx, ny, xi, eta)

    end do


    ! South boundary

    j = 1

    do i = 2, nx-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)-0.5d0
      eta = dble(j)

      ! Dirichlet boundary condition for the CDS scheme
      a(np,:) = 0.d0
      a(np,5) = 1.d0
      a(np,8) = 1.d0
      b(np)   = 2.0d0 * phiA(nx, ny, xi, eta)

    end do


    ! North boundary

    j = ny

    do i = 2, nx-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)-0.5d0
      eta = dble(j)-1.d0

      ! Dirichlet boundary condition for the CDS scheme
      a(np,:) = 0.d0
      a(np,5) = 1.d0
      a(np,2) = 1.d0
      b(np)   = 2.0d0 * phiA(nx, ny, xi, eta)

    end do


    ! SW corner
    i = 1
    j = 1
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    phic = phiA(nx, ny, 1.d0, 1.d0)
    a(np,:) = 0.d0
    a(np,5) = 1.d0
    a(np,9) = 1.d0
    b(np)   = 4.d0*phic-phi(npn)-phi(npe)


    ! SE corner
    i = nx
    j = 1
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    phic = phiA(nx, ny, nx-1.d0, 1.d0)
    a(np,:) = 0.d0
    a(np,5) = 1.d0
    a(np,7) = 1.d0
    b(np)   = 4.d0 * phic - phi(npn) - phi(npw)

    ! NW corner
    i = 1
    j = ny
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    phic = phiA(nx, ny, 1.d0, ny-1.d0)
    a(np,:) = 0.d0
    a(np,5) = 1.d0
    a(np,3) = 1.d0
    b(np)   = 4.d0*phic - phi(nps) - phi(npe)

    ! NE corner
    i = nx
    j = ny
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    phic = phiA(nx, ny, nx-1.d0, ny-1.d0)
    a(np,:) = 0.d0
    a(np,5) = 1.d0
    a(np,1) = 1.d0
    b(np)   = 4.d0 * phic - phi(nps) - phi(npw)

  end subroutine

  subroutine interpolate_on_all_faces_cds(nx, ny, phi, phie, phin)
    integer,                      intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny),    intent(in)    :: phi    !< Function phi at the center of volume P
    real(8), dimension(nx*ny),    intent(inout) :: phie   !< Function phi at the center of face east of volume P
    real(8), dimension(nx*ny),    intent(inout) :: phin   !< Function phi at the center of face north of volume P

    ! Inner variables
    integer :: i, j, np, npn, npe

    ! East faces of all volumes
    do i = 1, nx-1
      do j = 2, ny-1

        np  = nx * (j-1) + i
        npe = np + 1

        phie(np) = ( phi(npe) + phi(np) ) / 2.d0

      end do
    end do

    ! North faces of all volumes
    do i = 2, nx-1
      do j = 1, ny-1

        np  = nx * (j-1) + i
        npn = np + nx

        phin(np) = ( phi(npn) + phi(np) ) / 2.d0

      end do
    end do

  end subroutine


  subroutine write_grid(nx, ny)
    implicit none
    integer,                      intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                      intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta

    ! Inner variables
    integer :: i, j, np
    real(8) :: xi, eta


    ! Centers of the real volumes
    do i = 1, nx-1

      do j = 1, ny-1

        np  = nx * (j-1) + i

        ! Center of the volume
        xi  = dble(i)
        eta = dble(j)

        write(10,*) xA(nx, ny, xi, eta), yA(nx, ny, xi, eta)

      end do

      write(10,*)

    end do

    ! Centers of the real volumes
    do j = 1, ny-1

      do i = 1, nx-1

        np  = nx * (j-1) + i

        ! Center of the volume
        xi  = dble(i)
        eta = dble(j)

        write(10,*) xA(nx, ny, xi, eta), yA(nx, ny, xi, eta)

      end do

      write(10,*)

    end do

  end subroutine



  ! Extrapolates a function to fictitious using Dirichlet bc and CDS
  subroutine extrapolate_to_fictitious_cds(nx, ny, f, fp)
    implicit none
    integer,                   intent(in)    :: nx !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in)    :: ny !< Number of real partitions + 2 fictitious in eta
    real(8), dimension(nx*ny), intent(inout) :: fp !< Estimate for the unknown phi at the center of the volumes

    interface function
      real(8) function f(nx, ny, xi, eta)
        implicit none
        integer, intent(in) :: nx  !< Number of real partitions + 2 fictitious in xi
        integer, intent(in) :: ny  !< Number of real partitions + 2 fictitious in eta
        real(8), intent(in) :: xi  !< Coordinate in the numerical domain
        real(8), intent(in) :: eta !< Coordinate in the numerical domain
      end function
    end interface

    ! Inner variables
    integer :: i, j, np, nps, npn, npw, npe, npsw, npse, npnw, npne
    real(8) :: xi, eta, fpc

    ! West boundary

    i = 1

    do j = 2, ny-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)
      eta = dble(j)-0.5d0

      ! Dirichlet boundary condition for the CDS scheme
      fp(np) = 2.0d0 * f(nx, ny, xi, eta)-fp(np+1)

    end do


    ! East boundary

    i = nx

    do j = 2, ny-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)-1.d0
      eta = dble(j)-0.5d0

      ! Dirichlet boundary condition for the CDS scheme
      fp(np) = 2.0d0 * f(nx, ny, xi, eta)-fp(np-1)

    end do


    ! South boundary

    j = 1

    do i = 2, nx-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)-0.5d0
      eta = dble(j)

      ! Dirichlet boundary condition for the CDS scheme
      fp(np) = 2.0d0 * f(nx, ny, xi, eta)-fp(np+nx)

    end do


    ! North boundary

    j = ny

    do i = 2, nx-1

      np  = nx * (j-1) + i

      ! Center of the face
      xi  = dble(i)-0.5d0
      eta = dble(j)-1.d0

      ! Dirichlet boundary condition for the CDS scheme
      fp(np) = 2.0d0 * f(nx, ny, xi, eta)-fp(np-nx)

    end do



    ! SW corner
    i = 1
    j = 1
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npne = npn + 1

    fpc = f(nx, ny, 1.d0, 1.d0)

    fp(np) = 4.d0*fpc-fp(npn)-fp(npe)-fp(npne)


    ! SE corner
    i = nx
    j = 1
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npnw = npn - 1

    fpc = f(nx, ny, nx-1.d0, 1.d0)

    fp(np) = 4.d0 * fpc - fp(npn) - fp(npw) - fp(npnw)

    ! NW corner
    i = 1
    j = ny
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npse = nps + 1

    fpc = f(nx, ny, 1.d0, ny-1.d0)

    fp(np) = 4.d0*fpc - fp(nps) - fp(npe) - fp(npse)

    ! NE corner
    i = nx
    j = ny
    np   = nx * (j-1) + i
    nps  = np - nx
    npn  = np + nx
    npw  = np - 1
    npe  = np + 1
    npsw = nps - 1

    fpc = f(nx, ny, nx-1.d0, ny-1.d0)

    fp(np) = 4.d0 * fpc - fp(nps) - fp(npw) - fp(npsw)

  end subroutine
  
  subroutine check_metrics(nx, ny)
    implicit none
    integer,                   intent(in) :: nx     !< Number of real partitions + 2 fictitious in xi
    integer,                   intent(in) :: ny     !< Number of real partitions + 2 fictitious in eta

    integer :: i, j, np
    real(8) :: xi, eta
    real(8) :: norma = 0.d0
    real(8) :: normb = 0.d0
    real(8) :: normc = 0.d0
    real(8) :: normd = 0.d0

    real(8) :: ap(nx*ny)
    real(8) :: bp(nx*ny)
    real(8) :: cp(nx*ny)
    real(8) :: dp(nx*ny)

    do i = 2, nx-1
      do j = 2, ny-1

        np  = nx * (j-1) + i

        xi  = dble(i)-0.5d0
        eta = dble(j)-0.5d0

        x(np) = xA(nx, ny, xi, eta)
        y(np) = yA(nx, ny, xi, eta)
        ap(np) = alphaA(nx, ny, xi, eta)
        bp(np) = betaA(nx, ny, xi, eta)
        cp(np) = gammaA(nx, ny, xi, eta)
        dp(np) = JA(nx, ny, xi, eta)

      end do
    end do

    norma = normL2(nx, ny, alphap*Jp-ap*dp)
    normb = normL2(nx, ny, betap*Jp-bp*dp)
    normc = normL2(nx, ny, gammap*Jp-cp*dp)

    print*, norma, normb, normc

  end subroutine


end module
