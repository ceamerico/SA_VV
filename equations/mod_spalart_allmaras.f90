module mod_spallart_allmaras
    implicit none
    
    ! All equations here are based on Allmaras 2012 paper. 
    
    !********************************************************************************
    !+ ddt + div - lap -gradgrad phiphi + gradgrad psi phi - prodest = 0            *    
    !                                                                               *
    !*******************************************************************************
    !+ d(rho * hat_nu) / dt                                                        *
    !+ nabla•(ro * hat_nu * vec(u)                                                 *
    !- nabla•[ro(nu + hat_nu) nabla hat_nu] / sigma                                *
    !- cb2 * ro * (nabla hat_nu) ** 2.d0 / sigma                                   *
    !+ (nu + hat_nu) nabla ro • nabla hat_nu / sigma                               *
    !-ro (P - D) = 0.d0                                                            *
    !*******************************************************************************   
    
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    !equation =  &                                                                                                                                                                                                                      -
    !    &                                      ddt(nx, ny, dt, roa, hat_nua, ro)                                                                 & ! Numerical approach for the temporal component                                     -
    !    & .plus.                               div_uds(nx, ny, re, rn, rp, Jp, roe, ron, Uce, Vcn, hat_nu)                                       & ! Numerical approach for the advective component.                                   -
    !    & .minus. ( 1.d0 / sigma )          *  laplacian_cds(nx, ny, re, rn, rp, Jp, Je, Jn, alphae, betan, betae, gamman, real_psie, real_psin) & ! Numerical approach for the laplacian component with real_psi at faces.            -
    !    & .plus.  (( ro + hat_nu )/ sigma ) *  grad_dot_grad_cds(nx, ny, Jp, alphap, betap, gammap, roe, ron, hat_nue, hat_nun)                  & ! Numerical approach for the first scalar product component.                        -
    !    & .minus. ( cb2 * ro / sigma )      *  grad_dot_grad_cds(nx, ny, Jp, alphap, betap, gammap, hat_nue, hat_nun, hat_nue, hat_nun)          & ! Numerical approach for the second scalar product component.                       -
    !    & .minus.                              psi_prod_minus_dest(nx, ny, d, omg, ro, vlp, hat_nu)                                              & ! Numerical approach for the second attempt for the psi*(P-D) simplified component. -      
    !------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    !*********************************************************************************************************************************************************

    ! Constants of the Spalart-Allmaras model
    real(8), parameter :: cv1   =     7.1d0                                         ! Viscous constant 1   
    real(8), parameter :: cb1   =  0.1355d0                                         ! Basic constant 1
    real(8), parameter :: cb2   =   0.622d0                                         ! Basic constant 2
    real(8), parameter :: cw2   =     0.3d0                                         ! Wall constant 2
    real(8), parameter :: cw3   =     2.0d0                                         ! Wall constant 3   
    real(8), parameter :: ct3   =     1.2d0                                         ! Transport constant 3
    real(8), parameter :: ct4   =     0.5d0                                         ! Transport constant 4
    real(8), parameter :: kappa =    0.41d0                                         ! Von Kármán constant
    real(8), parameter :: r_lim =    10.0d0                                         ! constant 
    real(8), parameter :: sigma = 2.d0/3.d0                                         ! constant
    real(8), parameter :: cw1   = cb1 / ( kappa ** 2.d0 ) + (1.d0 + cb2) / sigma    ! Wall constant 1
    real(8), parameter :: Pr_t =      0.9d0                                         ! Turbulent Prandtl number
    
    contains
    
    !This SR calculates:
    !   • minimum distance between each centroid and the next wall. Once the geometry may not be always strait or rectangular - "d".
    !   • the initialization of "hat_nu" equals to 5 times the kinematic viscosity. Also initialize the work variable at each face.
    !   • sets the hat_nua as zero for the fist iteration (it is valid only if the SA model initializes in the very fist timestep).
    !   • sets fist values for turbulent kinematic viscosity, based on the hat_nu initialization.
    subroutine initialize_spalart_allmaras (x, y, xp, yp, nx, ny, ro, vlp, hat_nu, hat_nua, hat_nue, hat_nun, d, vt)
    implicit none
    
        real(8), dimension (nx*ny), intent(in )  :: x       ! Coord. x of the northest corner of volume P.
        real(8), dimension (nx*ny), intent(in )  :: y       ! Coord. y of the northest corner of volume P.       
        real(8), dimension (nx*ny), intent(in)   :: xp      ! Coord. x of the centroid of volume P.  
        real(8), dimension (nx*ny), intent(in)   :: yp      ! Coord. y of the centroid of volume P.
        integer, intent(in)                      :: nx      ! Number of volumes in csi direction (real+fictitious).
        integer, intent(in)                      :: ny      ! Number of volumes in eta direction (real+fictitious).
        real(8), dimension (nx*ny), intent(in )  :: ro      ! Specific mass (absolute density) at center of volumes (kg/m3)     
        real(8), dimension (nx*ny), intent(in )  :: vlp     ! Laminar viscosity at center of volume P [Pa.s] 
        real(8), dimension (nx*ny), intent(out)  :: hat_nu  ! Work varible
        real(8), dimension (nx*ny), intent(out)  :: hat_nua ! Work varible from previous iteration
        real(8), dimension (nx*ny), intent(out)  :: hat_nue ! Work varible at east boundary
        real(8), dimension (nx*ny), intent(out)  :: hat_nun ! Work varible ar north boundary        
        real(8), dimension (nx*ny), intent(out)  :: d       ! Minimum distance from each centroid to the next wall
        real(8), dimension (nx*ny), intent(out)  :: vt      ! Turbulent kinematic viscosity. 

        
        !internal functions
        integer                   :: i, j, k, np, npn, npe
        real(8)                   :: dist
        real(8), dimension(nx*ny) :: nu
        real(8), dimension(nx*ny) :: chi
        real(8), dimension(nx   ) :: xw         ! Auxilliary vector for x coord. of the northest
        real(8), dimension(nx   ) :: yw         ! Auxilliary vector for y coord.
        

        hat_nu  = 0.d0
        hat_nua = 0.d0
        hat_nue = 0.d0
        hat_nun = 0.d0
        d       = 0.d0
        vt      = 0.d0
        nu      = 0.d0
        chi     = 0.d0
        xw      = 0.d0
        yw      = 0.d0
                
        !Kinematic viscosity.
        nu = vlp / ro
        
        !Work variable, Eq. (8) Allmaras 2012
        hat_nu = 4.d0 * nu
                
        !Aux. Function, Eq. (1) Allmaras 2012 
        chi = hat_nu / nu
        
        !Turbulent viscosity, Eq. (1) Allmaras 2012 
        vt = hat_nu * ( chi ** 3.d0 ) / ( chi ** 3.d0 + cv1 ** 3.d0 )
        
        !hat_nu - east faces of all volumes.
        do i = 1, nx-1
          do j = 2, ny-1

            np  = nx * (j-1) + i
            npe = np + 1

            hat_nue(np) = ( hat_nu(npe) + hat_nu(np) ) / 2.d0

          end do
        end do

        !hat_nu - north faces of all volumes
        do i = 2, nx-1
          do j = 1, ny-1

            np  = nx * (j-1) + i
            npn = np + nx

            hat_nun(np) = ( hat_nu(npn) + hat_nu(np) ) / 2.d0

          end do
        end do
                
        !Minimum distante from the centroid to the wall calculation-> "d".
        j = ny-1
        do i = 2, nx-1
          
            np   = nx * (j-1) + i
         
            xw(i) = ( x(np) + x(np-1) ) * 0.5d0
            yw(i) = ( y(np) + y(np-1) ) * 0.5d0
         
        end do 
        
        do j = 2, ny-1
            do i = 2, nx-1
             
                np   = nx * (j-1) + i
             
                !Frequently a strait line is the minimum distance between the centroid and the wall. But, as curvilinear geometries may be involved,
                !for throat and inlet regions at nozzles, the real "dist" runs and checks this distance. If this value is smaller than vertical strait line
                !it becomes the minimum distance d for volume np.
                d(np) = dsqrt( ( xp(np) - xw(i) ) ** 2.0d0 + ( yp(np) - yw(i) ) ** 2.0d0  )
             
                do k = 2, nx-1
                 
                    dist = dsqrt( ( xp(np) - xw(k) ) ** 2.0d0 + ( yp(np) - yw(k) ) ** 2.0d0  )
                 
                    if (dist < d(np)) then

                        d(np) = dist
                     
                    end if
                
                end do    
                
            end do
        end do
           
    end subroutine
    
    
  subroutine get_omg( nx, ny, x, y, Jp, u, v, omg ) ! Vorticity calculation
      implicit none
      
      integer, intent(in)                    :: nx  ! Number of volumes in csi direction (real+fictitious)
      integer, intent(in)                    :: ny  ! Number of volumes in eta direction (real+fictitious)
      real(8), dimension (nx*ny), intent(in) :: x   ! Coord. x of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in) :: y   ! Coord. y of the northest corner of volume P
      real(8), dimension (nx*ny), intent(in) :: Jp  ! Jacobian at the center of volume P
      real(8), dimension (nx*ny), intent(in) :: u   ! Cartesian velocity of the last iteraction
      real(8), dimension (nx*ny), intent(in) :: v   ! Cartesian velocity of the last iteraction  
      real(8), dimension (nx*ny), intent(out) :: omg   !vorticity at the real volumes centers       
      
      ! Auxiliary variables
      integer :: i, j, np, nps, npw, npe, npsw, npn
      real(8) :: xkp, ykp, xep, yep

      xkp = 0.d0
      xep = 0.d0
      ykp = 0.d0
      yep = 0.d0
        
      do j = 2, ny-1
         do i = 2, nx-1

            np   = nx * (j-1) + i
            nps  = np  - nx
            npn  = np + nx
            npw  = np  - 1
            npe  = np + 1
            npsw = nps - 1


            xkp = 0.5d0 * ( x(np) + x(nps) - x(npw) - x(npsw) )
            ykp = 0.5d0 * ( y(np) + y(nps) - y(npw) - y(npsw) )
            xep = 0.5d0 * ( x(np) + x(npw) - x(nps) - x(npsw) )
            yep = 0.5d0 * ( y(np) + y(npw) - y(nps) - y(npsw) )
            
            omg(np) = 0.5d0 * Jp(np) * abs(    & 
               + ( v(npe) - v(npw) ) * yep &
               - ( v(npn) - v(nps) ) * ykp &
               + ( u(npe) - u(npw) ) * xep &
               - ( u(npn) - u(nps) ) * xkp )

         end do

      end do

  end subroutine get_omg
    
    !This function calculates:
    !   • real value for psi that enters the lapacian operator.
    subroutine get_real_psi(nx, ny, ro, hat_nu, vlp, real_psi) 
        integer,                      intent(in)    :: nx      !< Number of real partitions + 2 fictitious in xi
        integer,                      intent(in)    :: ny      !< Number of real partitions + 2 fictitious in eta
        real(8), dimension(nx*ny),    intent(in)    :: ro      !< Specific mass (absolute density) at center of volumes (kg/m3)
        real(8), dimension(nx*ny),    intent(in)    :: hat_nu  !< Work variable at the center of the volume
        real(8), dimension(nx*ny),    intent(in)    :: vlp     !<Laminar viscosity at center of volume P  
        real(8), dimension(nx*ny),    intent(out)    :: real_psi     !< Combination of hat_nu, kinematic viscosity and absolute density that enters the laplacian operator
        
        real_psi = ro * ( hat_nu + vlp / ro )
        
    endsubroutine
    
    !Calculates the kinematic eddy viscosity and refresh absolute viscosity.
    subroutine get_turbulent_viscosity(nx, ny, ro, hat_nu, vlp, vt) 
        integer,                      intent(in)    :: nx      !< Number of real partitions + 2 fictitious in xi
        integer,                      intent(in)    :: ny      !< Number of real partitions + 2 fictitious in eta
        real(8), dimension(nx*ny),    intent(in)    :: ro      !< Specific mass (absolute density) at center of volumes (kg/m3)
        real(8), dimension(nx*ny),    intent(in)    :: hat_nu  !< Work variable at the center of the volume
        real(8), dimension(nx*ny),    intent(inout) :: vlp     !< Laminar viscosity at center of volume P  
        real(8), dimension(nx*ny),    intent(out)   :: vt      !< Combination of hat_nu, kinematic viscosity and absolute density that enters the laplacian operator
        
        !Auxiliary vectors.
        real(8), dimension(nx*ny)                   :: chi
        real(8), dimension(nx*ny)                   :: nu
        
        vt = 0.d0
        
        nu = vlp / ro
        
        chi = hat_nu / nu
        
        vt = hat_nu * chi**3.d0 / ( chi**3.d0 + cv1 ** 3.d0 ) 
        
        vlp = vlp + vt * ro
                
    endsubroutine

    !Calculates the kinematic eddy thermal conductivity and refresh thermal conductivity.
    subroutine get_turbulent_thermal_conductivity(nx, ny, ro, hat_nu, vlp, vt, cp, kp, kt) 
        integer,                      intent(in)    :: nx      !< Number of real partitions + 2 fictitious in xi
        integer,                      intent(in)    :: ny      !< Number of real partitions + 2 fictitious in eta
        real(8), dimension(nx*ny),    intent(in)    :: ro      !< Specific mass (absolute density) at center of volumes (kg/m3)
        real(8), dimension(nx*ny),    intent(in)    :: hat_nu  !< Work variable at the center of the volume
        real(8), dimension(nx*ny),    intent(in)    :: vlp     !<Laminar viscosity at center of volume P  
        real(8), dimension(nx*ny),    intent(in)    :: vt      !< Combination of hat_nu, kinematic viscosity and absolute density that enters the laplacian operator        
        real(8), dimension(nx*ny),    intent(in)    :: cp      !< Specific heat at const pressure (J/(kg.K))
        real(8), dimension(nx*ny),    intent(inout) :: kp      !< Thermal conductivity at center of volume P (W/m.K)        
        real(8), dimension(nx*ny),    intent(out)   :: kt      !< Combination of hat_nu, kinematic viscosity and absolute density that enters the laplacian operator
        
        kt = ( vt * ro ) * cp / Pr_t
        
        kp = kp + kt
        
    endsubroutine
    

    
        !!This function calculates:
        !! - rho * (Production - Destruction + T) in Eq(9) Allmaras 2012. 
        !! No value or function was considered for the T component since its values are considered zero and the transition location,
        !! for actual code aplication, is irrelevant, so T = 0.d0. 
        function psi_prod_minus_dest(nx, ny, d, omg, ro, vlp, hat_nu) result(b)
        integer,                      intent(in)    :: nx     !< Number of real partitions + 2 fictitious in xi
        integer,                      intent(in)    :: ny     !< Number of real partitions + 2 fictitious in eta
        real(8), dimension(nx*ny),    intent(in)    :: d    !< d
        real(8), dimension(nx*ny),    intent(in)    :: omg    !< Omg
        real(8), dimension(nx*ny),    intent(in)    :: ro    !< Function ro at the center of the volume P
        real(8), dimension(nx*ny),    intent(in)    :: vlp    !< Function vlp at the center of the volume P    
        real(8), dimension(nx*ny),    intent(in)    :: hat_nu   !< Function hat_nu at the center of the volume P   
        real(8), dimension(nx*ny)                   :: b      !< Source term

        integer :: i, j, np
        real(8) :: Prod, Dest, chi, ft2, fw   
        real(8) :: nu
        real(8) :: tilde_s
        real(8) :: fv1
        real(8) :: fv2
        real(8) :: r
        real(8) :: g
    
        Prod   = 0.d0
        Dest   = 0.d0
        chi    = 0.d0
        ft2 = 0.d0
        fw  = 0.d0
        b   = 0.d0
        nu  = 0.d0
        tilde_s = 0.d0
        fv1 = 0.d0
        fv2 = 0.d0
        g = 0.d0
    
        ! Calculating the source term for real volumes
        do i = 2, nx-1
          do j = 2, ny-1
    
            np  = nx * (j-1) + i
        
            nu = vlp(np) / ro(np)
                
            chi = hat_nu(np) / nu
        
            ft2 = ct3 * exp( -ct4 * chi ** 2.d0)
        
            fv1 = chi ** 3 / (chi ** 3 + cv1 ** 3)
        
            fv2 = 1.d0 - chi / ( 1.0d0 + chi * fv1 )
        
            tilde_s = omg(np) + hat_nu(np) * fv2 / (kappa**2 * d(np)**2)
            
            !For mach2d code, uncomment lines 328-336 and delete line 337.
            
            ! if (tilde_s < 0.0d0) then
            !     write(*,*) "negative value for tilde_s in the volume", np, "value =", tilde_s
            !     write(*,*) "physical accuracy not guaranteed"
            !     write(*,*) "if this error persist, change the modified vorticity calculation to Eq (11) and Eq (12)."
            !     stop
            ! end if
    
            !r = min(hat_nu(np) / (tilde_s * kappa ** 2 * d(np) ** 2), r_lim)
            
            r = hat_nu(np) / (tilde_s * kappa ** 2 * d(np) ** 2)
        
            g = r + cw2 * (r ** 6 - r)
        
            fw = g * ( (1.0d0 + cw3**6) / (g**6 + cw3**6) ) ** (1.0d0/6.0d0)
            
            Prod = cb1 * ( 1.d0 - ft2 ) * tilde_s * hat_nu(np) 
            
            Dest = (cw1 * fw - cb1 * ft2 / kappa**2 ) * (hat_nu(np) / d(np)) ** 2 

            b(np) = - ro(np) * ( Prod - Dest )
    
          end do
        end do
    
        end function
        
        !for mach2d code, remember:
        !If tilde_s results negative values, calculate it by this method:
        
        !s_overbar = hat_nu(np) * fv2 / (kappa ** 2.d0 * 1.d0 ** 2.d0)
        !
        !if (s_overbar >= -0.7d0 * omg(np) ) then
        !
        !    tilde_s = omg(np) + s_overbar
        !
        !    else
        !
        !    tilde_s = omg(np) + omg(np) * (omg(np) * 0.7d0 ** 2 + 0.9d0 * s_overbar) / ( ( 0.9d0 - 2.0d0 * 0.7d0 ) * omg(np) - s_overbar )
        !
        !end if  
    
    end module