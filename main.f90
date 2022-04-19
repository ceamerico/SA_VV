program main
  
  use mod_spalart_allmaras_verification

  implicit none

    integer, parameter :: grid = 8
    integer, parameter :: nt = 10 * 2 **(grid-1)
    integer, parameter :: nx = 10 * 2 **(grid-1)
    integer, parameter :: ny = 20 * 2 **(grid-1)
    real(8), parameter :: dt = 1.d0

    call verify_SA_equation_code(nx, ny, dt)
    
end program