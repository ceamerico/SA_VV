module mtdma2d5
   !
   ! GOALS
   !
   ! Solves a linear system Ax = b, where A is a five diagonal matrix.
   ! TDMA method is applied line-by-line. Lines are swept first in the x-direction
   ! and then in the y-direction. So, TDMA is applied twice.
   ! Uses OpenMP for parallel computation.
   !
   !
   ! FUNCTIONS AND SUBROUTINES
   !
   ! 1) set_tdma2d5_num_threads
   ! 2) get_tdma2d_num_threads()
   ! 3) tdma2d5
   !
   integer, private :: tdma2d5_num_threads = 2

contains

   subroutine set_tdma2d5_num_threads(num)
      implicit none
      integer, intent(in) :: num

      tdma2d5_num_threads = num

   end subroutine set_tdma2d5_num_threads


   integer function get_tdma2d5_num_threads()
      implicit none

      get_tdma2d5_num_threads = tdma2d5_num_threads

   end function get_tdma2d5_num_threads


   subroutine tdma2d5(nx, ny, nitm, a, b, x)
      implicit none
      integer, intent(in) :: nx   ! Total number of volumes in the x direction
      integer, intent(in) :: ny   ! Total number of volumes in the y direction
      integer, intent(in) :: nitm ! Maximum number of iteractions
      real(8), dimension(nx*ny,5), intent(in)    :: a ! 5-diagonals matrix of coefficients
      real(8), dimension(nx*ny),   intent(in)    :: b ! Source
      real(8), dimension(nx*ny),   intent(inout) :: x ! Solution

      ! a(np,1) = S coefficient
      ! a(np,2) = W coefficient
      ! a(np,3) = P coefficient
      ! a(np,4) = E coefficient
      ! a(np,5) = N coefficient

      ! Auxiliary variables
      integer :: i, j, np, npn, nps, it
      real(8) :: D
      real(8), dimension(nx*ny) :: P
      real(8), dimension(nx*ny) :: Q


      do it = 1, nitm


         !======================= BEGIN PARALLEL COMPUTATION ===================================

         !$OMP PARALLEL &
            !$OMP DEFAULT(NONE) &
            !$OMP PRIVATE( i, j, np, nps, npn, D) &
            !$OMP SHARED( nx, ny, a, b, x, P, Q) &
            !$OMP NUM_THREADS(tdma2d5_num_threads)

         !------------------------------------------------------------------
         !
         ! Solving linear system (x direction implict, y direction explicit)
         !
         !------------------------------------------------------------------

         !$OMP DO
         do j = 1, ny
            !!$ write(*,*) "j: ", j, "TID: ", omp_get_thread_num()

            i = 1

            np  = nx * (j-1) + i
            nps = np - nx
            npn = np + nx

            if ( j == 1 ) then ! South boundary

               D =  b(np) - a(np,5) * x(npn)

            else if ( j == ny ) then ! North boundary

               D =  b(np) - a(np,1) * x(nps)

            else

               D =  b(np) - a(np,5) * x(npn) - a(np,1) * x(nps)

            end if

            P(np) = - a(np,4) / a(np,3)

            Q(np) = D / a(np,3)


            do i = 2, nx

               np  = nx * (j-1) + i
               nps = np - nx
               npn = np + nx

               if ( j == 1 ) then ! South boundary

                  D =  b(np) - a(np,5) * x(npn)

               else if ( j == ny ) then ! North boundary

                  D =  b(np) - a(np,1) * x(nps)

               else

                  D =  b(np) - a(np,5) * x(npn) - a(np,1) * x(nps)

               end if

               P(np) = - a(np,4) / ( a(np,3) + a(np,2) * P(np-1) )

               Q(np) = ( D - a(np,2) * Q(np-1) ) / ( a(np,3) + a(np,2) * P(np-1) )

            end do

         end do

         !$OMP END DO

         !$OMP DO
         do j = 1, ny
            !!$ write(*,*) "j: ", j, "TID: ", omp_get_thread_num()
            i = nx

            np  = nx * (j-1) + i

            x(np) = Q(np)

            do i = nx-1, 1, -1

               np  = nx * (j-1) + i

               x(np) = P(np) * x(np+1) + Q(np)

            end do

         end do
         !$OMP END DO

         !------------------------------------------------------------------
         !
         ! Solving linear system (x direction explicit, y direction implicit)
         !
         !------------------------------------------------------------------

         !$OMP DO
         do i = 1, nx
            !!$ write(*,*) "i: ", i, "TID: ", omp_get_thread_num()

            j = 1
            np  = nx * (j-1) + i

            if ( i == 1 ) then ! West boundary

               D = b(np) - a(np,4) * x(np+1)

            else if ( i == nx ) then ! East boundary

               D = b(np) - a(np,2) * x(np-1)

            else

               D = b(np) - a(np,2) * x(np-1) - a(np,4) * x(np+1)

            end if

            P(np) = - a(np,5) / a(np,3)

            Q(np) = D / a(np,3)


            do j = 2, ny

               np  = nx * (j-1) + i

               if ( i == 1 ) then ! West boundary

                  D = b(np) - a(np,4) * x(np+1)

               else if ( i == nx ) then ! East boundary

                  D = b(np) - a(np,2) * x(np-1)

               else

                  D = b(np) - a(np,2) * x(np-1) - a(np,4) * x(np+1)

               end if

               P(np) = - a(np,5) / ( a(np,3) + a(np,1) * P(np-nx) )

               Q(np) = ( D - a(np,1) * Q(np-nx) ) / ( a(np,3) + a(np,1) * P(np-nx) )

            end do

         end do
         !$OMP END DO

         !$OMP DO
         do i = 1, nx
            !!$ write(*,*) "i: ", i, "TID: ", omp_get_thread_num()

            j = ny

            np  = nx * (j-1) + i

            x(np) = Q(np)

            do j = ny-1, 1, -1

               np  = nx * (j-1) + i

               x(np) = x(np+nx) * P(np) + Q(np)

            end do

         end do

         !$OMP END DO
         !$OMP END PARALLEL
         !======================== END PARALLEL COMPUTATION ====================================

      end do

   end subroutine tdma2d5

end module mtdma2d5
