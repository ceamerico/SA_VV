module mtdma2d9
   !
   ! GOALS
   !
   ! Solves a linear system Ax = b, where A is a nine diagonal matrix.
   ! TDMA method is applied line-by-line. Lines are swept first in the x-direction
   ! and then in the y-direction. So, TDMA is applied twice.
   ! Uses OpenMP for parallel computation.
   !
   !
   ! FUNCTIONS AND SUBROUTINES
   !
   ! 1) set_tdma2d9_num_threads
   ! 2) get_tdma2d9_num_threads()
   ! 3) tdma2d9
   !
   integer, private :: tdma2d9_num_threads = 2

contains

   subroutine set_tdma2d9_num_threads(num)
      implicit none
      integer, intent(in) :: num

      tdma2d9_num_threads = num

   end subroutine set_tdma2d9_num_threads


   integer function get_tdma2d9_num_threads()
      implicit none

      get_tdma2d9_num_threads = tdma2d9_num_threads

   end function get_tdma2d9_num_threads


   subroutine tdma2d9(nx, ny, nitm, a, b, x)
      implicit none
      integer, intent(in) :: nx   ! Total number of volumes in the x direction
      integer, intent(in) :: ny   ! Total number of volumes in the y direction
      integer, intent(in) :: nitm ! Maximum number of iteractions
      real(8), dimension(nx*ny,9), intent(in)    :: a ! 9-diagonals matrix of coefficients
      real(8), dimension(nx*ny),   intent(in)    :: b ! Source
      real(8), dimension(nx*ny),   intent(inout) :: x ! Solution

      ! a(np,1) = SW coefficient
      ! a(np,2) = S  coefficient
      ! a(np,3) = SE coefficient
      ! a(np,4) = W  coefficient
      ! a(np,5) = P  coefficient
      ! a(np,6) = E  coefficient
      ! a(np,7) = NW coefficient
      ! a(np,8) = N  coefficient
      ! a(np,9) = NE coefficient

      ! Auxiliary variables
      integer :: i, j, it, npsw, nps, npse, npw, np, npe, npnw, npn, npne
      real(8) :: D
      real(8), dimension(nx*ny) :: P
      real(8), dimension(nx*ny) :: Q


      do it = 1, nitm

         !======================= BEGIN PARALLEL COMPUTATION ===================================

         !$OMP PARALLEL &
            !$OMP DEFAULT(NONE) &
            !$OMP PRIVATE( i, j, npsw, nps, npse, npw, np, npe, npnw, npn, npne, D) &
            !$OMP SHARED( nx, ny, a, b, x, P, Q) &
            !$OMP NUM_THREADS(tdma2d9_num_threads)

         !------------------------------------------------------------------
         !
         ! Solving linear system (x direction implicit, y direction explicit)
         !
         !------------------------------------------------------------------
         !$OMP DO
         do j = 1, ny
            !!$ write(*,*) "j: ", j, "TID: ", omp_get_thread_num()

            i = 1

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npse = nps + 1
            npne = npn + 1

            if ( j == 1 ) then ! South boundary

               D =  b(np) - a(np,8) * x(npn) - a(np,9) * x(npne)

            else if ( j == ny ) then ! North boundary

               D =  b(np) - a(np,2) * x(nps) - a(np,3) * x(npse)

            else

               D =  b(np) &
                  - a(np,2) * x(nps) - a(np,3) * x(npse) &
                  - a(np,8) * x(npn) - a(np,9) * x(npne)

            end if

            P(np) = - a(np,6) / a(np,5)

            Q(np) = D / a(np,5)

            do i = 2, nx-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npsw = nps - 1
               npse = nps + 1
               npnw = npn - 1
               npne = npn + 1

               if ( j == 1 ) then ! South boundary

                  D =  b(np) - a(np,7) * x(npnw) - a(np,8) * x(npn) - a(np,9) * x(npne)

               else if ( j == ny ) then ! North boundary

                  D =  b(np) - a(np,1) * x(npsw) - a(np,2) * x(nps) - a(np,3) * x(npse)

               else

                  D =  b(np) &
                     - a(np,1) * x(npsw) - a(np,2) * x(nps) - a(np,3) * x(npse) &
                     - a(np,7) * x(npnw) - a(np,8) * x(npn) - a(np,9) * x(npne)

               end if

               P(np) = - a(np,6) / ( a(np,5) + a(np,4) * P(np-1) )

               Q(np) = ( D - a(np,4) * Q(np-1) ) / ( a(np,5) + a(np,4) * P(np-1) )

            end do

            i = nx

            np   = nx * (j-1) + i
            nps  = np - nx
            npn  = np + nx
            npsw = nps - 1
            npnw = npn - 1

            if ( j == 1 ) then ! South boundary

               D =  b(np) - a(np,7) * x(npnw) - a(np,8) * x(npn)

            else if ( j == ny ) then ! North boundary


               D =  b(np) - a(np,1) * x(npsw) - a(np,2) * x(nps)

            else

               D =  b(np) &
                  - a(np,1) * x(npsw) - a(np,2) * x(nps) &
                  - a(np,7) * x(npnw) - a(np,8) * x(npn)

            end if

            P(np) = - a(np,6) / ( a(np,5) + a(np,4) * P(np-1) )

            Q(np) = ( D - a(np,4) * Q(np-1) ) / ( a(np,5) + a(np,4) * P(np-1) )

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
            np   = nx * (j-1) + i
            npn  = np + nx
            npw  = np - 1
            npe  = np + 1
            npnw = npn - 1
            npne = npn + 1

            if ( i == 1 ) then ! West boundary

               D = b(np) - a(np,6) * x(npe) - a(np,9) * x(npne)


            else if ( i == nx ) then ! East boundary

               D = b(np) - a(np,4) * x(npw) - a(np,7) * x(npnw)

            else

               D = b(np) &
                  - a(np,4) * x(npw) - a(np,7) * x(npnw) &
                  - a(np,6) * x(npe) - a(np,9) * x(npne)

            end if

            P(np) = - a(np,8) / a(np,5)

            Q(np) = D / a(np,5)

            do j = 2, ny-1

               np   = nx * (j-1) + i
               nps  = np - nx
               npn  = np + nx
               npw  = np - 1
               npe  = np + 1
               npsw = nps - 1
               npse = nps + 1
               npnw = npn - 1
               npne = npn + 1

               if ( i == 1 ) then ! West boundary

                  D = b(np) - a(np,3) * x(npse) - a(np,6) * x(npe) - a(np,9) * x(npne)

               else if ( i == nx ) then ! East boundary

                  D = b(np) - a(np,1) * x(npsw) - a(np,4) * x(npw) - a(np,7) * x(npnw)

               else

                  D = b(np) &
                     - a(np,1) * x(npsw) - a(np,4) * x(npw) - a(np,7) * x(npnw) &
                     - a(np,3) * x(npse) - a(np,6) * x(npe) - a(np,9) * x(npne)

               end if

               P(np) = - a(np,8) / ( a(np,5) + a(np,2) * P(np-nx) )

               Q(np) = ( D - a(np,2) * Q(np-nx) ) / ( a(np,5) + a(np,2) * P(np-nx) )

            end do

            j = ny

            np   = nx * (j-1) + i
            nps  = np - nx
            npw  = np - 1
            npe  = np + 1
            npsw = nps - 1
            npse = nps + 1

            if ( i == 1 ) then ! West boundary

               D = b(np) - a(np,3) * x(npse) - a(np,6) * x(npe)

            else if ( i == nx ) then ! East boundary

               D = b(np) - a(np,1) * x(npsw) - a(np,4) * x(npw)

            else

               D = b(np) &
                  - a(np,1) * x(npsw) - a(np,4) * x(npw) &
                  - a(np,3) * x(npse) - a(np,6) * x(npe)

            end if

            P(np) = - a(np,8) / ( a(np,5) + a(np,2) * P(np-nx) )

            Q(np) = ( D - a(np,2) * Q(np-nx) ) / ( a(np,5) + a(np,2) * P(np-nx) )

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

   end subroutine tdma2d9

end module mtdma2d9
