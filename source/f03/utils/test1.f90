      integer count,i
      real, dimension(50) :: r
      character(1), dimension(50) :: coins
      integer, dimension(1) :: seed
      real, dimension(2) :: a(2),b(2)

      call system_clock(count)
      seed = count

      do i=1,5
      call random_seed( put = seed )
      print*, seed
      call random_number( r)
      where (r < 0.5)
      coins = 'h'
      elsewhere
      coins = 't'
      end where
      print '(50a1)', coins
      call random_seed( get = seed) 
      end do

      a(1) = 2
      a(2) = 4
      b = 2*a
      print*,b

      end


