      program makepixel
!      
!  Purpose:
!  =======
!
!     Generate a pixel file for AMADEUS 3D
!
      integer :: mode
!
!  Program History:
!  ===============
!
!     26-May-2008  G.Hobler  creation
!     18-Mar-2011  G.Hobler  adapted for IonShaper
!
!  Start of Executable Code:
!  ------------------------
!
      print *, 'Enter pixel mode:'
      print *, '   (1) circular - spiral'
      print *, '   (2) circular - raster'
      print *, '   (3) circular - serpentine (not implemented yet)'
      print *, '   (4) triangular (symmetric)'
      print *, ' '
!
      do
         read *, mode
         if( mode >= 1 .and. mode <= 4 ) exit
      end do
!
      open( 11, file='pixelfile.pix', form='formatted', status='unknown' )

      select case( mode )
      case( 1 )
         call circular_spiral( )
!             ===============
      case( 2 )
         call circular_raster( )
!             ===============
      case( 3 )
         call circular_serpentine( )
!             ===================
      case( 4 )
         call triangular_symmetric( )
!             ====================
      case default
         stop 'invalid mode'
      end select
!
      close( 11 )
!
      end
!======================================================================================      
!
      subroutine circular_spiral( )
!                ===============
!  Purpose:
!
!     Write pixel information to file unit 11 for a spiral scan (in circles starting
!     from the center)
!
!  Arguments:
!
!     none
!
!  Local Variables:
!
      real, parameter :: PI=3.141592654
!
      integer :: ir, nr, ip, np, nptot
      real :: dr, rmax, tdwell, phi, x, y
!
!  Start of Executable Code:
!  ------------------------
!
      print *, 'Enter radius of circle [nm]'
      read *, rmax
!
      print *, 'Enter pixel spacing [nm]'
      read *, dr
!
      print *, 'Enter dwell time [a.u.]'
      read *, tdwell
!
      nr = nint( rmax/dr )
!
!  Determine number of points
!
      nptot = 0
      do ir = 0, nr
         np = max( 1, nint( 2*PI*ir ) )
         nptot = nptot + np         
      end do
!
!      write( 11, * ) nptot
      write( *, * ) nptot, ' pixels generated.'
!
!  Determine points
!
      do ir = 0, nr
         np = max( 1, nint( 2*PI*ir ) )
         do ip = 1, np
            phi = 2*PI*real(ip-1)/real(np)
            x = ir*dr*cos(phi)
            y = ir*dr*sin(phi)
            write( 11, * ) x, y, tdwell
         end do
      end do
!
      end subroutine circular_spiral
!
!=======================================================================
!
      subroutine circular_raster( )
!                ===============
!  Purpose:
!
!     Write pixel information to file unit 11 for a raster scan inside a circle
!
!  Arguments:
!
!     none
!
!  Local Variables:
!
      integer :: nr, ix, iy, nptot
      real :: rmax, dr, tdwell, x, y
!
!  Start of Executable Code:
!  ------------------------
!
      print *, 'Enter radius of circle [nm]'
      read *, rmax
!
      print *, 'Enter pixel spacing [nm]'
      read *, dr
!
      print *, 'Enter dwell time [a.u.]'
      read *, tdwell
!
      nr = rmax / dr
!
!  Determine number of points
!
      nptot = 0
      do ix = -nr, nr
         x = ix * dr
         do iy = -nr, nr
            y = iy * dr
            if( x**2 + y**2 > rmax**2 ) cycle
            nptot = nptot + 1
         end do
      end do
!
      write( 11, * ) nptot
!
!  Determine points
!
      do ix = -nr, nr
         x = ix * dr
         do iy = -nr, nr
            y = iy * dr
            if( x**2 + y**2 > rmax**2 ) cycle
            write( 11, * ) x, y, tdwell
         end do
      end do            
!
      end subroutine circular_raster
!
!=======================================================================
!
      subroutine circular_serpentine( )
!                ===================
!  Purpose:
!
!
!
!  Arguments:
!

!
!  Local Variables:
!

!
!  Start of Executable Code:
!  ------------------------
!

!
      end subroutine circular_serpentine
!
!=======================================================================
!
      subroutine triangular_symmetric( )
!                ====================
!  Purpose:
!
!     Write pixel information to file unit 11 for a spiral scan (in circles starting
!     from the center)
!
!  Arguments:
!
!     none
!
!  Local Variables:
!
      integer :: nx, ny, ix, iy, nptot, nym
      real :: dp, td, tdwell, x, y
!
!  Start of Executable Code:
!  ------------------------
!
      print *, 'Enter pixel spacing [nm]'
      read *, dp
!
      print *, 'Enter number of pixels in x-,y-direction (preferably an odd number)'
      read *, nx, ny
!
      print *, 'Enter maximum dwell time [a.u.]'
      read *, tdwell
!
!  Determine number of points
!
      nptot = nx * (ny-2)
!
      write( 11, * ) nptot
!
!  Determine points
!
      nym = (ny+1)/2
      do iy = 2, ny-1
         y = 0.1 + (iy-1)*dp
         td = tdwell * ( 1.0 - real(abs(iy-nym)) / real(nym-1) )
         do ix = 1, nx
            x = 0.1 + (ix-1)*dp
            write( 11, '(f7.5,a,f7.5,a,f7.5)' ) x, y, td
         end do
      end do
!
      end subroutine triangular_symmetric
!
!=======================================================================
!
         
