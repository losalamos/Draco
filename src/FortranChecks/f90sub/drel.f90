!----------------------------------*-F90-*----------------------------------
!
! file   FortranCheck/f90sub/drel.f90
! author Kelly Thompson
! date   Tuesday, Jun 12, 2012, 16:03 pm
! brief  Test F90 main linking against C++ library and calling a C++ function.
! note   Copyright (c) 2012-2014 Los Alamos National Security, LLC.
!        All rights reserved.
!---------------------------------------------------------------------------
! $Id$
!---------------------------------------------------------------------------

! Ref: http://gcc.gnu.org/onlinedocs/gfortran/Interoperable-Subroutines-and-Functions.html
!      http://fortranwiki.org/fortran/show/Generating+C+Interfaces

subroutine drelf90() bind(c)

  use iso_c_binding, only : c_size_t,C_NULL_CHAR 
  implicit none

  interface
     ! INCLUDE 'ds++/Release.hh'
     subroutine ec_release(release_string,maxlen) BIND(C, name="ec_release")
       use iso_c_binding, only: c_char,c_size_t
       implicit none
       character(kind=c_char,len=1), intent(out)  :: release_string
       integer(c_size_t), intent(in), value :: maxlen
     end subroutine ec_release
  end interface
  
  !----------------------------------------------------------------------
  ! Variable declarations

   integer, parameter :: maxlen = 80
   character(len=maxlen) :: release_string
   integer :: allpass
   
   !----------------------------------------------------------------------
   ! Initialization 

   allpass = 1
   release_string = repeat(' ',maxlen)
   release_string(maxlen:maxlen) = C_NULL_CHAR

   !----------------------------------------------------------------------
   ! Retrieve the version string from ds++
   call ec_release( release_string, len(release_string,kind=c_size_t) )
   
   print '(a)', trim(release_string)
   if( release_string(1:6) .eq. "Draco-" )then
      print '(a)', "Test: passed"
      print '(a)', "     Found 'Draco-' in release string."
   else
      print '(a)', "Test: failed"
      print '(a)', "     Did not find 'Draco-' in the release string."
      allpass = 0
   endif
   
   print '(a)', " "
   print '(a)', "*********************************************"
   if(allpass.eq.0)then
         print '(a)', "**** cppmain Test: FAILED."
      else
         print '(a)', "**** cppmain Test: PASSED."
      endif
   print '(a)', "*********************************************"

 end subroutine drelf90

!---------------------------------------------------------------------------
! end of drel.f90
!---------------------------------------------------------------------------
