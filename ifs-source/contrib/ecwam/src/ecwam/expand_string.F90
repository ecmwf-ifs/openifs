! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE EXPAND_STRING(     &
     &     myproc,                  &     ! %p
     &     nproc,                   &     ! %n
     &     timestep,                &     ! %t
     &     max_timestep,            &
     &     s,                       &     ! %s
     &     n)

!    S. SAARINEM   ECMWF   MAY 1996   MESSAGE PASSING

!*    PURPOSE.
!     --------

!             EXPAND STRING BY ADDING AN EXTENSION WHICH IS FUNCTION OF
!       THE PROCESS NUMBER, AND/OR THE TOTAL NUMBER OF PROCESSES,
!       AND/OR THE TIME STEP NUMBER, AND/OR THE INDEX OF THE STRING

!**   INTERFACE.
!     ----------
!     CALL EXPAND_STRING(MYPROC,NPROC,TIMESTEP,MAX_TIMESTEP,S,N)

!     *MYPROC*       INTEGER  PROCESS NUMBER (PE)
!     *NPROC*        INTEGER  TOTAL NUMBER OF PE'S
!     *TIMESTEP*     INTEGER  TIME STEP NUMBER (OR ANY OTHER INDEX)
!     *MAX_TIMESTEP* INTEGER  MAXIMUM TIME STEP NUMBER (OR ANY OTHER
!                             INDEX)
!     *S*            CHARACTER OR CHARACTER ARRAY THAT WILL BE EXPANDED
!     *N*            INTEGER  DIMENSION OF S 


!     METHOD.
!     -------
!     THE STRING S MUST CONTAIN THE FOLLOWING SUB-STRING THAT WILL BE
!     SUBSTITUTED BY THE CORRESPONDING EXPANDED SUB STRING, THAT STARTS
!     WITH % AND p, n, t, or s

!     %p   -->  myproc in character mode, left justified
!     %n   -->  nproc  in character mode, left justified
!     %t   -->  timestep in character mode, left justified
!     %s   -->  the array index of S, left justified

!     e.g.   s='name_%n.%p'
!            call(myproc,nproc,0,0,s,1)

!     will produce, if nproc=11

!     on PE1,  s = name_11.1 
!     on PE2,  s = name_11.2 
!     on PE3,  s = name_11.3 

!     ...

!     on PE10,  s = name_11.10 
!     on PE11,  s = name_11.11


!     EXTERNALS.
!     ---------


!     REFERENCE.
!     ----------

!       NONE

!-----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      implicit none
      integer(KIND=JWIM), intent(in)  :: myproc, nproc
      integer(KIND=JWIM), intent(in)  :: timestep, max_timestep
      integer(KIND=JWIM), intent(in)  :: n
      character(len=*), intent(inout) :: s(n)
! === END OF INTERFACE BLOCK ===
      character(len=len(s))   :: t
      character(len=2*len(s)) :: tt
      integer(KIND=JWIM) :: i, j, jj, loc_p, len_t
      integer(KIND=JWIM) :: ndigs(4), num(4)
      character(len=6)   :: fmt(4)

      if (n < 1) return

!*    Setup output formats
      num(1) = myproc
      num(2) = max(nproc,myproc)
      num(3) = n
      num(4) = max(max_timestep,timestep)

!*    Count number of digits in each integer
      do j=1,4
        ndigs(j) = 1
        if (num(j) /= 0) then
          ndigs(j) = 1 + log10(dble(abs(num(j))))
          if (num(j) < 0) ndigs(j) = ndigs(j) + 1 ! Room for minus sign
        endif
        ndigs(j) = min(ndigs(j),9)
   ! M  ax 9 digits supported; i.e. '999 999 999'
        write(fmt(j),'("(i",i1,")")') ndigs(j)
      enddo


!*    Expand fields '%s', '%p', '%n' and '%t' with their values


!*    A special treatment with the sequence numbering
      if (n>1) then
        loc_p = index(s(1),'%s')
        if (loc_p > 0) then
          s(2:) = s(1)
        endif
      endif

      do i=1,n
        t = adjustl(s(i))//' '
        loc_p = index(t,'%')

        if (loc_p > 0) then
          len_t = len_trim(t)
          j = loc_p
          tt(:j-1) = t(:j-1)
          tt(j:) = ' '
          jj = j-1

          do while (j <= len_t)
            if (t(j:j) == '%') then
              j = j + 1
              if (j <= len_t) then
                select case ( t(j:j) ) 
                case ( 'p' )   ! myproc
                write(tt(jj+1:jj+ndigs(1)),fmt(1)) myproc
                jj = jj + ndigs(1)
                case ( 'n' )   ! nproc
                write(tt(jj+1:jj+ndigs(2)),fmt(2)) nproc
                jj = jj + ndigs(2)
                case ( 's' )   ! sequence number i=[1..n]
                write(tt(jj+1:jj+ndigs(3)),fmt(3)) i
                jj = jj + ndigs(3)
                case ( 't' )   ! timestep
                write(tt(jj+1:jj+ndigs(4)),fmt(4)) timestep
                jj = jj + ndigs(4)
                case default
                tt(jj+1:jj+2) = '%'//t(j:j)
                jj = jj + 2
                end select
              else
                tt(jj+1:jj+1) = '%'
                jj = jj + 1
              endif
            else
              tt(jj+1:jj+1) = t(j:j)
              jj = jj + 1
            endif
            j = j + 1
          enddo

          t = adjustl(tt)

!*   Get also rid of any blanks in the middle of the string

          len_t = len_trim(t)
          j = 1
          do while (j < len_t)
            if (t(j:j) == ' ') then
              t(j:) = t(j+1:)
              len_t = len_trim(t)
            else
              j = j + 1
            endif
          enddo

        endif

        s(i) = t
      enddo

      END SUBROUTINE expand_string
