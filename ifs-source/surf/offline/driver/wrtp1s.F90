      SUBROUTINE WRTP1S
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *WRTP1S*  - Writing prognostic variables of the one-column surface model

!     Purpose.
!     --------
!     Write out prognostic variables

!**   Interface.
!     ----------
!        *CALL* *WRTP1S

!        Explicit arguments :
!        --------------------


!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation 
!        of the one-column surface model

!     Author.
!     -------
!        Jean-Francois Mahfouf and Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-03-22

!     ------------------------------------------------------------------


      USE YOMLUN1S , ONLY : NULOUT   ,NPOSGG   ,NPOSGGL
      USE YOMDYN1S , ONLY : NSTEP
      USE YOMLOG1S , ONLY : LDBGS1
      USE YOMGP1S0 , ONLY : TSLNU0   ,QLINU0   ,TILNU0   ,FSNNU0   ,&
     &            TSNNU0   ,ASNNU0   ,RSNNU0   ,TRENU0   ,WRENU0   ,&
     &            TLICENU0,TLMNWNU0,TLWMLNU0,TLBOTNU0,TLSFNU0,& ! ENDUTRA
     &            HLICENU0,HLMLNU0                              ! ENDUTRA 

      USE YOMDPHY, ONLY : NCSS
      USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
      IMPLICIT LOGICAL (L)
      INTEGER(KIND=JPIM) ::IYYMD,IHM
      REAL(KIND=JPRD)    ::ZJUL


!     ------------------------------------------------------------------

!*       1.   SCALING SOME VARIABLES.
!             -----------------------
!  [different unit that TM295]    
!     NOTE: UNITS OF FSNNU,WRENU IN THE SURF1D MODEL:            Kg/m2
!           UNITS OF SNOW AND INTERCEP RESER WATER IN THE
!                                     PHYSICS:                   Kg/m2
!           UNITS OF SNOW AND INTERCEP RESER WATER IN THE
!                                     INPUT AND OUTPUT FILES:    Kg/m2


!     NOTE: UNITS OF QLINU IN THE SURF1D MODEL:                  M3 M-3
!           UNITS OF SOIL WATER IN THE INPUT FILE:               M3 M-3
!           UNITS OF QLINU IN THE PHYSICS:                       Kg/m2/layer




!*       2.    WRITE TIME STEP.
!              ----------------

      IF (LDBGS1) THEN
        WRITE(NULOUT,*) 'TIME STEP'
        WRITE(NULOUT,'(I6)') NSTEP

        WRITE(NULOUT,*) 'PROGNOSTIC VARIABLES'

!     -----------------------------------------------------------------


!*       3.    WRITE SOIL VARIABLES.
!              ---------------------

        WRITE(NULOUT,*) ' SOIL TEMPERATURE - SOIL MOISTURE - ICE T'
        DO JSLEV=1,NCSS
          WRITE(NULOUT,'(3(2X,E12.6))') TSLNU0(1,JSLEV),QLINU0(1,JSLEV),&
     &                TILNU0(1,JSLEV)
        ENDDO


!*       4.    WRITE SKIN VARIABLES.
!              ---------------------

        WRITE(NULOUT,*) ' SKIN TEMP. - SKIN. RES. CONT. '
        WRITE(NULOUT,'(2(2X,E12.6))') TRENU0(1) , WRENU0(1)


!*       5.    WRITE SNOW VARIABLES
!              -----------------

        WRITE(NULOUT,*) 'SNOW_DEPTH   SNOW_T  SNOW_ALBEDO  SNOW_DENSITY'
        WRITE(NULOUT,'(E12.6)') FSNNU0(1,1),TSNNU0(1,1),ASNNU0(1),RSNNU0(1,1)
      ENDIF

!        7.1   WRITE IN SURFSOIL FILE.
!        -----------------------------

      IF (NSTEP.EQ.0) THEN
        WRITE(NPOSGG,'(a)') '#     juld     date time      step'&
     &  //'        tskin         srwc'&
     &  //'     tsoil(1)     tsoil(2)     tsoil(3)     tsoil(4)'&
     &  //'     qsoil(1)     qsoil(2)     qsoil(3)     qsoil(4)'&
     &  //'     snowmass        snowt        snowa        snowd'&
     &  //'      tice(1)      tice(2)      tice(3)      tice(4)'
        WRITE(NPOSGG,'(a)') '#    jjj.j yyyymmdd hhmm          '&
     &  //'            K            m'&
     &  //'            K            K            K            K'&
     &  //'        m3m-3        m3m-3        m3m-3        m3m-3'&
     &  //'            m            K            -        kgm-3'&
     &  //'            K            K            K            K'

        WRITE(NPOSGGL,'(a)') '#     juld     date time      step'&
     &  //'        tlice        tlmnw        tlwml        tlbot         tlsf'&
     &  //'        Hlice         Hlml'
        WRITE(NPOSGGL,'(a)') '#    jjj.j yyyymmdd hhmm          '&
     &  //'            K            K            K            K            K'&
     &  //'            m            m'

      ENDIF
!   time stuff
      CALL DATTIM(ZJUL,IYYMD,IHM)
! to avoid trouble reading very small value (1e-XXX)
IF (WRENU0(1)<1e-30_JPRB) WRENU0=0._JPRB
      
     WRITE(NPOSGG,'(f10.3,1X,I8,1X,I4,2X,I8,18(1X,E12.6E3))')&
     &                 ZJUL,IYYMD,IHM,NSTEP&
     &                ,TRENU0(1),WRENU0(1)&
     &                ,TSLNU0(1,1),TSLNU0(1,2),TSLNU0(1,3),TSLNU0(1,4)&
     &                ,QLINU0(1,1),QLINU0(1,2),QLINU0(1,3),QLINU0(1,4)&
     &                ,FSNNU0(1,1),TSNNU0(1,1),ASNNU0(1),RSNNU0(1,1) &
     &                ,TILNU0(1,1),TILNU0(1,2),TILNU0(1,3),TILNU0(1,4)

     WRITE(NPOSGGL,'(f10.3,1X,I8,1X,I4,2X,I8,18(1X,E12.6E3))')&
     &                 ZJUL,IYYMD,IHM,NSTEP,&
     &                 TLICENU0,TLMNWNU0,TLWMLNU0,TLBOTNU0,TLSFNU0,& 
     &                 HLICENU0,HLMLNU0                              
!     --------------------------------------------------------

      RETURN
      END SUBROUTINE WRTP1S
