! (C) Copyright 1989- Meteo-France.

      SUBROUTINE SUARCH

!     INITIALISE PARAMETERS TO BE SELECTED

!     PEDRO VITERBO         E.C.M.W.F.    03/94.
#include "param.h"
#include "comddh.h"
#include "comcha.h"
#include "namus.h"
C      CHARACTER CCDUM*2,CCD*1

!     1. ESTABLISH WORK TO BE DONE

      IF (CARCSUI.EQ.'A1') THEN
      ELSE IF (CARCSUI.EQ.'A1') THEN
      ELSE IF (CARCSUI.EQ.'A2') THEN
      ELSE IF (CARCSUI.EQ.'A3') THEN
      ELSE IF (CARCSUI.EQ.'B1') THEN
      ELSE IF (CARCSUI.EQ.'B2') THEN
      ENDIF

      DO J=1,15
        IF (PARAM(J).EQ.CPARDUM) exit
      ENDDO
      IPAR=J-1

      DO J=1,100
        IF (PROCESS(J).EQ.CPRODUM) exit
      enddo
      IPRO=J-1

      NTPL=0

      DO JPAR=1,IPAR
        PRODO: DO JPRO=1,IPRO
          DO J=1,ND
            IF (PARAM(JPAR) == CDPAR(J).AND.PROCESS(JPRO) == CDPRO(J))
     &       THEN
              NTPL=NTPL+1
              NINDEX(NTPL)=J
              cycle prodo
            ENDIF
          enddo
          WRITE(*,"(' COMBINATION OF PARAM ',A3,' AND PROCESS ',A4,
     &          ' NOT AVAILABLE ')") PARAM(JPAR),PROCESS(JPRO)
        enddo prodo
      enddo

      IF (NTPL == 0) THEN
        WRITE(*,*) ' NO WORK TO BE DONE '
        STOP
      ENDIF

      DO J=1,NTPL
        JI=NINDEX(J)
        CDWRI(J)=CDPAR(JI)//CDPRO(JI)
        WRITE(*,"(' PARAM ',A3,' AND PROCESS ',A4,
     &    ' SELECTED FOR ARCHIVING')") CDPAR(JI),CDPRO(JI)
      enddo

      DO J=1,15
        IF (NLVPTS(J).EQ.-999) exit
      enddo
      NLVTOT=J-1

      do j=1,100
        if (nudom(j) == -999) exit
      enddo
      nudtot=j-1
      write(*,*) 'Domains selected: ',(nudom(j),j=1,nudtot)

      RETURN
      END SUBROUTINE SUARCH
