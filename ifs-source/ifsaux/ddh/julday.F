! (C) Copyright 1989- Meteo-France.

      FUNCTION JULDAY(IYYY,MM,ID)
      PARAMETER (IGREG=15+31*(10+12*1582))
      IF (IYYY .EQ. 0) write(*,*) ' WARNING There is no year zero'
      IF (IYYY .LT. 0) IYYY=IYYY+1
      IF (MM .GT. 2) THEN
         JY=IYYY
         JM=MM+1
      ELSE
         JY=IYYY-1
         JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY)+INT(30.6001*JM)+ID+1720995
      IF (ID+31*(MM+12*IYYY) .GE. IGREG) THEN
         JA=INT(0.01*JY)
         JULDAY=JULDAY+2-JA+INT(0.25*JA)
      ENDIF
      RETURN
      END
