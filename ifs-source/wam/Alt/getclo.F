! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      INTEGER FUNCTION getclo(yaoptions, yaargument)

      INTEGER :: getcla

      CHARACTER(LEN=1) :: yolastarg
      CHARACTER(LEN=*) :: yaoptions, yaargument
      CHARACTER(LEN=120) :: arg

      INTEGER :: here, imorearg, ivarg
      DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
      DATA yolastarg / " " /
 
      arg=' '
      CALL getarg(here,arg)
 
      iol=len_trim(arg)
      IF (iol .EQ. 2 .AND. arg(1:1) .EQ. '-' .AND. ivarg .EQ. 0 ) THEN
        iol = len_trim(yaoptions)
        DO jl=1,iol
          getclo = 0
          IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
            getclo = ichar(arg(2:2))
            IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
              yolastarg=yaoptions(jl:jl)
              ivarg=1
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ELSEIF ( ivarg .EQ. 1 ) THEN
         WRITE(*,*) ' option -', yolastarg, ' requires arguments'
         getclo=-1
      ELSEIF (iol .EQ. 0) THEN
        getclo=0
      ELSE
         WRITE(*,*) 'illegal option: ',arg(1:iol)
         getclo=-1
      ENDIF
      here = here + 1
      RETURN

      ENTRY getcla(yaargument)

      getcla = 1
      CALL getarg(here,arg)
      IF ( arg (1:1) .NE. '-' ) THEN
        here = here + 1
        yaargument=arg
      ELSE
        IF (ivarg.EQ.1) THEN
          WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for',
     .    ' the option -',yolastarg
          getcla = -1
        ELSE
          getcla = 0
        ENDIF
      ENDIF
      ivarg=0

      RETURN
      END FUNCTION getclo
