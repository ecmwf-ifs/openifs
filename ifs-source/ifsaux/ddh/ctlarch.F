! (C) Copyright 1989- Meteo-France.

      SUBROUTINE CTLARCH

!     SCANS THROUGH THE PLOTTING OF THE LIMITED DOMAINS REQUIRED
                                                                               
!     PEDRO VITERBO         E.C.M.W.F.    12/94.                            

#include "param.h"
#include "comddh.h"
#include "comcha.h"
#include "comdat.h"

      REAL ZZZ(JPZLE),ZPP(JPZLE),ZWRITE(jplevp1,jpvr)
C ,ZZLV(JPTIM*15)
      REAL ZDP(JPZLE)
      CHARACTER*4 CDU1,CDU2
      CHARACTER*2 CCDUM*2
      CHARACTER*12 cdouf


!     1. INITIAL CONTROLS

      IF (NTPL.GT.JPVR) THEN
        WRITE(*,*) " INSUFFICIENT SIZE FOR THE NUMBER OF VARIABLES"
        WRITE(*,*) " TRY AGAIN WITH LESS VARIABLES"
        STOP
      ENDIF
      IF (AREATYP.EQ.'zon') THEN
         CCDUM='ZO'
      ELSE IF (AREATYP.EQ.'lar') THEN
         CCDUM='DL'
      ENDIF

      do jdom=1,nudtot
        iu=30+jdom
        ndom=nudom(jdom)
        write(cdouf,'(a,i3.3)') 'outarch_',ndom
        open(unit=iu,file=cdouf,status='REPLACE')
        call sucom2(.true.)
        if (lheader) then
          write(iu,"(a)")
     &'#  exp vp  i ty   west   north   east   south'//
     &               '   west2   north2   east2  south2'
          if (ntype == 4.or.ntype == 1) then
            write(iu,"('# ',a4,3(1x,i2),8(1x,f7.2))")
     &         cexp,nvpln,ndom,ntype,rwes,rnor
          else if (ntype == 3) then
            write(iu,"('# ',a4,3(1x,i2),8(1x,f7.2))")
     &         cexp,nvpln,ndom,ntype,rwes,rnor,reas,rsou
          else if (ntype == 2) then
            write(iu,"('# ',a4,3(1x,i2),8(1x,f7.2))")
     &         cexp,nvpln,ndom,ntype,rwes,rnor,reas,rsou,
     &         rwes2,rnor2,reas2,rsou2
          endif
          write(iu,
     &     "('#   idat itim     vdat  vtim  lev',1x,100(3x,a7,3x))")
     &        (cdwri(j),j=1,ntpl)
        endif
      enddo
      
!     2. RETRIEVES AND SCALES THE FIELDS

      DO JT=NSTEPB,NSTEPE-NSTEPI,NSTEPI
      
!     3. OPEN FILES

        if(jt == nstepb) then
          IF (CEXP(1:1) == '1'.OR.CEXP(1:1) == '2'.OR.CEXP(1:1) == '3'
     &    .OR.CEXP(1:1) == '4'.OR.CEXP(1:1) == '5'.OR.CEXP(1:1) == '6'
     &    .OR.CEXP(1:1) == '7'.OR.CEXP(1:1) == '8'.OR.CEXP(1:1) == '9')
     &   THEN
            READ(CEXP,"(I4)") IEXP
            WRITE(CDU2,"(I4.4)") IEXP
          ELSE
            CDU2=CEXP
          ENDIF
!     TEMPORARILY FOOLS OPONE
          NSTEP(2)=0
          NSTEP(1)=NSTEPB+NSTEPI
          WRITE(CDU1,"(I0)") NSTEP(1)
          CLFNAME(1)='DHF'//CCDUM//CDU2//'+'//CDU1
          WRITE(CDU1,"(I0)") NSTEP(2)
          CLFNAME(2)='DHF'//CCDUM//CDU2//'+'//CDU1
          CALL OPONE
          NSTEP(2)=NSTEPB
        else
          CALL SWAPTS(JT)
        endif

        domain_fill : do jdom=1,nudtot
          ndom=nudom(jdom)
          DO JD=1,NTPL
            INDEX=NINDEX(JD)
            CALL GETCOM1(NSTEP(2),INDEX,JPZLE,IVERT,IIP,
     &                 SECOND(1),SECOND(2),ZZZ,ZPP,ZDP)

!     4. Fill Write Fields
     
            DO JLEV=1,IVERT
              ZWRITE(jlev,JD)=ZZZ(JLEV)
            enddo
          enddo

!     5. FIND DATE, TIME

          if (jdom == 1) then
            call sucom2(.false.)
            if (CDwri(1)(4:7) == 'INIT') then
              isec=0
            else if ((cdwri(1)(1:3) == 'STL'.or.
     &                cdwri(1)(1:3) == 'STS'.or.
     &                cdwri(1)(1:3) == 'STT'.or.
     &                cdwri(1)(1:3) == 'STI'.or.
     &                cdwri(1)(1:3) == 'SSS'.or.
     &                cdwri(1)(1:3) == 'SII'.or.
     &                cdwri(1)(1:3) == 'SWL').and.
     & cdwri(1)(len_trim(cdwri(1)):len_trim(cdwri(1))) == '0') then
              isec=0
            else
              isec=int(second(1))
            endif
            call incdat(ndate,ntime*3600,isec,idatn,itimn)
            idn=mod(idatn,1000000)
            ztimewr=float(itimn)/(3600.)
          endif

!     6. Write out result file

          iu=30+jdom
          do jlev=1,ivert
            write(iu,'(2(i8,1x,f5.2,1x),i3,100(1x,e12.5))') 
     &        ndate,float(ntime),idatn,ztimewr,jlev,
     &        (zwrite(jlev,iv),iv=1,ntpl)
!      write(iu,'(1x,i10,1x,i5,1x,i3,100(1x,e12.5)') 
!     &         idn,itimn,jlev,(zwrite(jlev,iv),iv=1,ntpl)
          enddo
        enddo domain_fill
      enddo
      do jdom=1,nudtot
        iu=30+jdom
        close(iu)
      enddo
      RETURN
      END
