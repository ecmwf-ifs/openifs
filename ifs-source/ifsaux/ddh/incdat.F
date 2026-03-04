! (C) Copyright 1989- Meteo-France.

      subroutine incdat(kdat,ktim,kinc,kdatn,ktimn)
C
C     given a date/time and an increment, computes the new date/time
C
c     kdat   input   date (yyyymmdd)
c     ktim   input   time (seconds)
c     kinc   input   increment time (seconds positive/negative)
c     kdatn  output  final date (yyyymmdd)
c     ktimn  output  final time (seconds) 
C                                                                               
C     PEDRO VITERBO         E.C.M.W.F.    07/96.                            
C
      if (kinc.eq.0) then
        kdatn=kdat
        ktimn=ktim
      else
        iyyyy=kdat/10000
        imd=mod(kdat,10000)
        imm=imd/100
        idd=mod(imd,100)
c
        ijul=julday(iyyyy,imm,idd)
c
        isi=ktim+kinc
        idi=isi/86400
        ijul1=ijul+idi
        if (isi.lt.0.and.mod(isi,86400).ne.0) ijul1=ijul1-1
c
        call invjul(ijul1,iyyyy1,imm1,idd1)
        kdatn=iyyyy1*10000+imm1*100+idd1
        if (isi.ge.0.or.mod(isi,86400).eq.0) then
          ktimn=mod(isi,86400)
        else
          ktimn=86400+mod(isi,86400)
        endif
      endif
c
      return
      end
