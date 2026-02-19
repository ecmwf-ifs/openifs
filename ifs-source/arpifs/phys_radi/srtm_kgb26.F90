SUBROUTINE SRTM_KGB26

!     Originally by J.Delamere, Atmospheric & Environmental Research.
!     Revision: 2.4
!     BAND 26:  22650-29000 cm-1 (low - nothing; high - nothing)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     G.Mozdzynski March 2011 read constants from files

!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOESRTA26 , ONLY : SFLUXREF, RAYL 

!     ------------------------------------------------------------------

IMPLICIT NONE

! KURUCZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SRTM_KGB26',0,ZHOOK_HANDLE)

SFLUXREF = (/ &
 !  &     129.462_JPRB, 15*_ZERO_ /)
 & 29.0079_JPRB,  28.4088_JPRB,     20.3099_JPRB,  13.0283_JPRB &
 & ,  11.8619_JPRB,  9.95840_JPRB,     6.68696_JPRB,  5.38987_JPRB &
 & ,  3.49829_JPRB, 0.407693_JPRB,    0.299027_JPRB, 0.236827_JPRB &
 & , 0.188502_JPRB, 0.163489_JPRB, 4.64335E-02_JPRB, 2.72662E-03_JPRB /)  

!     Rayleigh extinction coefficient at all v 
RAYL = (/ &
 & 1.21263E-06_JPRB,1.43428E-06_JPRB,1.67677E-06_JPRB,1.93255E-06_JPRB &
 & , 2.19177E-06_JPRB,2.44195E-06_JPRB,2.66926E-06_JPRB,2.85990E-06_JPRB &
 & , 3.00380E-06_JPRB,3.06996E-06_JPRB,3.08184E-06_JPRB,3.09172E-06_JPRB &
 & , 3.09938E-06_JPRB,3.10456E-06_JPRB,3.10727E-06_JPRB,3.10818E-06_JPRB /)  

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRTM_KGB26',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_KGB26
