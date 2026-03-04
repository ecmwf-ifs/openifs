SUBROUTINE DATES_DEMO
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
! --------------------------------------------------------------
!
! Conseils a l'utilisateur:
!
! 1. VOUS COMPILEZ LES ENTIERS EN 32 BITS:
! Utilisez alors les routines
! - ecartds: Ecart en secondes entre deux dates.
! - ecartdj: Ecart en jours entre deux dates.
! - dapluss: Date dans n secondes.
! - daplusj: Date dans n jours.
! - qqmmaa: Conversion d'un entier type AAAAQQMM vers une date en clair.
! - ijoursem: Jour de la semaine de la date d'entree.
! - quant: quantieme de l'annee d'une date donnee.
! Ces routines sont compatibles avec des entiers 32 bits.
! En effet elles appelent les routines citees ci-dessous, mais avec
! les parametres subsequents assurant que seuls des entiers
! representables en 32 bits y soient utilises.
!
! 2. VOUS COMPILEZ LES ENTIERS EN 64 BITS:
! Vous pouvez alors utiliser toutes les routines ci-dessus
! plus les suivantes, qui traitent des formats de dates
! en entree/sortie en JOURS, HEURES, MINUTES ou SECONDES:
! - ecartd: Ecart entre deux dates.
! - gregod: Conversion Date > Ecart par rapport a une date fixe.
! - gregoi: Conversion Ecart par rapport a une date fixe > Date.
! - daplus: Quelle sera la date dans n jours (ou heures, etc...)?
! - amqhms_vers_dj: Conversion date grégorienne (en 5 entiers et un réel) > date julienne.
! - dj_vers_amqhms: Conversion date julienne > date grégorienne (en 5 entiers et un réel).
! - amqhmsree_vers_dj: Conversion date grégorienne (en un seul réel) > date julienne.
! - dj_vers_amqhmsree: Conversion date julienne > date grégorienne (en un seul réel).
!
! --------------------------------------------------------------
!
! Définition des dates employées ci-dessous:
!
! Date julienne DJ:
!       Elle est composée d'un réel.
!       R1: Ce réel croît de 1 tous les jours,
!               et vaut 2451545.0 le 1er janvier 2000 à 12 UTC.
!
! Date grégorienne "en clair" AMQHMS:
!       Elle est composée de 5 entiers et d'un réel.
!       E1: Année (4 chiffres!)
!       E2: Mois
!       E3: Jour
!       E4: Heure
!       E5: Minute
!       R1: Seconde
! --------------------------------------------------------------


IMPLICIT NONE
END SUBROUTINE DATES_DEMO
SUBROUTINE DATE_PLUS_ECH(KAN,KMO,KQU,PSSSSS,PSTATI,CDTIT)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! Ecriture en clair d'une date de type BASE 2000.01.15 00:00 +72H VALID 2000.01.18 15:00.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! kan,kmo,kqu,psssss,pstati
! En sortie:
! cdtit
! --------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: KAN,KMO,KQU,IHE,IMI,IMIV,IHEV,IQUV,IMOV,IANV,ILZE
REAL(KIND=8) :: PSSSSS,PSTATI
REAL(KIND=8) :: ZS
REAL(KIND=8) :: ZSSSSS,ZDJ,ZSV
REAL(KIND=8) :: ZECH
CHARACTER*200 CLZUE,CLZE,CLECH
CHARACTER *(*) CDTIT
!
!-------------------------------------------------
! Date de validité.
!-------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DATE_PLUS_ECH',0,ZHOOK_HANDLE)
ZS=0.
ZSSSSS=PSSSSS/3600.
IHE=INT(ZSSSSS) ! heure de la base.
ZSSSSS=(ZSSSSS-REAL(IHE))*60.
IMI=INT(ZSSSSS) ! minute de la base.
ZSSSSS=ZSSSSS-REAL(IMI)
CALL AMQHMS_VERS_DJ(KAN,KMO,KQU,IHE,IMI,ZS,ZDJ)
ZDJ=ZDJ+PSTATI/86400. ! date julienne de validité.
CALL DJ_VERS_AMQHMS(ZDJ,IANV,IMOV,IQUV,IHEV,IMIV,ZSV) ! date grégorienne de validité.
IF(PSTATI < 3600.) THEN
!
!-------------------------------------------------
! Echéance en minutes.
!-------------------------------------------------
!
        ZECH=PSTATI/60. ; CLZUE='mn'
;ELSEIF(PSTATI < 259200.) THEN
!
!-------------------------------------------------
! Echéance en heures.
!-------------------------------------------------
!
        ZECH=PSTATI/3600. ; CLZUE='h'
;ELSE
!
!-------------------------------------------------
! Echéance en jours.
!-------------------------------------------------
!
        ZECH=PSTATI/86400. ; CLZUE='j'
;ENDIF
!
! Affichage de l'echeance avec deux chiffres apres la virgule.
!
WRITE(CLZE,FMT='(f9.2)') ZECH
!
! Si l'echeance est voisine d'un entier a mieux que 10**-2 pres,
! on l'affiche au format entier.
!
IF(CLZE(LEN_TRIM(CLZE)-2:LEN_TRIM(CLZE)) == '.00') THEN
        CLZE=CLZE(1:LEN_TRIM(CLZE)-3)
ENDIF
CLZE=ADJUSTL(CLZE)
ILZE=LEN_TRIM(CLZE)
CLECH=CLZE(1:ILZE)//CLZUE
!
!-------------------------------------------------
! Titre 3, de type
! BASE 2000.01.15 00:00 +72H VALID 2000.01.18 15:00.
!-------------------------------------------------
!
IF(IMI == 0 .AND. IMIV == 0) THEN
!
!-------------------------------------------------
! Les minutes de base et validité sont nulles.
! On ne les affiche pas.
!-------------------------------------------------
!
        WRITE(CDTIT,FMT='(a,i2,a,i2.2,a,i4.4,a,i2.2,3a,i2,a,i2.2,a,i4.4,a,i2.2,a)')&
        &'BASE ',KQU,'.',KMO,'.',KAN,' ',IHE,'h UTC + ',CLECH(1:LEN_TRIM(CLECH))&
        &,', VALID ',IQUV,'.',IMOV,'.',IANV,' ',IHEV,'h UTC'
ELSE
!
!-------------------------------------------------
! Les minutes de base ou validité sont non nulles.
! On les affiche.
!-------------------------------------------------
!
        WRITE(CDTIT,FMT='(a,i2,a,i2.2,a,i4.4,a,i2.2,a,i2.2,3a,i2,a,i2.2,a,i4.4,a,i2.2,a,i2.2,a)')&
        &'BASE ',KQU,'.',KMO,'.',KAN,' ',IHE,':',IMI,' UTC + ',CLECH(1:LEN_TRIM(CLECH))&
        &,' VALID ',IQUV,'.',IMOV,'.',IANV,' ',IHEV,':',IMIV,' UTC'
ENDIF
IF (LHOOK) CALL DR_HOOK('DATE_PLUS_ECH',1,ZHOOK_HANDLE)
END SUBROUTINE DATE_PLUS_ECH
SUBROUTINE DATC(KAAAA,KMM,KQQ,KHH,KMI,KSS,KJS,CDJS,CDDT)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *datc* Date courante machine.
! --------------------------------------------------------------
! Sujet:
! ------
! Arguments explicites:
! ---------------------
! Arguments implicites:
! ---------------------
! Methode:
! --------
! Externes:
! ---------
! Auteur:   95-05, J.M. Piriou.
! -------
! Modifications:
! --------------------------------------------------------------
! En entree:
! En sortie:
! kaaaa      annee.
! kmm      mois.
! kqq      quantieme.
! khh      heure.
! kmi      minute.
! kss      seconde.
! kjs      jour de la semaine (0: dimanche, 6 samedi).
! cdjs      jour de la semaine sur 3 caracteres (Dim, Lun, etc...).
! cddt      date totale (19950301-Mer-16:56:32).
! --------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IDATAT(8)
INTEGER(KIND=4) :: KJS
INTEGER(KIND=4) :: KSS
INTEGER(KIND=4) :: KMI
INTEGER(KIND=4) :: KHH
INTEGER(KIND=4) :: KQQ
INTEGER(KIND=4) :: KMM
INTEGER(KIND=4) :: KAAAA
INTEGER(KIND=4) :: IAAAAMMQQ
INTEGER(KIND=4) :: IJOURSEM
REAL(KIND=8) :: ZS
CHARACTER*200 CLGOL1,CLGOL2,CLGOL3
CHARACTER*3 CDJS
CHARACTER*(*) CDDT
CHARACTER*3 CLJOUR(0:6)
DATA CLJOUR/'Dim','Lun','Mar','Mer','Jeu','Ven','Sam'/
!
!-------------------------------------------------
! Date courante à la f90.
!-------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DATC',0,ZHOOK_HANDLE)
CLGOL1=' '
CLGOL2=' '
CLGOL3=' '
CALL DATE_AND_TIME(CLGOL1,CLGOL2,CLGOL3,IDATAT)
!
!-------------------------------------------------
! clgol1 est du type "AAAAMMQQ".
!-------------------------------------------------
!
READ(CLGOL1,FMT='(i4,2i2)') KAAAA,KMM,KQQ
!
!-------------------------------------------------
! clgol2 est du type "HHMMSS.SSS".
!-------------------------------------------------
!
READ(CLGOL2,FMT='(2i2)') KHH,KMI
READ(CLGOL2(5:),FMT=*) ZS
KSS=NINT(ZS)
READ(CLGOL1,FMT='(i8)') IAAAAMMQQ
!
!-------------------------------------------------
! Jour de la semaine.
!-------------------------------------------------
!
KJS=IJOURSEM(IAAAAMMQQ)
CDJS=CLJOUR(KJS)
!
!-------------------------------------------------
! Date totale.
!-------------------------------------------------
!
WRITE(CDDT,FMT='(i4.4,a,2(i2.2,a),2a,i2.2,a,i2.2,a,i2.2)') &
&KAAAA,'_',KMM,'_',KQQ,'_',CDJS,'_',KHH,':',KMI,':',KSS
IF (LHOOK) CALL DR_HOOK('DATC',1,ZHOOK_HANDLE)
END SUBROUTINE DATC
SUBROUTINE AMQHMS_VERS_DJ(KAAAA,KMM,KQQ,KHH,KMN,PS,PDJ)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *amqhms_vers_dj*
! --------------------------------------------------------------------------
! Auteur:
! -------
! 1999-08-17, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! kaaaa année (4 chiffres!)
! kmm   mois
! kqq   quantième du mois
! khh   heure
! kmn   minute
! ps    seconde
! En sortie:
! pdj date julienne associée à la date grégorienne UTC d'entrée
! --------------------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IDATE1
INTEGER(KIND=4) :: IDATE2
INTEGER(KIND=4) :: IECART
INTEGER(KIND=4) :: KAAAA
INTEGER(KIND=4) :: KHH
INTEGER(KIND=4) :: KMM
INTEGER(KIND=4) :: KMN
INTEGER(KIND=4) :: KQQ
REAL(KIND=8) :: PDJ
REAL(KIND=8) :: PS

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AMQHMS_VERS_DJ',0,ZHOOK_HANDLE)
IDATE1=20000101
IDATE2=KAAAA*10000+KMM*100+KQQ
!
!-------------------------------------------------
! Nombre de jours écoulés entre la date
! d'entrée à 0h UTC et le 1er janvier 2000 à 0h UTC.
!-------------------------------------------------
!
CALL ECARTDJ(IDATE1,IDATE2,IECART)
!
!-------------------------------------------------
! Date julienne.
!-------------------------------------------------
!
PDJ=2451545.0- 0.5 +REAL(IECART)+REAL(KHH)/24. &
& +REAL(KMN)/1440.+PS/86400.
IF (LHOOK) CALL DR_HOOK('AMQHMS_VERS_DJ',1,ZHOOK_HANDLE)
END SUBROUTINE AMQHMS_VERS_DJ
SUBROUTINE DAPLUS(KDAT1,KOPT,KDELT,KDAT2)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *DAPLUS* Quelle sera la date dans n jours (ou heures, etc...)?
! --------------------------------------------------------------------------
! Auteur:
! -------
! 94-10-31, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! kdat1
! kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! (cf. GREGOD).
! kdelt duree a ajouter a kdat1, unite: celle imposee par kopt.
! En sortie:
! kdat2 date finale.
!
! --------------------------------------------------------------------------
! Exemple: call DAPLUS(19940503,1,456,ires) fournira
! dans ires la date au format AAAAMMQQ situee 456 jours apres
! le 3 mai 1994.
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IGRE
INTEGER(KIND=4) :: KDAT1
INTEGER(KIND=4) :: KDAT2
INTEGER(KIND=4) :: KDELT
INTEGER(KIND=4) :: KOPT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DAPLUS',0,ZHOOK_HANDLE)
CALL GREGOD(KDAT1,KOPT,IGRE)
IGRE=IGRE+KDELT
CALL GREGOI(IGRE,KOPT,KDAT2)
IF (LHOOK) CALL DR_HOOK('DAPLUS',1,ZHOOK_HANDLE)
END SUBROUTINE DAPLUS
SUBROUTINE DAPLUSJ(K1,KEC,K2)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *daplusj* Date dans n jours.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! k1 date de depart au format AAAAMMQQ.
! kec nombre de jours ecoules.
! En sortie:
! k2 date d'arrivee au format AAAAMMQQ.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! PRECISION:
! Cette routine est utilisable avec des entiers 32 bits ou 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! -------------------------------------------------
! Date d'arrivee au jour pres.
! -------------------------------------------------
!


IMPLICIT NONE
INTEGER(KIND=4) :: K1
INTEGER(KIND=4) :: K2
INTEGER(KIND=4) :: KEC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DAPLUSJ',0,ZHOOK_HANDLE)
CALL DAPLUS(K1,1,KEC,K2)
IF (LHOOK) CALL DR_HOOK('DAPLUSJ',1,ZHOOK_HANDLE)
END SUBROUTINE DAPLUSJ
SUBROUTINE DAPLUSS(CD1,KEC,CD2)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *dapluss* Date dans n secondes.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! cd1 date de depart au format 'AAAAMMQQHHNNSS'.
! kec nombre de secondes ecoulees.
! En sortie:
! cd2 date d'arrivee au format 'AAAAMMQQHHNNSS'.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! Cette routine est utilisable avec des entiers 32 bits,
! si l'ecart entre les deux dates est inferieur a 2**31 secondes,
! soit 68 ans!...
!
! Au-dela de cette duree, les entiers doivent etre 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IAMQ1
INTEGER(KIND=4) :: IAMQ2
INTEGER(KIND=4) :: IDELTA
INTEGER(KIND=4) :: IECJOURS
INTEGER(KIND=4) :: IH1
INTEGER(KIND=4) :: IH2
INTEGER(KIND=4) :: IM1
INTEGER(KIND=4) :: IM2
INTEGER(KIND=4) :: IRESTE
INTEGER(KIND=4) :: IS1
INTEGER(KIND=4) :: IS2
INTEGER(KIND=4) :: ISEC
INTEGER(KIND=4) :: KEC
CHARACTER*(*) CD1,CD2
!
! -------------------------------------------------
! On lit les dates sur des entiers.
! -------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DAPLUSS',0,ZHOOK_HANDLE)
READ(CD1,FMT='(i8,3i2)') IAMQ1,IH1,IM1,IS1
!
! -------------------------------------------------
! Calculs d'ecarts et de leur partition
! en multiples de 86400 et sous-multiples.
! -------------------------------------------------
!
ISEC=IH1*3600+IM1*60+IS1 ! nombre de secondes ecoulees depuis cd10h.
IDELTA=KEC+ISEC ! nombre de secondes entre cd10h et cd2.
IRESTE=MODULO(IDELTA,86400) ! nombre de secondes entre cd20h et cd2.
IECJOURS=(IDELTA-IRESTE)/86400 ! nombre de jours entre cd10h et cd20h.
!
! -------------------------------------------------
! Date d'arrivee au jour pres.
! -------------------------------------------------
!
CALL DAPLUS(IAMQ1,1,IECJOURS,IAMQ2)
!
! -------------------------------------------------
! Date finale a la seconde pres.
! -------------------------------------------------
!
IH2=IRESTE/3600
IRESTE=IRESTE-3600*IH2
IM2=IRESTE/60
IRESTE=IRESTE-60*IM2
IS2=IRESTE
WRITE(CD2,FMT='(i8,3i2.2)') IAMQ2,IH2,IM2,IS2
IF (LHOOK) CALL DR_HOOK('DAPLUSS',1,ZHOOK_HANDLE)
END SUBROUTINE DAPLUSS
SUBROUTINE DJ_VERS_AMQHMS(PDJ,KAAAA,KMM,KQQ,KHH,KMN,PS)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *dj_vers_amqhms*
! --------------------------------------------------------------------------
! Auteur:
! -------
! 1999-08-17, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! pdj date julienne associée à la date grégorienne UTC d'entrée
! En sortie:
! kaaaa année (4 chiffres!)
! kmm   mois
! kqq   quantième du mois
! khh   heure
! kmn   minute
! ps    seconde
! --------------------------------------------------------------------------
!
!-------------------------------------------------
! Nombre de jours entre le 1er janvier 2000 à 0 UTC
! et la date julienne courante.
!-------------------------------------------------
!


IMPLICIT NONE
INTEGER(KIND=4) :: IDATE1
INTEGER(KIND=4) :: IDATE2
INTEGER(KIND=4) :: IECART
INTEGER(KIND=4) :: KAAAA
INTEGER(KIND=4) :: KHH
INTEGER(KIND=4) :: KMM
INTEGER(KIND=4) :: KMN
INTEGER(KIND=4) :: KNOUV
INTEGER(KIND=4) :: KQQ
REAL(KIND=8) :: PDJ
REAL(KIND=8) :: PS
REAL(KIND=8) :: ZECART
REAL(KIND=8) :: ZFRAC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DJ_VERS_AMQHMS',0,ZHOOK_HANDLE)
ZECART=PDJ-2451544.5
!
!-------------------------------------------------
! Nombre entier de jours.
!-------------------------------------------------
!
ZFRAC=MODULO(ZECART, 1._8 )
IECART=NINT(ZECART-ZFRAC)
!
!-------------------------------------------------
! Date grégorienne associée.
!-------------------------------------------------
!
IDATE1=20000101
CALL DAPLUSJ(IDATE1,IECART,IDATE2)
KQQ=MOD(IDATE2,100)
KNOUV=IDATE2/100
KMM=MOD(KNOUV,100)
KAAAA=KNOUV/100
!
!-------------------------------------------------
! Calcul de des heure, minute et seconde.
!-------------------------------------------------
!
ZFRAC=(ZECART-REAL(IECART))*24.
KHH=INT(ZFRAC)
ZFRAC=(ZFRAC-REAL(KHH))*60.
KMN=INT(ZFRAC)
PS=(ZFRAC-REAL(KMN))*60.
IF (LHOOK) CALL DR_HOOK('DJ_VERS_AMQHMS',1,ZHOOK_HANDLE)
END SUBROUTINE DJ_VERS_AMQHMS
SUBROUTINE DJ_VERS_AMQHMSREE(PDJ,PGRER)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** **
! --------------------------------------------------------------------------
! Auteur:
! -------
! 2002-11, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! pdj date julienne
! En sortie:
! pgrer date en clair au format AAAAMMQQ.HHMMSS
! --------------------------------------------------------------------------
!

IMPLICIT NONE
REAL(KIND=8), INTENT(IN) :: PDJ
REAL(KIND=8), INTENT(OUT) :: PGRER
REAL(KIND=8) :: ZS
INTEGER(KIND=4) :: IAAAA,IMM,IQQ,IHH,IMN
!
!-------------------------------------------------
! Conversion grégorien julien; cible 5 entiers et un réel.
!-------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DJ_VERS_AMQHMSREE',0,ZHOOK_HANDLE)
CALL DJ_VERS_AMQHMS(PDJ,IAAAA,IMM,IQQ,IHH,IMN,ZS)
!
!-------------------------------------------------
! On passe de ces 5 entiers et un réel à un seul réel.
!-------------------------------------------------
!
PGRER=REAL(IAAAA)*10000.+REAL(IMM)*100. &
& + REAL(IQQ)+REAL(IHH)/100. &
& + REAL(IMN)/10000.+ZS/1.E+06
IF (LHOOK) CALL DR_HOOK('DJ_VERS_AMQHMSREE',1,ZHOOK_HANDLE)
END SUBROUTINE DJ_VERS_AMQHMSREE
SUBROUTINE AMQHMSREE_VERS_DJ(PGRER,PDJ)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** **
! --------------------------------------------------------------------------
! Auteur:
! -------
! 2002-11, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! pgrer date en clair au format AAAAMMQQ.HHMMSS
! En sortie:
! pdj date julienne associée à la date grégorienne
! --------------------------------------------------------------------------
!

IMPLICIT NONE
REAL(KIND=8), INTENT(OUT) :: PDJ
REAL(KIND=8), INTENT(IN) :: PGRER
REAL(KIND=8) :: ZS,ZLOC
INTEGER(KIND=4) :: IAAAA,IMM,IQQ,IHH,IMN,ILOC
!
!-------------------------------------------------
! On passe de cette date grégorienne donnée
! comme un seul réel à 5 entiers et un réel.
!-------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AMQHMSREE_VERS_DJ',0,ZHOOK_HANDLE)
ILOC=INT(PGRER)
IQQ=MOD(ILOC,100)
ILOC=ILOC/100
IMM=MOD(ILOC,100)
IAAAA=ILOC/100

ILOC=NINT((PGRER-REAL(INT(PGRER)))*1.E+06)
ZS=REAL(MOD(ILOC,100))
ILOC=ILOC/100
IMN=MOD(ILOC,100)
IHH=ILOC/100
!
!-------------------------------------------------
! Conversion grégorien julien; cible 5 entiers et un réel.
!-------------------------------------------------
!
CALL AMQHMS_VERS_DJ(IAAAA,IMM,IQQ,IHH,IMN,ZS,PDJ)
IF (LHOOK) CALL DR_HOOK('AMQHMSREE_VERS_DJ',1,ZHOOK_HANDLE)
END SUBROUTINE AMQHMSREE_VERS_DJ
SUBROUTINE ECARTD(KDAT1,KDAT2,KOPT,KGRE)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *ECART* Ecart entre deux dates.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 97-01-09, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree: kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! kdat1 et kdat2 dates au format ci-dessus.
! En sortie:
! - si kopt=1 : kgre nombre de jours    entre kdat1 et kdat2
! - si kopt=2 : kgre nombre d'heures    entre kdat1 et kdat2
! - si kopt=3 : kgre nombre de minutes  entre kdat1 et kdat2
! - si kopt=4 : kgre nombre de secondes entre kdat1 et kdat2
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! L'ecart entre les deux dates doit etre inferieur a
! - 2**31 heures si kopt=2
! - 2**31 jours si kopt=1
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! L'ecart entre les deux dates doit etre inferieur a
! - 2**63 secondes si kopt=4
! - 2**63 minutes si kopt=3
! - 2**63 heures si kopt=2
! - 2**63 jours si kopt=1
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IGRE1
INTEGER(KIND=4) :: IGRE2
INTEGER(KIND=4) :: KDAT1
INTEGER(KIND=4) :: KDAT2
INTEGER(KIND=4) :: KGRE
INTEGER(KIND=4) :: KOPT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ECARTD',0,ZHOOK_HANDLE)
CALL GREGOD(KDAT1,KOPT,IGRE1)
CALL GREGOD(KDAT2,KOPT,IGRE2)
KGRE=IGRE2-IGRE1
IF (LHOOK) CALL DR_HOOK('ECARTD',1,ZHOOK_HANDLE)
END SUBROUTINE ECARTD
SUBROUTINE ECARTDJ(K1,K2,KEC)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *ecartdj* Ecart en jours entre deux dates.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! k1 date de depart au format AAAAMMQQ.
! k2 date d'arrivee au format AAAAMMQQ.
! En sortie:
! kec: nombre de jours entre les deux dates.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! Cette routine est utilisable avec des entiers 32 bits,
! si l'ecart entre les deux dates est inferieur a 2**31 jours,
! soit 5879489 ans!...
!
! Au-dela de cette duree, les entiers doivent etre 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! -------------------------------------------------
! Ecart entre les deux dates au jour pres.
! -------------------------------------------------
!


IMPLICIT NONE
INTEGER(KIND=4) :: K1
INTEGER(KIND=4) :: K2
INTEGER(KIND=4) :: KEC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ECARTDJ',0,ZHOOK_HANDLE)
CALL ECARTD(K1,K2,1,KEC)
IF (LHOOK) CALL DR_HOOK('ECARTDJ',1,ZHOOK_HANDLE)
END SUBROUTINE ECARTDJ
SUBROUTINE ECARTDS(CD1,CD2,KEC)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *ecartds* Ecart en secondes entre deux dates.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   97-11, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! cd1 date de depart au format 'AAAAMMQQHHNNSS'.
! cd2 date d'arrivee au format 'AAAAMMQQHHNNSS'.
! En sortie:
! kec: nombre de secondes entre les deux dates.
! --------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! Cette routine est utilisable avec des entiers 32 bits,
! si l'ecart entre les deux dates est inferieur a 2**31 secondes,
! soit 68 ans!...
!
! Au-dela de cette duree, les entiers doivent etre 64 bits.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IAMQ1
INTEGER(KIND=4) :: IAMQ2
INTEGER(KIND=4) :: IH1
INTEGER(KIND=4) :: IH2
INTEGER(KIND=4) :: IM1
INTEGER(KIND=4) :: IM2
INTEGER(KIND=4) :: IS1
INTEGER(KIND=4) :: IS2
INTEGER(KIND=4) :: KEC
INTEGER(KIND=4) :: KECQ
CHARACTER*(*) CD1,CD2
!
! -------------------------------------------------
! On lit les dates sur des entiers.
! -------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ECARTDS',0,ZHOOK_HANDLE)
READ(CD1,FMT='(i8,3i2)') IAMQ1,IH1,IM1,IS1
READ(CD2,FMT='(i8,3i2)') IAMQ2,IH2,IM2,IS2
!
! -------------------------------------------------
! Ecart entre les deux dates au jour pres.
! -------------------------------------------------
!
CALL ECARTD(IAMQ1,IAMQ2,1,KECQ)
!
! -------------------------------------------------
! Ecart en secondes.
! -------------------------------------------------
!
KEC=KECQ*86400+(IH2-IH1)*3600+(IM2-IM1)*60+IS2-IS1
IF (LHOOK) CALL DR_HOOK('ECARTDS',1,ZHOOK_HANDLE)
END SUBROUTINE ECARTDS
SUBROUTINE GREGOD(KDAT,KOPT,KGRE)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *GREGOD *  - Conversion Date > Ecart par rapport a une date fixe.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 92-05-27, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree: kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! En sortie:
! - si kopt=1 : kgre nombre de jours    entre 19000101       et kdat
! - si kopt=2 : kgre nombre d'heures    entre 1900010100     et kdat
! - si kopt=3 : kgre nombre de minutes  entre 190001010000   et kdat
! - si kopt=4 : kgre nombre de secondes entre 19000101000000 et kdat
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IDEBM(12)
INTEGER(KIND=4) :: I0
INTEGER(KIND=4) :: IA100
INTEGER(KIND=4) :: IA4
INTEGER(KIND=4) :: IA400
INTEGER(KIND=4) :: IAAAA
INTEGER(KIND=4) :: IBISSEXT
INTEGER(KIND=4) :: ICONV
INTEGER(KIND=4) :: IFRJOUR
INTEGER(KIND=4) :: II
INTEGER(KIND=4) :: II1
INTEGER(KIND=4) :: IJOURP
INTEGER(KIND=4) :: IMM
INTEGER(KIND=4) :: IN
INTEGER(KIND=4) :: IN1
INTEGER(KIND=4) :: IN2
INTEGER(KIND=4) :: IQQ
INTEGER(KIND=4) :: KDAT
INTEGER(KIND=4) :: KGRE
INTEGER(KIND=4) :: KOPT
DATA IDEBM/0,31,59,90,120,151,181,212,243,273,304,334/
!
! --------------------------------------------------------------------------
! **      1. Calcul du nb de jours separant ki2 du 1er janv 1900
!
! *       1.1 Extraction des quantieme, mois et annee
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GREGOD',0,ZHOOK_HANDLE)
IF(KOPT == 1) THEN
  ! Date de type AAAAMMQQ
  ICONV=1
  IFRJOUR=0
  II=KDAT
ELSEIF(KOPT == 2) THEN
  ! Date de type AAAAMMQQHH
  ICONV=24
  IFRJOUR=MOD(KDAT,100)
  II=KDAT/100
ELSEIF(KOPT == 3) THEN
  ! Date de type AAAAMMQQHHMM
  ICONV=1440
  IFRJOUR=MOD(KDAT,100)
  II=KDAT/100
  IFRJOUR=IFRJOUR+MOD(II,100)*60
  II=II/100
ELSEIF(KOPT == 4) THEN
  ! Date de type AAAAMMQQHHMMSS
  ICONV=86400
  IFRJOUR=MOD(KDAT,100)
  II=KDAT/100
  IFRJOUR=IFRJOUR+MOD(II,100)*60
  II=II/100
  IFRJOUR=IFRJOUR+MOD(II,100)*3600
  II=II/100
ELSE
  ! Cas d'entree erronee de l'utilisateur.
  PRINT*,'GREGODR/ERREUR: argument kopt errone!...'
  PRINT*,KOPT
  STOP 'call abort'
ENDIF
IQQ=II-(II/100)*100
IN=(II-IQQ)/100
IMM=IN-(IN/100)*100
IAAAA=(IN-IMM)/100
! *       1.2 L'annee est-elle bissextile?
! Une annee est bissextile ssi elle est
! (mult de 4 et non mult de 100) ou (mult de 400)
IAAAA=IAAAA
IA400=400*(IAAAA/400)
IA100=100*(IAAAA/100)
IA4=4*(IAAAA/4)
IF((IAAAA == IA400).OR.((IAAAA == IA4).AND.(IAAAA /= IA100)))THEN
  IBISSEXT=1
ELSE
  IBISSEXT=0
ENDIF
IF ((IBISSEXT == 1).AND.(IMM > 2)) THEN
  IJOURP=1
ELSE
  IJOURP=0
ENDIF
! *       1.3 Nombre de jours ecoules depuis le 1er janv
IF(IMM > 12) THEN
  PRINT*,'GREGODR/ERREUR: mois errone.'
  PRINT*,IMM
  STOP 'call abort'
ENDIF
IN2=IDEBM(IMM)+IJOURP+IQQ-1
! *       1.4 Calcul du nb de jours separant les 1er janvier de ii et 1900
I0=1900
IN2=IN2+365*(IAAAA-I0)+INT((IAAAA-1)/4)-INT((I0-1)/4)&
&-INT((IAAAA-1)/100)+INT((I0-1)/100)&
&+INT((IAAAA-1)/400)-INT((I0-1)/400)
! --------------------------------------------------------------------------
! **      2. Calcul du nb de jours separant ii1 du 1er janv 1900
!
! *       2.1 Extraction des quantieme, mois et annee
II1=19000101
II=II1
IQQ=II-(II/100)*100
IN=(II-IQQ)/100
IMM=IN-(IN/100)*100
IAAAA=(IN-IMM)/100
! *       2.2 L'annee est-elle bissextile?
! Une annee est bissextile ssi elle est
! (mult de 4 et non mult de 100) ou (mult de 400)
IAAAA=IAAAA
IA400=400*(IAAAA/400)
IA100=100*(IAAAA/100)
IA4=4*(IAAAA/4)
IF((IAAAA == IA400).OR.((IAAAA == IA4).AND.(IAAAA /= IA100)))THEN
  IBISSEXT=1
ELSE
  IBISSEXT=0
ENDIF
IF ((IBISSEXT == 1).AND.(IMM > 2)) THEN
  IJOURP=1
ELSE
  IJOURP=0
ENDIF
! *       2.3 Nombre de jours ecoules depuis le 1er janv
IN1=IDEBM(IMM)+IJOURP+IQQ-1
! *       2.4 Calcul du nb de jours separant les 1er janvier de ii et 1900
I0=1900
IN1=IN1+365*(IAAAA-I0)+INT((IAAAA-1)/4)-INT((I0-1)/4)&
&-INT((IAAAA-1)/100)+INT((I0-1)/100)&
&+INT((IAAAA-1)/400)-INT((I0-1)/400)
! --------------------------------------------------------------------------
! **      3. Difference in2-in1
KGRE=(IN2-IN1)*ICONV+IFRJOUR
IF (LHOOK) CALL DR_HOOK('GREGOD',1,ZHOOK_HANDLE)
END SUBROUTINE GREGOD
SUBROUTINE GREGOI(KGRE,KOPT,KDAT)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *GREGOI *  - Conversion Ecart par rapport a une date fixe > Date.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 92-05-27, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree: kopt option de precision sur les dates:
! 1 : au jour pres
! 2 : a l'heure pres
! 3 : a la minute pres
! 4 : a la seconde pres
! - si kopt=1 : kgre nombre de jours    entre 19000101       et kdat
! - si kopt=2 : kgre nombre d'heures    entre 1900010100     et kdat
! - si kopt=3 : kgre nombre de minutes  entre 190001010000   et kdat
! - si kopt=4 : kgre nombre de secondes entre 19000101000000 et kdat
! En sortie:
! - si kopt=1 : kdat au format AAAAMMQQ
! - si kopt=2 : kdat au format AAAAMMQQHH
! - si kopt=3 : kdat au format AAAAMMQQHHMM
! - si kopt=4 : kdat au format AAAAMMQQHHMMSS
! --------------------------------------------------------------------------
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ATTENTION A LA PRECISION:
! 1. Vous compilez les entiers sur 32 bits:
! Vous devez alors vous limiter a kopt <= 2.
! 2. Vous compilez les entiers sur 64 bits:
! Vous pouvez utiliser toutes les valeurs de kopt.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


IMPLICIT NONE
INTEGER(KIND=4) :: IJOURS(12)
INTEGER(KIND=4) :: IA100
INTEGER(KIND=4) :: IA4
INTEGER(KIND=4) :: IA400
INTEGER(KIND=4) :: IAAAA
INTEGER(KIND=4) :: IBISSEXT
INTEGER(KIND=4) :: ICONV
INTEGER(KIND=4) :: IDAT
INTEGER(KIND=4) :: IEC
INTEGER(KIND=4) :: IECI
INTEGER(KIND=4) :: IGII2P
INTEGER(KIND=4) :: II2P
INTEGER(KIND=4) :: IMM
INTEGER(KIND=4) :: IMOD
INTEGER(KIND=4) :: IQQ
INTEGER(KIND=4) :: KDAT
INTEGER(KIND=4) :: KGRE
INTEGER(KIND=4) :: KOPT
REAL(KIND=8) :: ZARRDEC
DATA IJOURS/31,28,31,30,31,30,31,31,30,31,30,31/
! --------------------------------------------------------------------------
! **   On determine la date approximative d'arrivee en annees decimales
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GREGOI',0,ZHOOK_HANDLE)
IF(KOPT == 1) THEN
  ! Date de type AAAAMMQQ
  ICONV=1
ELSEIF(KOPT == 2) THEN
  ! Date de type AAAAMMQQHH
  ICONV=24
ELSEIF(KOPT == 3) THEN
  ! Date de type AAAAMMQQHHMM
  ICONV=1440
ELSEIF(KOPT == 4) THEN
  ! Date de type AAAAMMQQHHMMSS
  ICONV=86400
ELSE
  ! Cas d'entree erronee de l'utilisateur.
  PRINT*,'GREGOI/ERREUR: argument kopt errone!...'
  PRINT*,KOPT
  STOP 'call abort'
ENDIF
ZARRDEC=1900.+(REAL(KGRE)/REAL(ICONV)-5.)/365.2425
! --------------------------------------------------------------------------
! **   On determine la date en clair ii2p associee a la date decimale
!
IAAAA=INT(ZARRDEC)
ZARRDEC=12.*(ZARRDEC-REAL(IAAAA))
IMM=INT(ZARRDEC)+1
ZARRDEC=28.*(ZARRDEC-REAL(IMM-1))
IQQ=INT(ZARRDEC)+1
II2P=IQQ+IMM*100+IAAAA*10000
! --------------------------------------------------------------------------
! **   On calcule le nombre de jours separant 19000101 de ii2p
!
CALL GREGOD(II2P,1,IGII2P)
IMOD=MOD(KGRE,ICONV)
IF(IMOD < 0) IMOD=IMOD+ICONV
IEC=(KGRE-IMOD)/ICONV-IGII2P
! --------------------------------------------------------------------------
! **   On avance de iec jours par rapport a ii2p
!
! *       L'annee est-elle bissextile?
! Une annee est bissextile ssi elle est
! (mult de 4 et non mult de 100) ou (mult de 400)
IAAAA=IAAAA
IA400=400*(IAAAA/400)
IA100=100*(IAAAA/100)
IA4=4*(IAAAA/4)
IF((IAAAA == IA400).OR.((IAAAA == IA4).AND.(IAAAA /= IA100)))THEN
  IBISSEXT=1
ELSE
  IBISSEXT=0
ENDIF
! Si oui, 29 jours en fevrier
IF(IBISSEXT == 1) IJOURS(2)=29
! *       Boucle sur les jours
DO IECI=1,IEC
  IQQ=IQQ+1
  IF(IQQ > IJOURS(IMM)) THEN
    IQQ=1
    IMM=IMM+1
  ENDIF
  IF(IMM > 12) THEN
    IMM=1
    IAAAA=IAAAA+1
  ENDIF
ENDDO
! --------------------------------------------------------------------------
! **   On met en forme la date finale
!
IDAT=IQQ+IMM*100+IAAAA*10000
IF(KOPT == 2) THEN
  IMOD=MOD(KGRE,ICONV)
  IF(IMOD < 0) IMOD=IMOD+ICONV
  IDAT=IDAT*100+IMOD
ELSEIF(KOPT == 3) THEN
  IMOD=MOD(KGRE,ICONV)
  IF(IMOD < 0) IMOD=IMOD+ICONV
  IDAT=IDAT*100+IMOD/60
  IMOD=MOD(IMOD,60)
  IDAT=IDAT*100+IMOD
ELSEIF(KOPT == 4) THEN
  IMOD=MOD(KGRE,ICONV)
  IF(IMOD < 0) IMOD=IMOD+ICONV
  IDAT=IDAT*100+IMOD/3600
  IMOD=MOD(IMOD,3600)
  IDAT=IDAT*100+IMOD/60
  IMOD=MOD(IMOD,60)
  IDAT=IDAT*100+IMOD
ENDIF
KDAT=IDAT
IF (LHOOK) CALL DR_HOOK('GREGOI',1,ZHOOK_HANDLE)
END SUBROUTINE GREGOI
FUNCTION IJOURSEM(KDAT)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *IJOURSEM* Jour de la semaine de la date d'entree.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 94-10-31, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------
! En entree:
! kdat1 au format AAAAMMQQ
! En sortie:
! ijour=0 si dimanche, 1 lundi, ..., 6 samedi.
! --------------------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IDATDIM
INTEGER(KIND=4) :: IECART
INTEGER(KIND=4) :: IGRE
INTEGER(KIND=4) :: IGREDIM
INTEGER(KIND=4) :: KDAT
INTEGER(KIND=4) :: IJOURSEM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('IJOURSEM',0,ZHOOK_HANDLE)
CALL GREGOD(KDAT,1,IGRE)
IDATDIM=19941030 ! cette date etait un dimanche.
CALL GREGOD(IDATDIM,1,IGREDIM)
IECART=IGRE-IGREDIM
IJOURSEM=MODULO(IECART,7)
IF (LHOOK) CALL DR_HOOK('IJOURSEM',1,ZHOOK_HANDLE)
END FUNCTION IJOURSEM
SUBROUTINE QQMMAA(KDATD,CDRESD)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------------------
! **** *QQMMAA *  - Conversion d'un entier type AAAAQQMM vers une date en clair.
! --------------------------------------------------------------------------
! Auteur:
! -------
! 92-05-27, J.M. Piriou.
!
! Modifications:
! --------------
!
! --------------------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IAN
INTEGER(KIND=4) :: IGRE
INTEGER(KIND=4) :: ILOC
INTEGER(KIND=4) :: IMM
INTEGER(KIND=4) :: IQQ
INTEGER(KIND=4) :: KDATD
CHARACTER*(*) CDRESD
CHARACTER*03 CLJOUR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('QQMMAA',0,ZHOOK_HANDLE)
IQQ=MOD(KDATD,100)
ILOC=KDATD/100
IMM=MOD(ILOC,100)
IAN=ILOC/100
CALL GREGOD(KDATD,1,IGRE)
IGRE=MOD(IGRE,7)
IF(IGRE == 0) THEN
  CLJOUR='Lun'
ELSEIF(IGRE == 1) THEN
  CLJOUR='Mar'
ELSEIF(IGRE == 2) THEN
  CLJOUR='Mer'
ELSEIF(IGRE == 3) THEN
  CLJOUR='Jeu'
ELSEIF(IGRE == 4) THEN
  CLJOUR='Ven'
ELSEIF(IGRE == 5) THEN
  CLJOUR='Sam'
ELSEIF(IGRE == 6) THEN
  CLJOUR='Dim'
ENDIF
WRITE(CDRESD,FMT='(a3,a1,i2,a1,i2.2,a1,i4.4)')&
&CLJOUR,' ',IQQ,'.',IMM,'.',IAN
IF (LHOOK) CALL DR_HOOK('QQMMAA',1,ZHOOK_HANDLE)
END SUBROUTINE QQMMAA
SUBROUTINE QUANT(KDATE,KQUANT)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
! --------------------------------------------------------------
! **** *quant* Quantieme de l'annee d'une date donnee.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   1999-04, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! kdate date au format AAAAMMQQ.
! En sortie:
! quantieme de l'annee (1 le 1er janvier, 32 le 1er fevrier, etc...)
! --------------------------------------------------------------


IMPLICIT NONE
INTEGER(KIND=4) :: IBASE
INTEGER(KIND=4) :: IEC
INTEGER(KIND=4) :: KDATE
INTEGER(KIND=4) :: KQUANT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('QUANT',0,ZHOOK_HANDLE)
IBASE=10000*(KDATE/10000)+0101 ! 1er janvier de l'annee courante.
CALL ECARTDJ(IBASE,KDATE,IEC)
KQUANT=IEC+1
IF (LHOOK) CALL DR_HOOK('QUANT',1,ZHOOK_HANDLE)
END SUBROUTINE QUANT
