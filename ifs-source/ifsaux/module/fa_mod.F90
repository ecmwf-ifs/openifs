! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE FA_MOD
! Feb-2013 P. Marguinaud Define JNGEOM & JNEXPL
!                        Fix bug in cadre allocation
! Jan-2011 P. Marguinaud Interface to thread-safe FA
! Sep-2012 P. Marguinaud DrHook

USE PARKIND1, ONLY : JPIM, JPRB, JPRD, JPIA
USE YOMHOOK, ONLY : LHOOK, DR_HOOK, JPHOOK
USE LFIMOD, ONLY : LFICOM

USE LFI_PRECISION
IMPLICIT NONE

! Date index

INTEGER (KIND=JPLIKB), PARAMETER :: &
& JD_YEA =  1, JD_MON =  2, JD_DAY =  3, &  !  Base year, month, day
& JD_HOU =  4, JD_MIN =  5, JD_TUN =  6, &  !  Base hour, minute, term unit (1=hour, 2=day)
& JD_THO =  7, JD_GR8 =  8, JD_IAN =  9, &  !  Term hour or day, ?, forecast/initialised analysis flag
& JD_CU1 = 10, JD_CU2 = 11,              &  !  Time of event1 (hour or day), time of event2 (hour or day)
& JD_DEX = 12,              JD_SEM = 14, &  !  Extended date flag,                       base seconds
& JD_SET = 15, JD_CE1 = 16, JD_CE2 = 17, &  !  Term in seconds, seconds of event1, seconds of event2
& JD_TST = 18, JD_FMT = 19,              &  !  Time step, time format
& JD_SIZ = 22

INTEGER (KIND=JPLIKB), PARAMETER :: JPPRCM = JPDBLD / JPDBLR

INTEGER (KIND=JPLIKB), PARAMETER :: JNGEOM =  18_JPLIKB
INTEGER (KIND=JPLIKB), PARAMETER :: JNEXPL =   8_JPLIKB

INTEGER (KIND=JPLIKB), PARAMETER :: JPNIIL =-999_JPLIKB
INTEGER (KIND=JPLIKB), PARAMETER :: JPXNOM =  16_JPLIKB
INTEGER (KIND=JPLIKB), PARAMETER :: JPXPRF =   8_JPLIKB
INTEGER (KIND=JPLIKB), PARAMETER :: JPXSUF = JPXNOM+JPXPRF 

INTEGER (KIND=JPLIKB), PARAMETER :: NUNDEF = JPNIIL
REAL (KIND=JPDBLR),    PARAMETER :: XUNDEF = -99._JPDBLR
LOGICAL,               PARAMETER :: LUNDEF = .FALSE.
CHARACTER,             PARAMETER :: CUNDEF = CHAR (0)
!****
!           !-----------------------------------------------!
!           ! Sous-programme du Logiciel de Fichiers ARPEGE !
!           !-----------------------------------------------!
!
!       - Version originale du logiciel: Mars 1990, auteur:
!           Jean CLOCHARD, Meteorologie Nationale, FRANCE.
!
!       - Juin 90: Ajout du "type de transformation horizontale"
!         aux elements definissant un cadre
!         (1=transformation de Frank-Schmidt "pure", 2=ARPEGE complete),
!         et ajout de la notion d'"Identificateur de Fichier".
!
!**
!---------------- VARIABLES "SIMPLES" GLOBALES -------------------------
!
!     NFIOUV = Nombre de fichiers ouverts simultanement
!     NCADEF =   "    "  cadres definis        "
!     LFAMUL = Option "utilisation du logiciel en mode multi-taches"
!     LFAMOP = Option "mode mise au point du logiciel"
!     VRGLAS = Verrou global du lociciel (en mode multi-taches)
!     NIMSGA = Niveau global de messagerie du logiciel
!     NRFAGA = Niveau de filtrage global des erreurs fatales detectees ?
!     NBIPDG = Nombre de bits IMPLICITE par valeur point-de-grille
!     NBICSP =    "   "   "   IMPLICITE par coefficient spectral
!     NSTROI = Sous-troncature non compactee IMPLICITE (coef. spectraux)
!     NPUILA = "Puissance de laplacien" IMPLICITE pour coeffi. spectraux
!     NMIDPL = Degre de modulation IMPLICITE de puissance de laplacien
!     NIGRIB = Niveau IMPLICITE de codage GRIB
!              (  0 ==> pas de codage, 1 ==> "standard", 2 => "Arpege" )
!              ( -1 ==> pas de codage (rangt coeff spec selon modele)  )
!              (  3 ==> "GRIBEX"      (rangt coeff spec selon modele)  )
!     NCPCAD = Nombre maximum de caracteres possibles par nom de cadre.
!     SPSMIN = Pression sol minimum   "     "    "     "      "    "
!     SPSMAX = Pression sol maximum   "     "    "     "      "    "
!     MPRESX = Pression maximum pour designer un niveau isobare
!     NBIMAC = Nombre de bits par mot machine
!     NBIMAX = Nombre maximum de bits par valeur elementaire compactee
!     NTYPTX = Type maximum reconnu de transformation horizontale
!     LIGARD = Option de conservation des cadres definis dynamiquement
!
!        Ci-dessous, limites en termes de dimensions reglables au niveau
!     usager, initialisees aux valeurs maximum du logiciel:
!
!     NXNIVV = Nombre maximum de niveaux verticaux au niveau usager
!     NXTRON = Troncature maximum                  "    "      "
!     NXLATI = Nombre maximum de latitudes de pole a pole "    "
!     NXLONG = Nombre maximum de longitudes par parallele "    "
!
!     CHAINC = Impression de substitut aux variables non reconnues
!**
!------------------------- TABLEAUX GLOBAUX ----------------------------
!
!     NULIND = Table d'indirection pour les autres tables "fichiers"
!     NCAIND =   "   "      "       "    "     "      "   "cadres"
!**
!--------- DESCRIPTIF DES ELEMENTS CONCERNANT UNE UNITE LOGIQUE --------
!
!     NULOGI = Numero de l'unite logique
!     NUCADR = Rang du cadre auquel le fichier est rattache
!     CIDENT = Identification (utilisateur) du fichier
!     LNOMME = Vrai si l'unite logique a un nom utilisateur (BOF...)
!     NIVOMS = Niveau de la messagerie defini par l'utilisateur
!     LERRFA = Vrai si toute erreur detectee doit etre fatale
!     LCREAF = Vrai si le fichier est en mode creation
!     NBFPDG = Nombre de bits EFFECTIF par valeur point-de-grille
!     NBFCSP =    "   "   "   EFFECTIF par coefficient spectral
!     NPUFLA = "Puissance de laplacien" EFFECTIVE pour coeffi. spectraux
!     NSTROF = Sous-troncature non compactee EFFECTIVE (coef. spectraux)
!     NMFDPL = Degre de modulation EFFECTIF de la puissance de laplacien
!     NFGRIB = Niveau EFFECTIF de codage GRIB
!     VRFICH = Verrou du fichier (en mode multi-taches)
!
!     MADATE = Contenu de l'article 'DATE'
!
!     CTNPRF = Prefixes de champ (re)connus.
!     NIVDSC = Descripteurs associes aux prefixes (re)connus.
!     NRASHO = Nombre de champs de coeff. spectraux ARPEGE ranges
!              horizontalement (donc a l'inverse du modele) detectes au
!              cours des lectures/ecritures
!     NRASVE = Nombre de champs de coeff. spectraux ARPEGE ranges
!              verticalement (donc comme le modele), detectes au cours
!              des lectures/ecritures
!
!**
!----------- DESCRIPTIF DES ELEMENTS CONSTITUANT UN CADRE --------------
!
!                     Elements DEFINISSANT un cadre:
!
!     CNOMCA = Nom du cadre (conventionnel; non ecrit sur le fichier)
!     NLCCAD = Longueur en caracteres du nom de cadre
!     NULCAD = Nombre d'unites logiques rattachees au cadre
!     MTRONC = Troncature associee aux donnees contenues dans le fichier
!     NNIVER = Nombre de niveaux verticaux
!     NLATIT =   "    "  latitudes de pole a pole
!     NXLOPA =   "    maximum de longitudes par cercle de latitude
!     NTYPTR = Type de transformation horizontale
!     SSLAPO = Sinus   de la latitude  du pole d'interet
!     SCLOPO = Cosinus  "   "     "     "    "      "
!     SSLOPO = Sinus    "   "     "     "    "      "
!     SCODIL = Coefficient de dilatation
!     SPREFE = Pression de reference (facteur multiplicatif de la premi-
!              ere fonction de la coordonnee hybride)
!     NLOPAR = Nombre de longitudes par parallele de l'hemisphere "nord"
!     NOZPAR =    "   d'onde zonal maximum "    "  "     "        "
!     SINLAT = Sinus des latitudes de la grille (*hemisphere nord seul*)
!     SFOHYB = Valeurs des fonctions "A" et "B" de la coordonnee hybride
!              (Definies aux LIMITES DE COUCHES)
!
!                Autres elements rattaches a un cadre:
!
!     NVAPDG = Nombre total de points de grille dans un champ horizontal
!     NGARDE = Option de conservation lors de la fermeture du dernier
!              fichier s'y rattachant. ( 0=Non, 1=selon LIGARD, 2=Oui )
!     LIMLAM = IMPLICIT SWITCH FOR SETTING OF GLOBAL/LAM CASE
!     NSFLAM = TOTAL NUMBER OF SPECTRAL COEFFICIENTS FOR LAM CASE
!

TYPE FACADR
  INTEGER (KIND=JPLIKB)          :: MTRONC         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NNIVER         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NLATIT         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NXLOPA         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NULCAD         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NLCCAD         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NGARDE         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NVAPDG         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NTYPTR         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NSFLAM         = NUNDEF
  REAL (KIND=JPDBLR)             :: SSLAPO         = XUNDEF
  REAL (KIND=JPDBLR)             :: SCLOPO         = XUNDEF
  REAL (KIND=JPDBLR)             :: SPREFE         = XUNDEF
  REAL (KIND=JPDBLR)             :: SSLOPO         = XUNDEF
  LOGICAL                        :: LIMLAM         = LUNDEF
  CHARACTER*(JPXNOM)             :: CNOMCA         = CUNDEF
  REAL (KIND=JPDBLR)             :: SCODIL         = XUNDEF
  REAL (KIND=JPDBLR),    POINTER :: SINLAT (:)     => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NLOPAR (:)     => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NOZPAR (:)     => NULL ()
  REAL (KIND=JPDBLR),    POINTER :: SFOHYB (:,:)   => NULL ()

!      EN-TETE POUR GRIBEX (SECTIONS 1 ET 2) 

  INTEGER (KIND=JPLIKB), POINTER :: NSEC2SP(:)     => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NSEC2LL(:)     => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NSEC2GG(:)     => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NSEC2LA(:)     => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NSEC2AL(:)     => NULL ()
  REAL (KIND=JPDBLR),    POINTER :: XSEC2(:)       => NULL ()
  LOGICAL                        :: LISEC2         = .TRUE.

  INTEGER (KIND=JPLIKB), POINTER :: NOMPAR(:)      => NULL ()

  INTEGER (KIND=JPLIKB)          :: NSEFRE         = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NSMAX          = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NMSMAX         = NUNDEF
  INTEGER (KIND=JPLIKB), POINTER :: NCPL4M  (:)    => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NISMAX  (:)    => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NISNAX  (:)    => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NDIM0GG (:)    => NULL ()
  INTEGER (KIND=JPLIKB)          :: IADDPK         = 0_JPLIKB

END TYPE FACADR

TYPE FAFICH
  INTEGER (KIND=JPLIKB)          :: NULOGI       = JPNIIL   ! From farine.F90
  INTEGER (KIND=JPIA)            :: NFILEP       = 0_JPIA   ! Auxiliary FILE pointer
  INTEGER (KIND=JPLIKB)          :: NOFFST       = 0_JPLIKB ! Auxiliary FILE offset
  INTEGER (KIND=JPLIKB)          :: NUCADR       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NIVOMS       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NBFPDG       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NBFCSP       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NPUFLA       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NFGRIB       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NSTROF       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NMFDPL       = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NRASHO       = 0_JPLIKB ! From farine.F90
  INTEGER (KIND=JPLIKB)          :: NRASVE       = 0_JPLIKB ! From farine.F90
  INTEGER (KIND=JPLIKB), POINTER :: MADATE (:)   => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: MADATX (:)   => NULL ()
  LOGICAL                        :: LERRFA       = LUNDEF
  REAL (KIND=JPDBLR)             :: VRFICH       = XUNDEF
  LOGICAL                        :: LCREAF       = LUNDEF
  CHARACTER*(JPXNOM)             :: CIDENT       = CUNDEF
  LOGICAL                        :: LIFLAP       = .FALSE.  ! From farine.F90
  LOGICAL                        :: LNOMME       = LUNDEF
  INTEGER (KIND=JPLIKB)          :: NSEC1 (2:21) 
  INTEGER (KIND=JPLIKB), POINTER :: NSC2ALF(:)   => NULL ()
  LOGICAL                        :: LISEC1       = .TRUE.   ! From farine.F90
  LOGICAL                        :: LISC2F       = .TRUE.   ! From farine.F90
  INTEGER (KIND=JPLIKB)          :: NCOGRIF(12) 
  REAL (KIND=JPDBLR),    POINTER :: FLAP1D(:)    => NULL ()
  REAL (KIND=JPDBLR),    POINTER :: FLAP1DA(:)   => NULL ()
  INTEGER (KIND=JPLIKB)          :: NCPLSIZE     = NUNDEF
  INTEGER (KIND=JPLIKB)          :: NCPLBITS     = NUNDEF
  INTEGER (KIND=JPLIKB)          :: IOPTGRSX2O   = NUNDEF
  INTEGER (KIND=JPLIKB)          :: IOPTGRSN2O   = NUNDEF
  CHARACTER (LEN=64)             :: CMODEL       = ''
  INTEGER (KIND=JPLIKB)          :: NIDCEN       = 85_JPLIKB
END TYPE FAFICH

TYPE FAGR1TAB

  CHARACTER(LEN=JPXPRF) :: CIPREF = CUNDEF
  CHARACTER(LEN=JPXSUF) :: CISUFF = CUNDEF
  INTEGER (KIND=JPLIKB) :: NCODPA(8) = (/ NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF /)
  LOGICAL               :: LFNIVA = .FALSE.  
  REAL (KIND=JPDBLR)    :: FMULTI = XUNDEF  ! Facteur multiplicatif avant encodage
  LOGICAL               :: LMULTI = .FALSE. ! Appliquer FMULTI ou non
  
END TYPE FAGR1TAB

TYPE FA_COM

  TYPE(LFICOM), POINTER :: LFI => NULL ()

  INTEGER (KIND=JPLIKB), POINTER :: NULIND (:) => NULL ()
  INTEGER (KIND=JPLIKB), POINTER :: NCAIND (:) => NULL ()

  INTEGER (KIND=JPLIKB) NFIOUV, NCADEF, NIMSGA, NRFAGA
  INTEGER (KIND=JPLIKB) NBIPDG, NBICSP, NPUILA
  INTEGER (KIND=JPLIKB) NIGRIB, NCPCAD, NSTROI, NMIDPL
  INTEGER (KIND=JPLIKB) NBIMAC, NBIMAX, MPRESX
  INTEGER (KIND=JPLIKB) NXNIVV, NXTRON, NXLATI, NXLONG, NTYPTX
!
  INTEGER (KIND=JPLIKB), POINTER :: NIVDSC (:,:) => NULL ()
!
  REAL (KIND=JPDBLR) SPSMIN, SPSMAX, VRGLAS
!
!
  LOGICAL LFAMUL, LFAMOP, LIGARD
!
  CHARACTER*(JPXNOM) CHAINC
  CHARACTER(LEN=8), POINTER :: CTNPRF (:) => NULL ()
!
!**
!----- DESCRIPTION DES "PARAMETER" DU LOGICIEL DE "FICHIER ARPEGE" -----
!
!     JPNXFA = Nombre maximum de fichiers ouverts "simultanement"
!     JPNXCA =    "      "    "  cadres definissables "simultanement"
!     JPXNIV =    "      "    "  niveaux verticaux (champs d'altitude)
!     JPXTRO = Troncature maximum gerable
!     JPXLAT = Nombre maximum de latitudes de pole a pole
!     JPXLON = Nombre maximum de longitudes par parallele
!     JPLDAT = Longueur de l'article 'DATE', en mots
!     JPLB1P = Longueur du tableau "Bloc 1" pour sous-programmes GRIB
!     JPLB2P = Longueur du tableau "Bloc 2" pour sous-programmes GRIB
!     JPNIIL = Code "valeur absente" du logiciel pour les entiers
!     JPXCSP = Dimension maxi d'un champ en coefficients spectraux
!     JPXPDG = Dimension maxi d'un champ en points de grille
!     JPXCHA = Dimension maxi d'un champ ( maximum de JPXCSP et JPXPDG )
!     JPXPAH = Nombre maximum de latitudes par hemisphere
!     JPXIND = DIMENSIONING OF NOZPAR()
!     JPXGEO = DIMENSIONING OF SINLAT()
!     JPNVER = Numero de version du logiciel (qui est le contenu de
!              l'article dont le nom est l'identificateur du fichier)
!     JPUILA = Puissance de laplacien maximum pour laquelle les tableaux
!              servant a calculer laplacien et inverse sont precalcules
!     JPXNOM = Nombre maximum de caracteres par NOM d'article LFI.
!     JPXPRF =   "       "    "      "      par PReFixe de champ.
!     JPXSUF = JPXPRF+JPXNOM.
!     JPTNIV = Nombre de types de niveaux verticaux (re)connus.
!     CPDATE = Nom de l'article DATE
!
!         Noms des articles contenant les differentes parties du CADRE:
!
!     CPCADI = "Dimensions" (MTRONC, NNIVER, NLATIT, NXLOPA)
!     CPCAFS = Parametres de la transformation ARPEGE
!              (SSLAPO, SCLOPO, SSLOPO, SCODIL)
!     CPCARP = Tableaux lies a la reduction des points pres des poles
!     CPCASL = Tableau des sinus des latitudes
!     CPCACH = Valeurs des fonctions "A" et "B" de la coordonnee hybride
!     JPCADI et JPCAFS sont les longueurs des 2 premiers de ces articles
!
  INTEGER (KIND=JPLIKB) JPNXFA, JPNXCA, JPLDAT
  INTEGER (KIND=JPLIKB) JPXNIV, JPXTRO, JPXLAT
  INTEGER (KIND=JPLIKB) JPUILA, JPXAU1, JPXLON
  INTEGER (KIND=JPLIKB) JPXAU2, JPXPAH, JPXIND, JPXGEO
  INTEGER (KIND=JPLIKB) JPXCSP, JPXCHA, JPLB1P
  INTEGER (KIND=JPLIKB) JPLB2P, JPCADI, JPCAFS, JPNVER
  INTEGER (KIND=JPLIKB) JPXPDG, JPXNOM, JPXPRF, JPXSUF, JPTNIV
!
  CHARACTER CPCADI*(16), CPCAFS*(16), CPCARP*(16), CPCACH*(16)
  CHARACTER CPCASL*(16), CPDATE*(16), CPDATX*(16)
!


!*
!      FAMODU - MODULE POUR LE LOGICIEL FA

!        D. PARADIS       METEO FRANCE     21/7/00

!  ----------------------------------------------------------------------
!  1 - COEFFICIENTS MULTIPLICATIFS DES COEFF SPECTRAUX
!
!     XLAP1D  = Coefficients elementaires pour "laplacien" [n*(n+1)],
!              ainsi que leurs inverses.
!     XLAP2D  = Cf. XLAP1D, mais etale dans un champ complet, et calcules
!              pour les puissances 1 a JPUILA, ainsi que les inverses.
!     FLAP1D  = Tableau des "puissances de laplacien" pour l'ecriture
!                     (n*(n+1))**NPUFLA()
!              ( Ce dernier tableau n'est calcule que pour une
!                puissance de laplacien non nulle )
!     XLAP1DA = ALADIN version de XLAP1D
!     XLAP2DA = ALADIN version de XLAP2D
!     FLAP1DA = ALADIN version de FLAP1D
!     LIXLAP  = .T. s'il faut initialiser XLAPxDx
!     LIFLAP() = .T. s'il faut initialiser FLAP1Dx()
!
!  2 - EN-TETE POUR GRIBEX (SECTIONS 1 ET 2)
!
!     JPSEC1   = taille du tableau contenant les elements de la section 1
!     JPSEC2   = taille du tableau contenant les elements de la section 2
!     JPSEC4   = taille du tableau contenant les elements de la section 4
!     NSEC1    = tableau contenant les elts 2:21 de la section 1 de GRIBEX
!     LISEC1() = .T. s'il faut initialiser le tableau NSEC1 ci-dessus
!
!      Tableaux contenant des elements de la section 2 de GRIBEX:
!     NSEC2SP  = cas de la representation spectrale ARPEGE avec niveaux modele
!     NSEC2GG  = cas de la grille de Gauss
!     NSEC2LL  = cas de la grille latitude-longitude
!     NSEC2LA  = cas de la grille Lambert conforme (type general de Aladin)
!     NSEC2AL  = cas de la representation spectrale Aladin deguisee lat-lon
!     NSC2ALF  = supplement a NSEC2AL dependant du fichier
!     XSEC2    = coordonnee verticale
!     LISEC2() = .T. s'il faut initialiser les tableaux ci-dessus sauf NSC2ALF
!     LISC2F() = .T. s'il faut initialiser NSC2ALF
!
!  3 - PARAMETRES DEFINISSANT LE CODAGE GRIBEX
!
!     NCODGRI  = tableau contenant les parametres de codage implicites
!     NCOGRIF  = tableau contenant les parametres de codage pour chaque fichier
!     CIPREF   = tableau contenant les prefixes des noms des champs connus
!     CISUFF   = tableau contenant les suffixes des noms des champs connus
!     NCODPA   = tableau contenant les 6 descripteurs GRIB des champs connus
!                NCODPA(1)= numero de version de la table de code parametre,
!                           KSEC1(1)
!                NCODPA(2)= indicateur de parametre, KSEC1(6)
!                NCODPA(3)= indicateur de type de niveau, KSEC1(7)
!                NCODPA(4)= niveau (premier niveau de la couche), KSEC1(8)
!                NCODPA(5)= deuxieme niveau de la couche, KSEC1(9)
!                NCODPA(6)= indicateur de type de champ, KSEC1(18)
!                           (0 sauf si min/max dans le temps => 2,
!                                ou si cumul dans le temps   => 4,
!                                ou cumul depuis le debut    => 8 )
!                NCODPA(7)= facteur d'echelle decimal, KSEC1(23)
!     JPXPAR   = nombre maximal de champs connus pour CISUFF, CIPREF et NCODPA
!     NBPARC   = nombre de champs connus pour CISUFF, CIPREF et NCODPA
!
!  4 - TABLEAU POUR CODER LES COEF SPECTRAUX ALADIN AVEC GRIBEX
!
!     NOMPAR   = tableau decrivant la position des coeff spectraux (CSP) ALADIN
!                en fonction du nombre d'onde zonal dans un tableau de donnees
!                stocke dans FA avec la methode de codage 3 ou -1 (rangt modele).
!                NOMPAR est l'equivalent de NOZPAR mais pour un rangement vertical
!                des CSP: NOMPAR(1)=NSMAX, NOMPAR(2)=NMSMAX (ces 2 valeurs de
!                NOMPAR sont speciales commme pour NOZPAR) et pour un JM donne,
!                compris entre 0 et NMSMAX, NOMPAR(2*JM+3) donne l'indice du
!                premier CSP associe a JM qui est contenu dans un champ spectral
!                et NOMPAR(2*JM+4) donne l'indice du dernier CSP associe a JM
!                qui est contenu dans un champ spectral.
!                
!  ----------------------------------------------------------------------

!  1 - COEFFICIENTS MULTIPLICATIFS DES COEFF SPECTRAUX

  REAL (KIND=JPDBLR), POINTER :: XLAP1D(:,:)     => NULL ()
  REAL (KIND=JPDBLR), POINTER :: XLAP1DA(:,:)    => NULL ()
  REAL (KIND=JPDBLR), POINTER :: XLAP2D(:,:,:)   => NULL ()
  REAL (KIND=JPDBLR), POINTER :: XLAP2DA(:,:,:)  => NULL ()

  LOGICAL LIXLAP
 
!  2 - EN-TETE POUR GRIBEX (SECTIONS 1 ET 2)

  INTEGER (KIND=JPLIKB) :: JPSEC1, JPSEC2
  INTEGER (KIND=JPLIKB) :: JPSEC4

  INTEGER (KIND=JPLIKB) :: IOPTGRSX2O   = NUNDEF
  INTEGER (KIND=JPLIKB) :: IOPTGRSN2O   = NUNDEF
 
!  3 - PARAMETRES DEFINISSANT LE CODAGE GRIBEX

  INTEGER (KIND=JPLIKB) :: NCODGRI(12) = (/ NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF, &
                                          & NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF, NUNDEF /)
  INTEGER (KIND=JPLIKB) JPXPAR
  INTEGER (KIND=JPLIKB) NBPARC
  INTEGER (KIND=JPLIKB) :: NIDCEN = 85_JPLIKB
  
  
  TYPE (FAGR1TAB), POINTER :: YGR1TAB (:) => NULL ()

!  4 - TABLEAU POUR CODER LES COEF SPECTRAUX ALADIN AVEC GRIBEX


  LOGICAL :: FACADE_LLPREA = .TRUE.
  LOGICAL :: FACAGE_LLPREA = .TRUE.
  LOGICAL :: FACIES_LLPREA = .TRUE.
  LOGICAL :: FACTUM_LLPREA = .TRUE.
  LOGICAL :: FAGIOT_LLPREA = .TRUE.
  LOGICAL :: FALIMU_LLPREA = .TRUE.
  LOGICAL :: FAMISO_LLPREA = .TRUE.
  LOGICAL :: FANERG_LLPREA = .TRUE.
  LOGICAL :: FANMSG_LLPREA = .TRUE.
  LOGICAL :: FANUCA_LLPREA = .TRUE.
  LOGICAL :: FANUMU_LLPREA = .TRUE.
  LOGICAL :: FAREGI_LLPREA = .TRUE.
  LOGICAL :: FARFLU_LLPREA = .TRUE.
  LOGICAL :: FARINE_LLPREA = .TRUE.
  LOGICAL :: FAVORI_LLPREA = .TRUE.
  LOGICAL :: FAXION_LLPREA = .TRUE.
  LOGICAL :: FARINE_LLDEFM = .FALSE.
  INTEGER (KIND=JPLIKB) :: FAXION_ISCALX
  REAL (KIND=JPDBLR) FAXION_ZEPSIL

  INTEGER (KIND=JPLIKB) :: NULOUT = 0
  LOGICAL :: LOPENMP = .TRUE.

  INTEGER (KIND=JPLIKB) :: JPLSPX = 6
  INTEGER (KIND=JPLIKB) :: JPLMES = 1024

  TYPE(FACADR), POINTER :: CADRE   (:) => NULL ()
  TYPE(FAFICH), POINTER :: FICHIER (:) => NULL ()

END TYPE FA_COM

INTEGER, SAVE :: NGRIB2_GLO_SH = NUNDEF  ! Spherical harmonics
INTEGER, SAVE :: NGRIB2_GLO_GP = NUNDEF  ! Gaussian grid
INTEGER, SAVE :: NGRIB2_LAM_GP = NUNDEF  ! Grid-point LAM 
INTEGER, SAVE :: NGRIB2_LAM_BF = NUNDEF  ! LAM bi-Fourier
INTEGER, SAVE :: NGRIB1_LATLON = NUNDEF  
INTEGER, SAVE :: NGRIB2_LATLON = NUNDEF 
LOGICAL, SAVE :: LGRIB2_LAM_EX = .FALSE. ! Extension zone metadata are available
LOGICAL, SAVE :: LGRIB2_LAM_BF = .FALSE. ! LAM bi-Fourier is available
LOGICAL, SAVE :: LGRIB2_INIT   = .FALSE.

TYPE(FA_COM), SAVE, TARGET :: FA_COM_DEFAULT
LOGICAL, SAVE :: FA_COM_DEFAULT_INIT = .FALSE.

CONTAINS

SUBROUTINE NEW_CADRE (CA, KTYPTR, KPXLAT, KPXTRO, KPXNIV)
TYPE (FACADR) :: CA
INTEGER (KIND=JPLIKB), INTENT (IN) :: KTYPTR, KPXLAT, KPXTRO, KPXNIV
LOGICAL :: LLMLAM

INTEGER (KIND=JPLIKB) :: &
&      IPXAU1, IPXLON, IPXAU2, IPXPAH, &
&      IPXIND, IPXGEO, IPXCSP, IPXPDG, &
&      IPXCHA
INTEGER (KIND=JPLIKB) :: INPAHE
INTEGER (KIND=JPLIKB) :: JM, JN, IPOS

LLMLAM = KTYPTR .LE. 0

CALL CPARAMS (KPXLAT, KPXTRO, IPXAU1,  &
&      IPXLON, IPXAU2, IPXPAH, IPXIND, &
&      IPXGEO, IPXCSP, IPXPDG, IPXCHA)
   
IF (.NOT. LLMLAM) THEN
  INPAHE=(1+KPXLAT)/2
  ALLOCATE (                     &
 &    CA%NLOPAR (INPAHE),        &
 &    CA%NOZPAR (INPAHE),        &
 &    CA%SINLAT (INPAHE),        &
 &    CA%NOMPAR (2*KPXTRO+4))
ELSE
  ALLOCATE (                     &
 &    CA%NLOPAR (JNEXPL),        &
 &    CA%NOZPAR (2*KPXTRO+4),    &
 &    CA%SINLAT (JNGEOM),        &
 &    CA%NOMPAR (2*MAX (-KTYPTR, KPXTRO)+4))
ENDIF

ALLOCATE (                       &
&  CA%SFOHYB  (2,0:KPXNIV),      &
&  CA%NSEC2SP (22),              &
&  CA%NSEC2LL (22),              &
&  CA%NSEC2GG (22+KPXLAT),       &
&  CA%NSEC2LA (22),              &
&  CA%NSEC2AL (22),              &
&  CA%XSEC2   (10+2*(KPXNIV+1)))

IF (LLMLAM) THEN

  CA%NSMAX  = KPXTRO
  CA%NMSMAX = - KTYPTR

  ALLOCATE (CA%NISNAX (0:CA%NMSMAX), CA%NISMAX (0:CA%NSMAX), &
          & CA%NCPL4M (0:CA%NMSMAX), CA%NDIM0GG (0:CA%NMSMAX))

  CALL ELLIPS64 (CA%NSMAX, CA%NMSMAX, CA%NISNAX, CA%NISMAX)

  CA%NSEFRE = 0
  DO JM = 0, CA%NMSMAX
    CA%NSEFRE = CA%NSEFRE + 4*(CA%NISNAX(JM)+1)
  ENDDO

  DO JM = 0, CA%NMSMAX
    CA%NCPL4M(JM) = 4*(CA%NISNAX(JM)+1)
  ENDDO

  IPOS = 1
  DO JM = 0, CA%NMSMAX
    CA%NDIM0GG (JM) = IPOS
    IPOS = IPOS + CA%NCPL4M (JM)
  ENDDO

ELSE
  CA%NSMAX  = KPXTRO
  CA%NMSMAX = KPXTRO
  CA%NSEFRE = (CA%NSMAX+1)*(CA%NSMAX+1)
  ALLOCATE (CA%NDIM0GG (0:CA%NSMAX))
  
  IPOS = 1
  DO JN = 0, CA%NSMAX
     CA%NDIM0GG (JN) = IPOS
     IPOS = IPOS + (CA%NSMAX+1-JN) * 2
  ENDDO

ENDIF

END SUBROUTINE NEW_CADRE

SUBROUTINE FREE_CADRE (CA)

TYPE (FACADR) :: CA

TYPE (FACADR) :: CADUM

IF (ASSOCIATED (CA%NLOPAR )) DEALLOCATE (CA%NLOPAR )
IF (ASSOCIATED (CA%NOZPAR )) DEALLOCATE (CA%NOZPAR )
IF (ASSOCIATED (CA%SINLAT )) DEALLOCATE (CA%SINLAT )
IF (ASSOCIATED (CA%SFOHYB )) DEALLOCATE (CA%SFOHYB )
IF (ASSOCIATED (CA%NSEC2SP)) DEALLOCATE (CA%NSEC2SP)
IF (ASSOCIATED (CA%NSEC2LL)) DEALLOCATE (CA%NSEC2LL)
IF (ASSOCIATED (CA%NSEC2GG)) DEALLOCATE (CA%NSEC2GG)
IF (ASSOCIATED (CA%NSEC2LA)) DEALLOCATE (CA%NSEC2LA)
IF (ASSOCIATED (CA%NSEC2AL)) DEALLOCATE (CA%NSEC2AL)
IF (ASSOCIATED (CA%XSEC2  )) DEALLOCATE (CA%XSEC2  )
IF (ASSOCIATED (CA%NOMPAR )) DEALLOCATE (CA%NOMPAR )

IF (ASSOCIATED (CA%NCPL4M )) DEALLOCATE (CA%NCPL4M )
IF (ASSOCIATED (CA%NISMAX )) DEALLOCATE (CA%NISMAX )
IF (ASSOCIATED (CA%NISNAX )) DEALLOCATE (CA%NISNAX )
IF (ASSOCIATED (CA%NDIM0GG)) DEALLOCATE (CA%NDIM0GG)

NULLIFY (CA%NLOPAR )
NULLIFY (CA%NOZPAR )
NULLIFY (CA%SINLAT )
NULLIFY (CA%SFOHYB )
NULLIFY (CA%NSEC2SP)
NULLIFY (CA%NSEC2LL)
NULLIFY (CA%NSEC2GG)
NULLIFY (CA%NSEC2LA)
NULLIFY (CA%NSEC2AL)
NULLIFY (CA%XSEC2  )
NULLIFY (CA%NOMPAR )

NULLIFY (CA%NCPL4M )
NULLIFY (CA%NISMAX )
NULLIFY (CA%NISNAX )
NULLIFY (CA%NDIM0GG)

CA = CADUM

END SUBROUTINE FREE_CADRE

SUBROUTINE NEW_FICHIER (FA, FI, KPLDAT, KPXTRO, KTYPTR)

TYPE (FA_COM) :: FA
TYPE (FAFICH) :: FI
INTEGER (KIND=JPLIKB), INTENT (IN) :: KPLDAT, KPXTRO, KTYPTR
LOGICAL :: LLMLAM

LLMLAM = KTYPTR .LE. 0

ALLOCATE (FI%MADATE (KPLDAT), FI%MADATX (KPLDAT))

FI%MADATE  = NUNDEF
FI%MADATX  = NUNDEF
FI%NCOGRIF = NUNDEF
FI%NSEC1   = NUNDEF
FI%NIDCEN  = FA%NIDCEN

IF (LLMLAM) THEN
  ALLOCATE (FI%NSC2ALF (MAX (-KTYPTR, KPXTRO)-1))
ELSE
  ALLOCATE (FI%NSC2ALF (KPXTRO-1))
ENDIF

END SUBROUTINE NEW_FICHIER

SUBROUTINE FREE_FICHIER (FI)

TYPE (FAFICH) :: FI
TYPE (FAFICH) :: FIDUM

IF (ASSOCIATED (FI%MADATE )) DEALLOCATE (FI%MADATE )
IF (ASSOCIATED (FI%MADATX )) DEALLOCATE (FI%MADATX )
IF (ASSOCIATED (FI%NSC2ALF)) DEALLOCATE (FI%NSC2ALF)
IF (ASSOCIATED (FI%FLAP1D )) DEALLOCATE (FI%FLAP1D )
IF (ASSOCIATED (FI%FLAP1DA)) DEALLOCATE (FI%FLAP1DA)

FI = FIDUM

END SUBROUTINE FREE_FICHIER

SUBROUTINE NEW_FA_DEFAULT ()
USE LFIMOD, ONLY : LFICOM_DEFAULT, NEW_LFI_DEFAULT
INTEGER :: IERR
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('FA_COM:NEW_FA_DEFAULT',0,ZHOOK_HANDLE)

CALL NEW_LFI_DEFAULT
IF (.NOT. FA_COM_DEFAULT_INIT) THEN
  CALL NEW_FA (FA_COM_DEFAULT, IERR)
  FA_COM_DEFAULT_INIT = .TRUE.
  FA_COM_DEFAULT%LFI => LFICOM_DEFAULT
ENDIF

IF (LHOOK) CALL DR_HOOK ('FA_COM:NEW_FA_DEFAULT',1,ZHOOK_HANDLE)

END SUBROUTINE NEW_FA_DEFAULT

SUBROUTINE NEW_FA (FA, KERR, KPXTRO, KPXLAT,  &
&                  KPXNIV, KPNXFA, KPNXCA)
TYPE(FA_COM) :: FA
INTEGER, INTENT(OUT) :: KERR
INTEGER, OPTIONAL, INTENT(IN) :: KPXTRO
INTEGER, OPTIONAL, INTENT(IN) :: KPXLAT
INTEGER, OPTIONAL, INTENT(IN) :: KPXNIV
INTEGER, OPTIONAL, INTENT(IN) :: KPNXFA
INTEGER, OPTIONAL, INTENT(IN) :: KPNXCA
REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('FA_COM:NEW_FA',0,ZHOOK_HANDLE)


FA%JPXNOM = JPXNOM
FA%JPXPRF = JPXPRF
FA%JPXSUF = JPXSUF

!
! Reglage de la troncature maximum gerable (JPXTRO)
! et du nombre maximum de niveaux verticaux (JPXNIV)
!
#if defined ( HIGHRES )
!
!     Setup high resolution parameters 
!
FA%JPXTRO = 4000
FA%JPXLAT = 4000
FA%JPXNIV = 200 
!
#else
! 
!     Setup low resolution parameters to save memory
!
FA%JPXTRO = 599 
FA%JPXLAT = 1200 
FA%JPXNIV = 200 
!
#endif
!
!
FA%JPNXFA=20 
FA%JPNXCA=20 

IF (PRESENT (KPXTRO)) FA%JPXTRO = INT (KPXTRO, JPLIKB)
IF (PRESENT (KPXLAT)) FA%JPXLAT = INT (KPXLAT, JPLIKB)
IF (PRESENT (KPXNIV)) FA%JPXNIV = INT (KPXNIV, JPLIKB)
IF (PRESENT (KPNXFA)) FA%JPNXFA = INT (KPNXFA, JPLIKB)
IF (PRESENT (KPNXCA)) FA%JPNXCA = INT (KPNXCA, JPLIKB)

FA%JPLDAT=11 
FA%JPUILA=3 
FA%JPTNIV=14
CALL CPARAMS (FA%JPXLAT, FA%JPXTRO, &
&   FA%JPXAU1, FA%JPXLON, FA%JPXAU2, FA%JPXPAH, &
&   FA%JPXIND, FA%JPXGEO, FA%JPXCSP, FA%JPXPDG, &
&   FA%JPXCHA)
FA%JPLB1P=19 
FA%JPLB2P=17 
FA%JPCADI=5 
FA%JPCAFS=4 
FA%JPNVER=1 
FA%CPCADI='CADRE-DIMENSIONS' 
FA%CPCAFS='CADRE-FRANKSCHMI'
FA%CPCARP='CADRE-REDPOINPOL' 
FA%CPCACH='CADRE-FOCOHYBRID'
FA%CPCASL='CADRE-SINLATITUD' 
FA%CPDATE='DATE-DES-DONNEES' 
FA%CPDATX='DATX-DES-DONNEES' 

FA%JPSEC1=37
FA%JPSEC2=22+MAX(FA%JPXTRO-1,FA%JPXLAT)
FA%JPSEC4=100
FA%JPXPAR=500

FA%XLAP1D  => NULL () 
FA%XLAP1DA => NULL ()
FA%XLAP2D  => NULL () 
FA%XLAP2DA => NULL ()

ALLOCATE (                                                       &
& FA%CTNPRF (FA%JPTNIV),        FA%NIVDSC (0:4,0:FA%JPTNIV),     &
& FA%CADRE  (FA%JPNXCA),        FA%FICHIER (0:FA%JPNXFA),        &
& FA%NULIND (FA%JPNXFA),        FA%NCAIND (FA%JPNXCA),           &
& STAT = KERR )
IF (KERR /= 0) GOTO 999

FA%FACADE_LLPREA = .TRUE.
FA%FACAGE_LLPREA = .TRUE.
FA%FACIES_LLPREA = .TRUE.
FA%FACTUM_LLPREA = .TRUE.
FA%FAGIOT_LLPREA = .TRUE.
FA%FALIMU_LLPREA = .TRUE.
FA%FAMISO_LLPREA = .TRUE.
FA%FANERG_LLPREA = .TRUE.
FA%FANMSG_LLPREA = .TRUE.
FA%FANUCA_LLPREA = .TRUE.
FA%FANUMU_LLPREA = .TRUE.
FA%FAREGI_LLPREA = .TRUE.
FA%FARFLU_LLPREA = .TRUE.
FA%FARINE_LLPREA = .TRUE.
FA%FAVORI_LLPREA = .TRUE.
FA%FAXION_LLPREA = .TRUE.
FA%FARINE_LLDEFM = .FALSE.
FA%NULOUT        = 0
FA%LOPENMP       = .TRUE.

999 CONTINUE

IF (LHOOK) CALL DR_HOOK ('FA_COM:NEW_FA',1,ZHOOK_HANDLE)

END SUBROUTINE NEW_FA

SUBROUTINE FREE_FA (FA, KERR)
TYPE(FA_COM) :: FA
INTEGER, INTENT(OUT) :: KERR

INTEGER (KIND=JPLIKB) :: ICAD, IFIC

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('FA_COM:FREE_FA',0,ZHOOK_HANDLE)

IF (ASSOCIATED (FA%XLAP1D  )) DEALLOCATE (FA%XLAP1D  )
IF (ASSOCIATED (FA%XLAP1DA )) DEALLOCATE (FA%XLAP1DA )
IF (ASSOCIATED (FA%XLAP2D  )) DEALLOCATE (FA%XLAP2D  )
IF (ASSOCIATED (FA%XLAP2DA )) DEALLOCATE (FA%XLAP2DA )

NULLIFY (FA%XLAP1D  )
NULLIFY (FA%XLAP1DA )
NULLIFY (FA%XLAP2D  )
NULLIFY (FA%XLAP2DA )

DEALLOCATE (               &
& FA%NIVDSC,  FA%CTNPRF,   &
& STAT = KERR )
IF (KERR /= 0) GOTO 999

NULLIFY (FA%NIVDSC,  FA%CTNPRF)


IF (ASSOCIATED (FA%YGR1TAB)) THEN
  DEALLOCATE (FA%YGR1TAB, STAT = KERR)
  IF (KERR /= 0) GOTO 999
  NULLIFY (FA%YGR1TAB)
ENDIF

DO ICAD = 1, INT (UBOUND (FA%CADRE, 1), JPLIKB)
  CALL FREE_CADRE (FA%CADRE (ICAD))
ENDDO

DO IFIC = 0, INT (UBOUND (FA%FICHIER, 1), JPLIKB)
  CALL FREE_FICHIER (FA%FICHIER (IFIC))
ENDDO

DEALLOCATE (FA%CADRE, FA%FICHIER, STAT = KERR)
NULLIFY (FA%CADRE, FA%FICHIER)

999 CONTINUE

IF (LHOOK) CALL DR_HOOK ('FA_COM:FREE_FA',1,ZHOOK_HANDLE)

END SUBROUTINE FREE_FA


SUBROUTINE CPARAMS (KPXLAT, KPXTRO, &
&   KPXAU1, KPXLON, KPXAU2, KPXPAH, &
&   KPXIND, KPXGEO, KPXCSP, KPXPDG, &
&   KPXCHA)

INTEGER (KIND=JPLIKB), INTENT (IN)  :: KPXLAT, KPXTRO
INTEGER (KIND=JPLIKB), INTENT (OUT) :: &
&      KPXAU1, KPXLON, KPXAU2, KPXPAH, &
&      KPXIND, KPXGEO, KPXCSP, KPXPDG, &
&      KPXCHA

KPXAU1=(1+KPXLAT)/2 
KPXLON=2*KPXLAT 
KPXAU2=(2*KPXTRO)+4 
KPXPAH=(8*(8/KPXAU1)+KPXAU1*(KPXAU1/8))    &
&          /((8/KPXAU1)+(KPXAU1/8)) 
KPXIND=(KPXAU1*(KPXAU1/KPXAU2)+KPXAU2*     &
&         (KPXAU2/KPXAU1))                 &
&         /((KPXAU1/KPXAU2)+(KPXAU2/KPXAU1)) 
KPXGEO=(12*(12/KPXAU1)+KPXAU1*(KPXAU1/12)) &
&          /((12/KPXAU1)+(KPXAU1/12)) 
!KPXCSP=(1+KPXTRO)*(2+KPXTRO) 
! that formulation is valid for triangular truncation (global model) only
! the next one, valid for square truncation (limited area model) returns a bigger value
KPXCSP=(2*KPXTRO+1)*(2*KPXTRO+1)
KPXPDG=KPXLON*KPXLAT 
KPXCHA=(KPXCSP*(KPXCSP/KPXPDG)+            &
&         KPXPDG*(KPXPDG/KPXCSP))          &
&         /((KPXCSP/KPXPDG)+(KPXPDG/KPXCSP)) 
END SUBROUTINE CPARAMS

END MODULE FA_MOD

