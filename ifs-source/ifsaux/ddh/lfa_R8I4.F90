!OPTIONS NODOUBLE

      SUBROUTINE LFAAFFC(KLBOUC,KUL,CDNA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaaffc* Affichage de caracteres.
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
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KLBOUC,KUL,ILONG,IERR,JLONG,ILDTA
      CHARACTER*(*) CDNA
      CHARACTER*400 CLDTA(KLBOUC)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAAFFC',0,ZHOOK_HANDLE)
      CALL LFALECC(KUL,CDNA,KLBOUC,CLDTA,ILONG,IERR)
      DO JLONG=1,ILONG
        CALL LONC(CLDTA(JLONG),ILDTA)
        WRITE(*,'(a)') CLDTA(JLONG)(1:ILDTA)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAAFFC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAAFFC
!
!
      SUBROUTINE LFAAFFI(KLBOUC,KUL,CDNA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaaffi* Affichage d'entiers.
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
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KLBOUC,KUL,ILONG,IERR,JLONG
      CHARACTER*(*) CDNA
      INTEGER(KIND=JPINTUSR) IDTA(KLBOUC)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAAFFI',0,ZHOOK_HANDLE)
      CALL LFALECI(KUL,CDNA,KLBOUC,IDTA,ILONG,IERR)
      DO JLONG=1,ILONG
        PRINT*,IDTA(JLONG)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAAFFI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAAFFI
!
!
      SUBROUTINE LFAAFFR(KLBOUC,KUL,CDNA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaaffr* Affichage de reels.
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
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KLBOUC,KUL,ILONG,IERR,JLONG
      CHARACTER*(*) CDNA
      REAL(KIND=JPREEUSR) ZREEL(KLBOUC)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAAFFR',0,ZHOOK_HANDLE)
      CALL LFALECR(KUL,CDNA,KLBOUC,ZREEL,ILONG,IERR)
      DO JLONG=1,ILONG
        PRINT*,ZREEL(JLONG)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAAFFR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAAFFR
!
!
      SUBROUTINE LFAAVAN(KUL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAAVAN* Saute l'article courant dans un fichier LFA.
      ! **** *LFAAVAN* Step over current article in an LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul               unite logique du fichier.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul               logical unit of the LFA file.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAAVAN',0,ZHOOK_HANDLE)
      IF(LGPOINT(KUL)) THEN
        !
        ! Le pointeur est avant des donnees.
        ! On les saute.
        !
        READ(KUL)
        !
        ! Position du pointeur.
        !
        LGPOINT(KUL)=.FALSE.
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAAVAN',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAAVAN
!
!
      SUBROUTINE LFACAS(KUL,CDNA,CDTYPE,KLONG,KERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFACAS* Renseignements sur un article de fichier LFA.
      ! **** *LFACAS* Get documentation about a LFA article.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul               unite logique du fichier.
      ! cdna: si cdna=' ' on recherche l'article suivant.
      ! .         cdna est alors en entree/sortie,
      ! .         et en sortie il vaudra le nom de l'article suivant
      ! .         (si cet article existe).
      ! .         kerr...retour de recherche: 0 si OK,
      ! .                1 si fin de fichier.
      ! .     si cdna<>' ' cdna est le nom de l'article cherche.
      ! .          Il est alors en entree seulement.
      ! .         kerr...retour de recherche: 0 si OK,
      ! .                1 si article inexistant.
      ! En sortie:
      ! cdtype            type d'article: 'R4', 'I8', 'C '.
      ! klong             nombre d'elements de cet article.
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul               file logical unit.
      ! cdna: if cdna=' ' on looks for nbext article.
      ! .         cdna is then in input/output
      ! .         and in output it will receive next article name
      ! .         (if this article exists).
      ! .         kerr...return from search: 0 if OK,
      ! .                1 if end of file.
      ! .     if cdna<>' ' cdna is the name from required article.
      ! .          It is then in input olny.
      ! .         kerr...return from search: 0 if OK,
      ! .                1 if non-existant article.
      ! Output:
      ! cdtype            article type: 'R4', 'I8', 'C '.
      ! klong             numbre of elements in this article.
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KLONG,KERR,IMES,ITYPE
      LOGICAL LLERF
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      CHARACTER*(*) CDNA
      CHARACTER*(*) CDTYPE
      !
      ! On rend temporairement silencieux et tolerant le logiciel.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFACAS',0,ZHOOK_HANDLE)
      IMES=NMES(KUL)
      NMES(KUL)=0
      LLERF=LGERF(KUL)
      LGERF(KUL)=.FALSE.
      !
      ! Recherche de l'article suivant.
      !
      CALL LFAIPOS(KUL,CDNA,KERR,ITYPE,KLONG)
      IF(KERR.EQ.0) THEN
        !
        ! On n'etait pas en fin de fichier.
        ! On lit le contenu de la documentation
        ! pour fournir
        ! le type de d'article (reel, entier, etc...)
        !
        CALL LFAITYPE(ITYPE,CDTYPE)
      ENDIF
      !
      ! On remet le niveau de messagerie et de tolerance a l'initial.
      !
      NMES(KUL)=IMES
      LGERF(KUL)=LLERF
      IF (LHOOK) CALL DR_HOOK('LFACAS',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFACAS
!
!
      SUBROUTINE LFACOP(KULE,CDNAE,CDNAS,KULS)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFACOP* Copie d'un article d'un fichier LFA a un autre.
      ! **** *LFACOP* Copy one article from a LFA file to another.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kule unite logique du fichier LFA d'entree.
      ! cdnae nom de l'article a lire.
      ! cdnas nom sous lequel l'article est recopie.
      ! kuls unite logique du fichier LFA de sortie.
      ! En sortie:
      ! Le fichier d'unite logique kuls est augmente d'un article.
      ! --------------------------------------------------------------
      ! Input:
      ! kule logical unit of input LFA file.
      ! cdnae article name to be read.
      ! cdnas article name to be written out.
      ! kuls logical unit of output LFA file.
      ! Output:
      ! The file which logical unit is kuls receives one more article.
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KULE,KULS,ILONG,IERR
      CHARACTER*(*) CDNAE, CDNAS
      CHARACTER*2 CLTYPE
      !
      ! Renseignements sur l'article cdnae.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFACOP',0,ZHOOK_HANDLE)
      CALL LFACAS(KULE,CDNAE,CLTYPE,ILONG,IERR)
      IF(IERR.EQ.0) THEN
        !
        ! L'article existe.
        !
        IF(CLTYPE(1:1).EQ.'R') THEN
          !
          ! Article de type reel.
          !
          CALL LFAICOPR(KULE,CDNAE,CDNAS,ILONG,KULS)
        ELSEIF(CLTYPE(1:1).EQ.'I') THEN
          !
          ! Article de type entier.
          !
          CALL LFAICOPI(KULE,CDNAE,CDNAS,ILONG,KULS)
        ELSEIF(CLTYPE(1:1).EQ.'C') THEN
          !
          ! Article de type caractere.
          !
          CALL LFAICOPC(KULE,CDNAE,CDNAS,ILONG,KULS)
        ELSE
          PRINT*,'LFACOP/ATTENTION: type de donnee inconnu!...'
          PRINT*,CLTYPE
        ENDIF
      ELSE
        PRINT*,'LFACOP/ATTENTION: article ',CDNAE,' inexistant!...'
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFACOP',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFACOP
!
!
      SUBROUTINE LFAECRC(KUL,CDNA,CDCAR,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAECRC* Ecriture de caracteres sur fichier LFA.
      ! **** *LFAECRC* Write character data on LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul              unite logique du fichier.
      ! cdna             nom de l'article a ecrire.
      ! cdcar(1,klong)   caracteres a ecrire.
      ! klong            longueur de l'article a ecrire.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul              logical unit of LFA file.
      ! cdna             name of article to write.
      ! cdcar(1,klong)   characters to write.
      ! klong            length of article to write.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KLONG,ILNOMF,ILNA,INBCARE &
      ,JLONG,ILLOC,ITYPE,IPROD
      CHARACTER*(*) CDNA
      CHARACTER*(*) CDCAR(KLONG)
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAECRC',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfaecrc: ecriture de l''article ',CDNA
      ENDIF
      IF(CGTYPO(KUL).EQ.'R') THEN
        IF(CLLANG().EQ.'FRA') THEN
          PRINT* &
          ,'LFAECRC/ERREUR: ecriture sur fichier ouvert en lecture!...'
          PRINT*,'Unite logique: ',KUL
          CALL LONC(CGNOMF(KUL),ILNOMF)
          PRINT*,'Fichier ',CGNOMF(KUL)(1:ILNOMF)
          CALL LONC(CDNA,ILNA)
          PRINT*,'Article ',CDNA(1:ILNA)
          STOP 1
        ELSE
          PRINT* &
          ,'LFAECRC/ERROR: write on file opened in read!...'
          PRINT*,'Logical unit: ',KUL
          CALL LONC(CGNOMF(KUL),ILNOMF)
          PRINT*,'File ',CGNOMF(KUL)(1:ILNOMF)
          CALL LONC(CDNA,ILNA)
          PRINT*,'Article ',CDNA(1:ILNA)
          STOP 1
        ENDIF
      ENDIF
      !
      ! On determine le nombre maximal de caracteres par element
      ! du tableau d'entree.
      !
      INBCARE=0
      DO JLONG=1,KLONG
        CALL LONC(CDCAR(JLONG),ILLOC)
        INBCARE=MAX(INBCARE,ILLOC)
      ENDDO
      !
      ! Ecriture de l'autodocumentation de l'article.
      !
      ITYPE=-INBCARE ! type de donnee (< 0 pour les donnees caracteres).
      CALL LFAIDOC(KUL,ITYPE,KLONG,CDNA)
      !
      ! Ecriture de l'article caractere
      ! sous forme d'entiers.
      !
      IPROD=INBCARE*KLONG ! nombre total de caracteres sur l'ensemble du tableau.
      CALL LFAIECRCLOC(KUL,CDCAR,INBCARE,KLONG,IPROD)
      IF (LHOOK) CALL DR_HOOK('LFAECRC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAECRC
!
!
      SUBROUTINE LFAECRI(KUL,CDNA,KENTIER,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAECRI* Ecriture d'entiers sur fichier LFA.
      ! **** *LFAECRI* Write integer data of LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul                  unite logique du fichier.
      ! cdna                  nom de l'article a ecrire.
      ! kentier(1,klong)      entiers a ecrire.
      ! klong            longueur de l'article a ecrire.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul                  logical unit of LFA file.
      ! cdna                 name of article to write.
      ! kentier(1,klong)     integers to write.
      ! klong                length of article to write.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KLONG,ILNOMF,ILNA,ITYPE,IPREC,JLONG
      CHARACTER*(*) CDNA
      INTEGER(KIND=JPINTUSR) KENTIER(KLONG)
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAECRI',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfaecri: ecriture de ',CDNA
      ENDIF
      IF(CGTYPO(KUL).EQ.'R') THEN
        IF(CLLANG().EQ.'FRA') THEN
          PRINT* &
          ,'LFAECRI/ERREUR: ecriture sur fichier ouvert en lecture!...'
          PRINT*,'Unite logique: ',KUL
          CALL LONC(CGNOMF(KUL),ILNOMF)
          PRINT*,'Fichier ',CGNOMF(KUL)(1:ILNOMF)
          CALL LONC(CDNA,ILNA)
          PRINT*,'Article ',CDNA(1:ILNA)
          STOP 1
        ELSE
          PRINT* &
          ,'LFAECRI/ERROR: write on file opened in read!...'
          PRINT*,'Logical unit: ',KUL
          CALL LONC(CGNOMF(KUL),ILNOMF)
          PRINT*,'File ',CGNOMF(KUL)(1:ILNOMF)
          CALL LONC(CDNA,ILNA)
          PRINT*,'Article ',CDNA(1:ILNA)
          STOP 1
        ENDIF
      ENDIF
      !
      ! Ecriture de l'autodocumentation de l'article.
      !
      IF(NPRECI(KUL).EQ.8) THEN
        ITYPE=4
      ELSEIF(NPRECI(KUL).EQ.4) THEN
        ITYPE=2
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAECRI/ERREUR: type non prevu!...'
          PRINT*,NPRECI(KUL)
          STOP 1
        ELSE
          PRINT*,'LFAECRI/ERROR: type unexpected!...'
          PRINT*,NPRECI(KUL)
          STOP 1
        ENDIF
      ENDIF
      CALL LFAIDOC(KUL,ITYPE,KLONG,CDNA)
      !
      ! Ecriture des entiers dans le LFA.
      !
      !
      IPREC=NPRECI(KUL)
      IF(IPREC.EQ.JPINTUSR) THEN
        !
        ! La precision des entiers a ecrire
        ! est celle des entiers passes en argument.
        ! Ecriture directe donc.
        !
        WRITE(KUL) (KENTIER(JLONG),JLONG=1,KLONG)
      ELSEIF(IPREC.EQ.4) THEN
        !
        ! La precision des entiers a ecrire est de 4 octets.
        !
        CALL LFAIECRI4(KUL,KLONG,KENTIER)
      ELSEIF(IPREC.EQ.8) THEN
        !
        ! La precision des entiers a ecrire est de 8 octets.
        !
        CALL LFAIECRI8(KUL,KLONG,KENTIER)
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAECRI/ERREUR: precision de sortie impossible: ',IPREC
          STOP 1
        ELSE
          PRINT*,'LFAECRI/ERROR: ouput precision: ',IPREC
          STOP 1
        ENDIF
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAECRI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAECRI
!
!
      SUBROUTINE LFAECRR(KUL,CDNA,PREEL,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAECRR* Ecriture de reels sur fichier LFA.
      ! **** *LFAECRR* Write real data on LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul              unite logique du fichier.
      ! cdna             nom de l'article a ecrire.
      ! preel(1,klong)   reels a ecrire.
      ! klong            longueur de l'article a ecrire.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul              logical unit of LFA file.
      ! cdna             name of article to write.
      ! preel(1,klong)   real data to write.
      ! klong            length of article to write.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KLONG,ILNOMF,ILNA,ITYPE,JLONG,IPREC
      CHARACTER*(*) CDNA
      REAL(KIND=JPREEUSR) PREEL(KLONG)
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAECRR',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfaecrr: ecriture de l''article ',CDNA
      ENDIF
      IF(CGTYPO(KUL).EQ.'R') THEN
        IF(CLLANG().EQ.'FRA') THEN
          PRINT* &
          ,'LFAECRR/ERREUR: ecriture sur fichier ouvert en lecture!...'
          PRINT*,'Unite logique: ',KUL
          CALL LONC(CGNOMF(KUL),ILNOMF)
          PRINT*,'Fichier ',CGNOMF(KUL)(1:ILNOMF)
          CALL LONC(CDNA,ILNA)
          PRINT*,'Article ',CDNA(1:ILNA)
          STOP 1
        ELSE
          PRINT* &
          ,'LFAECRR/ERROR: write on file opened in read!...'
          PRINT*,'Logical unit: ',KUL
          CALL LONC(CGNOMF(KUL),ILNOMF)
          PRINT*,'File ',CGNOMF(KUL)(1:ILNOMF)
          CALL LONC(CDNA,ILNA)
          PRINT*,'Article ',CDNA(1:ILNA)
          STOP 1
        ENDIF
      ENDIF
      !
      ! Ecriture de l'autodocumentation de l'article.
      !
      IF(NPRECR(KUL).EQ.8) THEN
        ITYPE=1
      ELSEIF(NPRECR(KUL).EQ.4) THEN
        ITYPE=3
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAECRR/ERREUR: type non prevu!...'
          PRINT*,NPRECR(KUL)
          STOP 1
        ELSE
          PRINT*,'LFAECRR/ERROR: type unexpected!...'
          PRINT*,NPRECR(KUL)
          STOP 1
        ENDIF
      ENDIF
      CALL LFAIDOC(KUL,ITYPE,KLONG,CDNA)
      !
      ! Ecriture de l'article reel.
      !
      IPREC=NPRECR(KUL)
      IF(IPREC.EQ.JPREEUSR) THEN
        !
        ! La precision des reels a ecrire
        ! est celle des reels passes en argument.
        ! Ecriture directe donc.
        !
        WRITE(KUL) (PREEL(JLONG),JLONG=1,KLONG)
      ELSEIF(IPREC.EQ.4) THEN
        !
        ! La precision des reels a ecrire est de 4 octets.
        !
        CALL LFAIECRR4(KUL,KLONG,PREEL)
      ELSEIF(IPREC.EQ.8) THEN
        !
        ! La precision des reels a ecrire est de 8 octets.
        !
        CALL LFAIECRR8(KUL,KLONG,PREEL)
      ELSEIF(IPREC.EQ.16) THEN
        !
        ! La precision des reels a ecrire est de 16 octets.
        !
        CALL LFAIECRR16(KUL,KLONG,PREEL)
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAECRR/ERREUR: precision de sortie impossible: ',IPREC
          STOP 1
        ELSE
          PRINT*,'LFAECRR/ERROR: output precision: ',IPREC
          STOP 1
        ENDIF
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAECRR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAECRR
!
!
      SUBROUTINE LFAERF(KUL,LDERF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAERF* Niveau d'erreur tolere par le logiciel LFA.
      ! **** *LFAERF* Choose error level of LFA software.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul               unite logique du fichier.
      ! lderf             .true. si toute erreur doit etre fatale,
      ! .false. si aucune ne doit l'etre.
      ! En sortie:
      ! lgerf             .true. si toute erreur est fatale,
      ! .false. si aucune ne l'est.
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul               logical unit of LFA file.
      ! lderf             .true. if any error has to be fatal.
      ! .false. si none has to be.
      ! Output:
      ! lgerf             .true. if any error has to be fatal.
      ! .false. si none has to be.
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL
      LOGICAL LDERF
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAERF',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfaerf: lgerf(',KUL,' mis a ',LDERF
      ENDIF
      LGERF(KUL)=LDERF
      IF (LHOOK) CALL DR_HOOK('LFAERF',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAERF
!
!
      SUBROUTINE LFAFER(KUL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAFER* Fermeture de fichier LFA.
      ! **** *LFAFER* Close a LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul        unite logique du fichier.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul        logical unit of LFA file.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAFER',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfafer: fermeture de l''unite logique ',KUL
      ENDIF
      CLOSE(KUL)
      IF (LHOOK) CALL DR_HOOK('LFAFER',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAFER
!
!
      SUBROUTINE LFAFORVLC(KULE,KULS,CDTYPE,CDNOMA,KNOMAL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaforvlc* Formatte vers article LFA de caracteres.
      ! **** *lfaforvlc* Formatted data to LFA character article.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-02, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KULE,KULS,KNOMAL,JLONG
      CHARACTER*400 CLVAL(KNOMAL)
      CHARACTER*(*) CDTYPE,CDNOMA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAFORVLC',0,ZHOOK_HANDLE)
      DO JLONG=1,KNOMAL
        READ(KULE,FMT='(a)') CLVAL(JLONG)
      ENDDO
      CALL LFAECRC(KULS,CDNOMA,CLVAL,KNOMAL)
      IF (LHOOK) CALL DR_HOOK('LFAFORVLC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAFORVLC
!
!
      SUBROUTINE LFAFORVLI(KULE,KULS,CDTYPE,CDNOMA,KNOMAL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaforvli* Formatte vers article LFA d'entiers.
      ! **** *lfaforvli* Formatted data to LFA integer article.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-02, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KULE,KULS,KNOMAL,JLONG
      INTEGER(KIND=JPINTUSR) IVAL(KNOMAL)
      CHARACTER*(*) CDTYPE,CDNOMA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAFORVLI',0,ZHOOK_HANDLE)
      DO JLONG=1,KNOMAL
        READ(KULE,FMT=*) IVAL(JLONG)
      ENDDO
      CALL LFAECRI(KULS,CDNOMA,IVAL,KNOMAL)
      IF (LHOOK) CALL DR_HOOK('LFAFORVLI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAFORVLI
!
!
      SUBROUTINE LFAFORVLR(KULE,KULS,CDTYPE,CDNOMA,KNOMAL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaforvlr* Formatte vers article LFA de reels.
      ! **** *lfaforvlr* Formatted data to LFA real article.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-02, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KULE,KULS,KNOMAL,JLONG
      REAL(KIND=JPREEUSR) ZVAL(KNOMAL)
      CHARACTER*(*) CDTYPE,CDNOMA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAFORVLR',0,ZHOOK_HANDLE)
      DO JLONG=1,KNOMAL
        READ(KULE,FMT=*) ZVAL(JLONG)
      ENDDO
      CALL LFAECRR(KULS,CDNOMA,ZVAL,KNOMAL)
      IF (LHOOK) CALL DR_HOOK('LFAFORVLR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAFORVLR
!
!
      SUBROUTINE LFAICOPC(KULE,CDNAE,CDNAS,KLONG,KULS)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaicopc* Copie d'un article caracteres d'un fichier LFA a un autre.
      ! **** *lfaicopc* Copy a character LFA article from one file to another.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kule unite logique du fichier LFA d'entree.
      ! cdnae nom de l'article a lire.
      ! cdnas nom sous lequel l'article est recopie.
      ! klong longeur de l'article a copier.
      ! kuls unite logique du fichier LFA de sortie.
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! kule logical unit of LFA input file.
      ! cdnae name of article to read.
      ! cdnas name of article to write.
      ! klong length of article to copy.
      ! kuls logical unit of LFA output file.
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KULE,KLONG,KULS,ILONG,IERR
      CHARACTER*(*) CDNAE, CDNAS
      CHARACTER*400 CLNAS
      CHARACTER*400 CLBOUC(KLONG)
      !
      ! Lecture de l'article d'entree.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAICOPC',0,ZHOOK_HANDLE)
      CALL LFALECC(KULE,CDNAE,KLONG,CLBOUC,ILONG,IERR)
      !
      ! Nom de l'article de sortie.
      !
      IF(CDNAS.EQ.' ') THEN
        CLNAS=CDNAE
      ELSE
        CLNAS=CDNAS
      ENDIF
      !
      ! Ecriture sur le fichier de sortie.
      !
      CALL LFAECRC(KULS,CLNAS,CLBOUC,ILONG)
      IF (LHOOK) CALL DR_HOOK('LFAICOPC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAICOPC
!
!
      SUBROUTINE LFAICOPI(KULE,CDNAE,CDNAS,KLONG,KULS)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaicopi* Copie d'un article entier d'un fichier LFA a un autre.
      ! **** *lfaicopi* Copy an integer LFA article from one file to another.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kule unite logique du fichier LFA d'entree.
      ! cdnae nom de l'article a lire.
      ! cdnas nom sous lequel l'article est recopie.
      ! klong longeur de l'article a copier.
      ! kuls unite logique du fichier LFA de sortie.
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! kule logical unit of LFA input file.
      ! cdnae name of article to read.
      ! cdnas name of article to write.
      ! klong length of article to copy.
      ! kuls logical unit of LFA output file.
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KULE,KLONG,KULS,ILONG,IERR,IPREC
      CHARACTER*(*) CDNAE, CDNAS
      CHARACTER*200 CLNAS
      INTEGER(KIND=JPINTUSR) IBOUC(KLONG)
      CHARACTER*3 CLLANG
      !
      ! cgdoc: documentation du type d'article LFP ou LFA lu (LIG, LFP, R8-, etc...).
      ! Si cgdoc est écrit par lfalecr ou lfaleci, il peut valoir "R4-", "R8-", "I4-", etc..
      ! Si cgdoc est écrit par lfplecr ou lfpleci, il peut valoir "R4-", "R8-", "I4-", mais aussi "LIG", "LFP", etc..
      !
      CHARACTER*3 CGDOC
      COMMON/LFADOC/CGDOC
      !
      ! Lecture de l'article d'entree.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAICOPI',0,ZHOOK_HANDLE)
      CALL LFALECI(KULE,CDNAE,KLONG,IBOUC,ILONG,IERR)
      !
      ! Nom de l'article de sortie.
      !
      IF(CDNAS.EQ.' ') THEN
        CLNAS=CDNAE
      ELSE
        CLNAS=CDNAS
      ENDIF
      !
      ! Choix de la precision de sortie, en fonction
      ! de celle lue en entree.
      !
      IF(CGDOC.EQ.'I4-') THEN
        IPREC=4
      ELSEIF(CGDOC.EQ.'I8-') THEN
        IPREC=8
      ELSEIF(CGDOC.EQ.'IG-') THEN
        IPREC=16
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAICOPI/ERREUR: precision inconnue!...'
          PRINT*,CGDOC
          STOP 1
        ELSE
          PRINT*,'LFAICOPI/ERROR: unknown precision!...'
          PRINT*,CGDOC
          STOP 1
        ENDIF
      ENDIF
      CALL LFAPRECI(KULS,IPREC)
      !
      ! Ecriture sur le fichier de sortie.
      !
      CALL LFAECRI(KULS,CLNAS,IBOUC,ILONG)
      IF (LHOOK) CALL DR_HOOK('LFAICOPI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAICOPI
!
!
      SUBROUTINE LFAICOPR(KULE,CDNAE,CDNAS,KLONG,KULS)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaicopr* Copie d'un article reel d'un fichier LFA a un autre.
      ! **** *lfaicopr* Copy a real LFA article from one file to another.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kule unite logique du fichier LFA d'entree.
      ! cdnae nom de l'article a lire.
      ! cdnas nom sous lequel l'article est recopie.
      ! klong longeur de l'article a copier.
      ! kuls unite logique du fichier LFA de sortie.
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! kule logical unit of LFA input file.
      ! cdnae name of article to read.
      ! cdnas name of article to write.
      ! klong length of article to copy.
      ! kuls logical unit of LFA output file.
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KULE,KLONG,KULS,ILONG,IERR,IPREC
      CHARACTER*(*) CDNAE, CDNAS
      CHARACTER*200 CLNAS
      REAL(KIND=JPREEUSR) ZBOUC(KLONG)
      CHARACTER*3 CLLANG
      !
      ! cgdoc: documentation du type d'article LFP ou LFA lu (LIG, LFP, R8-, etc...).
      ! Si cgdoc est écrit par lfalecr ou lfaleci, il peut valoir "R4-", "R8-", "I4-", etc..
      ! Si cgdoc est écrit par lfplecr ou lfpleci, il peut valoir "R4-", "R8-", "I4-", mais aussi "LIG", "LFP", etc..
      !
      CHARACTER*3 CGDOC
      COMMON/LFADOC/CGDOC
      !
      ! Lecture de l'article d'entree.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAICOPR',0,ZHOOK_HANDLE)
      CALL LFALECR(KULE,CDNAE,KLONG,ZBOUC,ILONG,IERR)
      !
      ! Nom de l'article de sortie.
      !
      IF(CDNAS.EQ.' ') THEN
        CLNAS=CDNAE
      ELSE
        CLNAS=CDNAS
      ENDIF
      !
      ! Choix de la precision de sortie, en fonction
      ! de celle lue en entree.
      !
      IF(CGDOC.EQ.'R4-') THEN
        IPREC=4
      ELSEIF(CGDOC.EQ.'R8-') THEN
        IPREC=8
      ELSEIF(CGDOC.EQ.'RG-') THEN
        IPREC=16
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAICOPR/ERREUR: precision inconnue!...'
          PRINT*,CGDOC
          STOP 1
        ELSE
          PRINT*,'LFAICOPR/ERROR: unknown precision!...'
          PRINT*,CGDOC
          STOP 1
        ENDIF
      ENDIF
      CALL LFAPRECR(KULS,IPREC)
      !
      ! Ecriture sur le fichier de sortie.
      !
      CALL LFAECRR(KULS,CLNAS,ZBOUC,ILONG)
      IF (LHOOK) CALL DR_HOOK('LFAICOPR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAICOPR
!
!
      SUBROUTINE LFAIDOC(KUL,KTYPE,KLONG,CDNA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *lfaidoc* Ecriture de l'autodocumentation de l'article: type d'article, longueur et nom.
      ! **** *lfaidoc* Write autodocumentation of a LFA article.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul: unite logique
      ! ktype: type d'article.
      ! klong: longueur de l'article.
      ! cdna: nom de l'article.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul: logical unit of LFA file.
      ! ktype: article type.
      ! klong: article length.
      ! cdna: article name.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JNA
      INTEGER(KIND=JPINTUSR) KTYPE,KLONG,ILNA
      INTEGER(KIND=JPINTESB) ITYPEESB,ILONGESB,ILNAESB,ICHAR
      CHARACTER*(*) CDNA
      !
      ! Conversion vers entiers 4 bits.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIDOC',0,ZHOOK_HANDLE)
      ITYPEESB=KTYPE
      ILONGESB=KLONG
      !
      ! Longueur de la chaine de caracteres.
      !
      CALL LONC(CDNA,ILNA)
      ILNAESB=ILNA
      WRITE(KUL) ITYPEESB,ILONGESB,ILNAESB &
      ,(ICHAR(CDNA(JNA:JNA)),JNA=1,ILNAESB)
      IF (LHOOK) CALL DR_HOOK('LFAIDOC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIDOC
!
!
      SUBROUTINE LFAIECRCLOC(KUL,CDCAR,KNBCARE,KLONG,KPROD)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiecrcloc* Ecriture de caracteres un LFA.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KNBCARE,KLONG,KPROD,JLONG
      CHARACTER*(*) CDCAR(KLONG)
      CHARACTER*(KNBCARE) CLCAR(KLONG)
      INTEGER(KIND=JPINTUSR) IASCII(KPROD)
      !
      ! -------------------------------------------------
      ! On porte le tableau d'entree sur un tableau
      ! dont les elements sont plus courts (de longueur knbcare).
      ! -------------------------------------------------
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIECRCLOC',0,ZHOOK_HANDLE)
      DO JLONG=1,KLONG
        CLCAR(JLONG)=CDCAR(JLONG)
      ENDDO
      !
      ! -------------------------------------------------
      ! Ecriture des caracteres sur le fichier LFA.
      ! -------------------------------------------------
      !
      WRITE(KUL) CLCAR
      IF (LHOOK) CALL DR_HOOK('LFAIECRCLOC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIECRCLOC
!
!
      SUBROUTINE LFAIECRI4(KUL,KLONG,KENTIER)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiecrr4* Ecriture d'entiers la precision imposee 4 octets.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      INTEGER(KIND=4) IENTIER(KLONG)
      INTEGER(KIND=JPINTUSR) KENTIER(KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIECRI4',0,ZHOOK_HANDLE)
      DO JLONG=1,KLONG
        IENTIER(JLONG)=KENTIER(JLONG)
      ENDDO
      !
      ! Ecriture du tableau affecte.
      !
      WRITE(KUL) (IENTIER(JLONG),JLONG=1,KLONG)
      IF (LHOOK) CALL DR_HOOK('LFAIECRI4',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIECRI4
!
!
      SUBROUTINE LFAIECRI8(KUL,KLONG,KENTIER)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiecrr8* Ecriture d'entiers la precision imposee 8 octets.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      INTEGER(KIND=8) IENTIER(KLONG)
      INTEGER(KIND=JPINTUSR) KENTIER(KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIECRI8',0,ZHOOK_HANDLE)
      DO JLONG=1,KLONG
        IENTIER(JLONG)=KENTIER(JLONG)
      ENDDO
      !
      ! Ecriture du tableau affecte.
      !
      WRITE(KUL) (IENTIER(JLONG),JLONG=1,KLONG)
      IF (LHOOK) CALL DR_HOOK('LFAIECRI8',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIECRI8
!
!
      SUBROUTINE LFAIECRR16(KUL,KLONG,PREEL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiecrr16* Ecriture de reels a la precision imposee 16 octets.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
#ifdef REALHUGE
      REAL(KIND=16) ZREEL(KLONG)
#else
      REAL(KIND=8) ZREEL(KLONG)
#endif
      REAL(KIND=JPREEUSR) PREEL(KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIECRR16',0,ZHOOK_HANDLE)
      DO JLONG=1,KLONG
        ZREEL(JLONG)=PREEL(JLONG)
      ENDDO
      !
      ! Ecriture du tableau affecte.
      !
      WRITE(KUL) (ZREEL(JLONG),JLONG=1,KLONG)
      IF (LHOOK) CALL DR_HOOK('LFAIECRR16',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIECRR16
!
!
      SUBROUTINE LFAIECRR4(KUL,KLONG,PREEL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiecrr4* Ecriture de reels a la precision imposee 4 octets.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      REAL(KIND=4) ZREEL(KLONG)
      REAL(KIND=JPREEUSR) PREEL(KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIECRR4',0,ZHOOK_HANDLE)
      DO JLONG=1,KLONG
        ZREEL(JLONG)=PREEL(JLONG)
      ENDDO
      !
      ! Ecriture du tableau affecte.
      !
      WRITE(KUL) (ZREEL(JLONG),JLONG=1,KLONG)
      IF (LHOOK) CALL DR_HOOK('LFAIECRR4',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIECRR4
!
!
      SUBROUTINE LFAIECRR8(KUL,KLONG,PREEL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiecrr8* Ecriture de reels a la precision imposee 8 octets.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      REAL(KIND=8) ZREEL(KLONG)
      REAL(KIND=JPREEUSR) PREEL(KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIECRR8',0,ZHOOK_HANDLE)
      DO JLONG=1,KLONG
        ZREEL(JLONG)=PREEL(JLONG)
      ENDDO
      !
      ! Ecriture du tableau affecte.
      !
      WRITE(KUL) (ZREEL(JLONG),JLONG=1,KLONG)
      IF (LHOOK) CALL DR_HOOK('LFAIECRR8',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIECRR8
!
!
      SUBROUTINE LFAILECCLOC(KUL,CDCAR,KDIMB,KNBCARE,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaileccloc* Lecture d'entiers sur un LFA, convertis en caracteres.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KLONG,JLONG,KDIMB,KNBCARE
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      CHARACTER*(*) CDCAR(KDIMB)
      CHARACTER*(KNBCARE) CLCAR(KLONG)
      CHARACTER*3 CLLANG
      !
      ! Test de la position du pointeur.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAILECCLOC',0,ZHOOK_HANDLE)
      IF(.NOT.LGPOINT(KUL)) THEN
        !
        ! Le pointeur est avant une autodocumentation.
        ! C'est qu'il y a eu un probleme en amont.
        !
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAILECCLOC/ERREUR: pointeur non positionne' &
          ,' avant des donnees!...'
          STOP 1
        ELSE
          PRINT*,'LFAILECCLOC/ERROR: pointer location not before' &
          ,' data!...'
          STOP 1
        ENDIF
      ENDIF
      !
      ! -------------------------------------------------
      ! Lecture des entiers sur le fichier LFA.
      ! -------------------------------------------------
      !
      READ(KUL) CLCAR
      DO JLONG=1,KLONG
        CDCAR(JLONG)=CLCAR(JLONG)
      ENDDO
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      IF (LHOOK) CALL DR_HOOK('LFAILECCLOC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAILECCLOC
!
!
      SUBROUTINE LFAILECCLOC8(KUL,CDCAR,KDIMB,KNBCARE,KLONG,KPROD)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaileccloc8* Lecture d'entiers sur un LFA, convertis en caracteres.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KDIMB,KNBCARE,KLONG,JLONG &
      ,KPROD,JPROD,JNA,IPOS
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      CHARACTER*(*) CDCAR(KDIMB)
      CHARACTER*(KNBCARE) CLCAR(KLONG)
      INTEGER(KIND=JPINTESB) IASCIIESB(KPROD)
      CHARACTER*3 CLLANG
      !
      ! Test de la position du pointeur.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAILECCLOC8',0,ZHOOK_HANDLE)
      IF(.NOT.LGPOINT(KUL)) THEN
        !
        ! Le pointeur est avant une autodocumentation.
        ! C'est qu'il y a eu un probleme en amont.
        !
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAILECCLOC8/ERREUR: pointeur non positionne' &
          ,' avant des donnees!...'
          STOP 1
        ELSE
          PRINT*,'LFAILECCLOC8/ERROR: pointer location not before' &
          ,' data!...'
          STOP 1
        ENDIF
      ENDIF
      !
      ! -------------------------------------------------
      ! Lecture des entiers sur le fichier LFA.
      ! -------------------------------------------------
      !
      READ(KUL) (IASCIIESB(JPROD),JPROD=1,KPROD)
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      !
      ! -------------------------------------------------
      ! Conversion entiers > caracteres.
      ! -------------------------------------------------
      !
      DO JLONG=1,KLONG
        !
        ! On initialise l'element a une chaine blanche.
        !
        CDCAR(JLONG)=' '
        DO JNA=1,KNBCARE
          !
          ! On convertit le code ASCII en un caractere.
          !
          IPOS=(JLONG-1)*KNBCARE+JNA
          CDCAR(JLONG)(JNA:JNA)=CHAR(IASCIIESB(IPOS))
        ENDDO
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAILECCLOC8',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAILECCLOC8
!
!
      SUBROUTINE LFAILECI4(KUL,KLONG,KENTIER)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaileci4* Lecture d'entiers a 4 octets sur fichier.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      INTEGER(KIND=4) IENTIER(KLONG)
      INTEGER(KIND=JPINTUSR) KENTIER(KLONG)
      !
      ! Ecriture du tableau affecte.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAILECI4',0,ZHOOK_HANDLE)
      READ(KUL) (IENTIER(JLONG),JLONG=1,KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      DO JLONG=1,KLONG
        KENTIER(JLONG)=IENTIER(JLONG)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAILECI4',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAILECI4
!
!
      SUBROUTINE LFAILECI8(KUL,KLONG,KENTIER)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaileci8* Lecture d'entiers a 8 octets sur fichier.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      INTEGER(KIND=8) IENTIER(KLONG)
      INTEGER(KIND=JPINTUSR) KENTIER(KLONG)
      !
      ! Ecriture du tableau affecte.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAILECI8',0,ZHOOK_HANDLE)
      READ(KUL) (IENTIER(JLONG),JLONG=1,KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      DO JLONG=1,KLONG
        KENTIER(JLONG)=IENTIER(JLONG)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAILECI8',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAILECI8
!
!
      SUBROUTINE LFAILECR16(KUL,KLONG,PREEL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfailecr16* Lecture de reels 16 octets sur fichier.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
#ifdef REALHUGE
      REAL(KIND=16) ZREEL(KLONG)
#else
      REAL(KIND=8) ZREEL(KLONG)
#endif
      REAL(KIND=JPREEUSR) PREEL(KLONG)
      !
      ! Ecriture du tableau affecte.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAILECR16',0,ZHOOK_HANDLE)
      READ(KUL) (ZREEL(JLONG),JLONG=1,KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      DO JLONG=1,KLONG
        PREEL(JLONG)=ZREEL(JLONG)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAILECR16',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAILECR16
!
!
      SUBROUTINE LFAILECR4(KUL,KLONG,PREEL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfailecr4* Lecture de reels 4 octets sur fichier.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      REAL(KIND=4) ZREEL(KLONG)
      REAL(KIND=JPREEUSR) PREEL(KLONG)
      !
      ! Ecriture du tableau affecte.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAILECR4',0,ZHOOK_HANDLE)
      READ(KUL) (ZREEL(JLONG),JLONG=1,KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      DO JLONG=1,KLONG
        PREEL(JLONG)=ZREEL(JLONG)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAILECR4',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAILECR4
!
!
      SUBROUTINE LFAILECR8(KUL,KLONG,PREEL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfailecr8* Lecture de reels 8 octets sur fichier.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,JLONG,KLONG
      REAL(KIND=8) ZREEL(KLONG)
      REAL(KIND=JPREEUSR) PREEL(KLONG)
      !
      ! Ecriture du tableau affecte.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAILECR8',0,ZHOOK_HANDLE)
      READ(KUL) (ZREEL(JLONG),JLONG=1,KLONG)
      !
      ! On affecte un tableau dans l'autre
      ! afin que le changement de precision s'opere.
      !
      DO JLONG=1,KLONG
        PREEL(JLONG)=ZREEL(JLONG)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAILECR8',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAILECR8
!
!
      SUBROUTINE LFAIMINMC(KUL,CDNA,CDTYPE,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiminmc* Extrema d'un article de caracteres de fichier LFA.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kul unite logique du fichier LFA d'entree.
      ! cdna nom de l'article a lire.
      ! cdtype type d'article (reel, entier, etc...).
      ! klong longueur de cet article.
      ! En sortie:
      ! Extrema sur output standard.
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KLONG,ILONG,IERR,JLONG &
      ,ILNOM,ILMIN,ILMAX,IGOL18
      CHARACTER*400 CLDONMIN,CLDONMAX
      CHARACTER*6 CLLONG
      CHARACTER*400 CLDON(KLONG)
      CHARACTER*80 CDNA
      CHARACTER*(*) CDTYPE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIMINMC',0,ZHOOK_HANDLE)
      CALL LFALECC(KUL,CDNA,KLONG,CLDON,ILONG,IERR)
      CLDONMIN=CLDON(1)
      CLDONMAX=CLDON(1)
      DO JLONG=1,ILONG
        IF(CLDON(JLONG).LT.CLDONMIN) CLDONMIN=CLDON(JLONG)
        IF(CLDON(JLONG).GT.CLDONMAX) CLDONMAX=CLDON(JLONG)
      ENDDO
      CALL LONC(CDNA,ILNOM)
      IGOL18=18
      ILNOM=MAX(IGOL18,ILNOM)
      ILMIN=12
      ILMAX=ILMIN
      WRITE(CLLONG,FMT='(i6)') ILONG
      WRITE(*,'(50a)') CDNA(1:ILNOM),'|',CDTYPE,'| l=' &
      ,CLLONG(1:6),', min=',CLDONMIN(1:ILMIN),' max=',CLDONMAX(1:ILMAX)
      IF (LHOOK) CALL DR_HOOK('LFAIMINMC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIMINMC
!
!
      SUBROUTINE LFAIMINMI(KUL,CDNA,CDTYPE,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiminmi* Extrema d'un article d'entiers de fichier LFA.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kul unite logique du fichier LFA d'entree.
      ! cdna nom de l'article a lire.
      ! cdtype type d'article (reel, entier, etc...).
      ! klong longueur de cet article.
      ! En sortie:
      ! Extrema sur output standard.
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      INTEGER(KIND=JPINTUSR) KUL,KLONG,ILONG,IERR,IMIN,IMAX &
      ,JLONG,ILNOM,ILMOY,ILRCM,ILMIN,ILMAX,IOPT,INC
      REAL(KIND=JPREEUSR) ZMOY,ZRCM,ZVALTMP
      CHARACTER*80 CLMIN,CLMAX,CLMOY,CLRCM,CLLONG
      CHARACTER*80 CDNA
      CHARACTER*(*) CDTYPE
      INTEGER(KIND=JPINTUSR) IDON(KLONG),IGOL18
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIMINMI',0,ZHOOK_HANDLE)
      CALL LFALECI(KUL,CDNA,KLONG,IDON,ILONG,IERR)
      IMIN=IDON(1)
      IMAX=IDON(1)
      ZMOY=0.
      ZRCM=0.
      DO JLONG=1,ILONG
        IF(IDON(JLONG).LT.IMIN) IMIN=IDON(JLONG)
        IF(IDON(JLONG).GT.IMAX) IMAX=IDON(JLONG)
        ZVALTMP=REAL(IDON(JLONG))
        ZMOY=ZMOY+ZVALTMP
        ZRCM=ZRCM+ZVALTMP*ZVALTMP
      ENDDO
      ZMOY=ZMOY/REAL(ILONG)
      ZRCM=SQRT(ZRCM/REAL(ILONG))
      CALL LONC(CDNA,ILNOM)
      IGOL18=18
      ILNOM=MAX(IGOL18,ILNOM)
      WRITE(CLMIN,FMT='(4x,i8)') IMIN
      WRITE(CLMAX,FMT='(4x,i8)') IMAX
      IOPT=-2
      INC=3
      CALL REECAR(ZMOY,IOPT,INC,CLMOY,ILMOY)
      CALL REECAR(ZRCM,IOPT,INC,CLRCM,ILRCM)
      ILMIN=12
      ILMAX=ILMIN
      WRITE(CLLONG,FMT='(i6)') ILONG
      IF(LGLANG) THEN
        WRITE(*,'(50a)') CDNA(1:ILNOM),'|',CDTYPE,'| l=' &
        ,CLLONG(1:6),', min=',CLMIN(1:ILMIN),' max=' &
        ,CLMAX(1:ILMAX),' moy= ',CLMOY(1:ILMOY),' rcm= ',CLRCM(1:ILRCM)
      ELSE
        WRITE(*,'(50a)') CDNA(1:ILNOM),'|',CDTYPE,'| l=' &
        ,CLLONG(1:6),', min=',CLMIN(1:ILMIN),' max=' &
        ,CLMAX(1:ILMAX),' mea= ',CLMOY(1:ILMOY),' rms= ',CLRCM(1:ILRCM)
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAIMINMI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIMINMI
!
!
      SUBROUTINE LFAIMINMR(KUL,CDNA,CDTYPE,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *lfaiminmr* Extrema d'un article de reels de fichier LFA.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kul unite logique du fichier LFA d'entree.
      ! cdna nom de l'article a lire.
      ! cdtype type d'article (reel, entier, etc...).
      ! klong longueur de cet article.
      ! En sortie:
      ! Extrema sur output standard.
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      INTEGER(KIND=JPINTUSR) KUL,KLONG,ILONG,IERR,JLONG,ILNOM &
      ,ILMOY,ILRCM,ILMIN,ILMAX,IGOL18,IOPT,INC
      REAL(KIND=JPREEUSR) ZMIN,ZMAX,ZMOY,ZRCM
      CHARACTER*80 CLMIN,CLMAX,CLMOY,CLRCM,CLLONG
      CHARACTER*80 CDNA
      CHARACTER*(*) CDTYPE
      REAL(KIND=JPREEUSR) ZDON(KLONG)
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIMINMR',0,ZHOOK_HANDLE)
      CALL LFALECR(KUL,CDNA,KLONG,ZDON,ILONG,IERR)
      ZMIN=ZDON(1)
      ZMAX=ZDON(1)
      ZMOY=0.
      ZRCM=0.
      DO JLONG=1,ILONG
        IF(ZDON(JLONG).LT.ZMIN) ZMIN=ZDON(JLONG)
        IF(ZDON(JLONG).GT.ZMAX) ZMAX=ZDON(JLONG)
        ZMOY=ZMOY+ZDON(JLONG)
        ZRCM=ZRCM+ZDON(JLONG)*ZDON(JLONG)
      ENDDO
      ZMOY=ZMOY/REAL(ILONG)
      ZRCM=SQRT(ZRCM/REAL(ILONG))
      CALL LONC(CDNA,ILNOM)
      IGOL18=18
      ILNOM=MAX(IGOL18,ILNOM)
      IOPT=-2
      INC=3
      CALL REECAR(ZMIN,IOPT,INC,CLMIN,ILMIN)
      CALL REECAR(ZMAX,IOPT,INC,CLMAX,ILMAX)
      CALL REECAR(ZMOY,IOPT,INC,CLMOY,ILMOY)
      CALL REECAR(ZRCM,IOPT,INC,CLRCM,ILRCM)
      WRITE(CLLONG,FMT='(i6)') ILONG
      IF(LGLANG) THEN
        WRITE(*,'(50a)') CDNA(1:ILNOM),'|',CDTYPE,'| l=' &
        ,CLLONG(1:6),', min= ',CLMIN(1:ILMIN),' max= ',CLMAX(1:ILMAX) &
        ,' moy= ',CLMOY(1:ILMOY),' rcm= ',CLRCM(1:ILRCM)
      ELSE
        WRITE(*,'(50a)') CDNA(1:ILNOM),'|',CDTYPE,'| l=' &
        ,CLLONG(1:6),', min= ',CLMIN(1:ILMIN),' max= ',CLMAX(1:ILMAX) &
        ,' mea= ',CLMOY(1:ILMOY),' rms= ',CLRCM(1:ILRCM)
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAIMINMR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIMINMR
!
!
      SUBROUTINE LFAINOMA(KASCII,KNOMA,KNA,CDNA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAINOMA* Passage suite d'entiers ASCII > chaine de caracteres.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kascii: tableau contenant la suite des entiers ASCII.
      ! knoma: dimension physique du tableau kascii.
      ! kna: nombre d'entiers reellement ecrits sur kascii.
      ! En sortie:
      ! cdna: nom de l'article.
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KNOMA,KNA,JNA
      INTEGER(KIND=JPINTUSR) KASCII(KNOMA)
      CHARACTER*(*) CDNA
      !
      ! -------------------------------------------------
      ! On cree le nom de l'article a partir
      ! du tableau d'entiers kascii.
      ! -------------------------------------------------
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAINOMA',0,ZHOOK_HANDLE)
      CDNA=' '
      DO JNA=1,KNA
        CDNA(JNA:JNA)=CHAR(KASCII(JNA))
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAINOMA',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAINOMA
!
!
      SUBROUTINE LFAIPOS(KUL,CDNA,KERR,KTYPE,KLONG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *lfaipos* Recherche de la position d'un article dans un fichier LFA.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! 98-03, Luc Gerard: protection de la boucle sur des indefs,
      ! intervenant en fin de fichier sur certains logiciels (1 HP!...).
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul          unite logique du fichier.
      ! cdna         nom de l'article.
      ! - si cdna chaine non blanche:
      ! On recherche l'article de nom cdna.
      ! - si cdna chaine blanche:
      ! On recherche l'article suivant.
      ! --------------------------------------------------------------------------
      ! En sortie:
      ! kerr indicateur d'erreur:
      ! +----------+--------------------------------------------------+
      ! | Valeur   |             Signification                        |
      ! +----------+--------------------------------------------------+
      ! | kerr=  0 | Tout est OK                                      |
      ! |          | Le pointeur sequentiel est alors positionne      |
      ! |          | sur l'article demande.                           |
      ! | kerr= -1 | Article inexistant                               |
      ! | kerr=-10 | Fin de fichier en recherche de l'article suivant |
      ! +----------+--------------------------------------------------+
      ! ktype type d'article (entier, reel, etc...)
      ! klong longueur de l'article.
      ! kerr         indicateur: 0 si article rencontre avant la fin du fichier.
      ! Le pointeur du fichier sequentiel est
      ! alors positionne sur l'article demande.
      ! ktype type d'article (entier, reel, etc...)
      ! klong longueur de l'article.
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      INTEGER(KIND=JPINTUSR) KUL,KERR,ILDNA,ILGNA &
      ,KTYPE,KLONG,ILNA,IVERSION,IASCII(JPNOMA),JNAUSR,INOMA
      INTEGER(KIND=JPINTESB) ITYPEESB,ILONGESB,ILNAESB &
      ,IASCIIESB(JPNOMA),IVERSIONESB,JNA
      LOGICAL LLREW
      CHARACTER*(*) CDNA
      CHARACTER*(JPNOMA) CLNA
      CHARACTER*3 CLLANG
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAIPOS',0,ZHOOK_HANDLE)
      CALL LONC(CDNA,ILDNA)
      !
      ! -------------------------------------------------
      ! Impressions preliminaires.
      ! -------------------------------------------------
      !
      IF(NMES(KUL).EQ.2) THEN
        PRINT*,'++ lfaipos: recherche de ',CDNA(1:ILDNA),'...'
      ENDIF
      !
      ! -------------------------------------------------
      ! Initialisations.
      ! -------------------------------------------------
      !
      KERR=-1
      LLREW=.FALSE. ! vrai si fichier rembobine, faux sinon.
      !
      ! -------------------------------------------------
      ! Test de la position du pointeur.
      ! -------------------------------------------------
      !
      IF(LGPOINT(KUL)) THEN
        !
        ! Le pointeur est avant un article de donnees.
        ! Ces donnees ne seraient-elles pas justement
        ! celles qu'il faut trouver?
        !
        ! Ce cas est peu probable dans le cadre
        ! d'une lecture des articles du fichier
        ! dans le desordre, mais est le cas general
        ! lors d'une lecture sequentielle.
        !
        CALL LONC(CGNA(KUL),ILGNA)
        IF(CGNA(KUL)(1:ILGNA).EQ.CDNA(1:ILDNA)) THEN
          !
          ! Le pointeur est positionne avant l'article
          ! de donnees que l'utilisateur
          ! souhaite lire. Il n'y a donc
          ! aucune action de lecture du fichier
          ! a effectuer.
          !
          KERR=0
          KTYPE=NTYPE(KUL)
          KLONG=NLONG(KUL)
          IF (LHOOK) CALL DR_HOOK('LFAIPOS',1,ZHOOK_HANDLE)
          RETURN
        ELSE
          !
          ! On est positionne sur un article de donnees,
          ! or il fauit aller chercher l'article suivant.
          ! Il faut donc sauter les donnees pour permettre
          ! de lire l'autodocumentation suivante.
          !
          READ(KUL)
          !
          ! Position du pointeur.
          !
          LGPOINT(KUL)=.FALSE.
        ENDIF
      ELSE
        !
        ! Le pointeur est avant un article d'autodocumentation.
        ! Il faut lire cette autodocumentation, ce qui est fait ci-apres!...
        !
      ENDIF
      !
      ! -------------------------------------------------
      ! Lecture de l'autodocumentation.
      ! -------------------------------------------------
      !
  100 CONTINUE
      READ(KUL,END=200) ITYPEESB,ILONGESB,ILNAESB &
  ,(IASCIIESB(JNA),JNA=1,MIN(ILNAESB,JPNOMA))
      !
      ! -------------------------------------------------
      ! Le min ci-dessus sert a proteger le logiciel
      ! du traitement errone de la fin de fichier par certains
      ! codes fortran (1 HP pour ne pas le nommer!...),
      ! lesquels bien qu'en situation de fin de fichier
      ! affectent des indefs a ilnaesb avant de sortir
      ! via l'etiquette "end=", et donc tournent sur une
      ! boucle jna=1,ilnaesb parfois gigantesque avant de sortir!...
      ! -------------------------------------------------
      !
      IF(ILNAESB.GT.JPNOMA) THEN
        !
        ! -------------------------------------------------
        ! Le present cas peut arriver si le code servant
        ! a lire un fichier a ete compile avec une valeur
        ! de jpnoma plus faible que celle du code
        ! ayant servit a l'ecrire, et que l'article
        ! ecrit avait une taille physique plus grande
        ! que le premier jpnoma.
        ! Ce cas est plus que rarissime: il n'est carrement encore
        ! jamais survenu!...
        ! -------------------------------------------------
        !
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAIPOS/ERREUR: ilnaesb > jpnoma!...'
          PRINT*,ILNAESB,JPNOMA
          STOP 1
        ELSE
          PRINT*,'LFAIPOS/ERROR: ilnaesb > jpnoma!...'
          PRINT*,ILNAESB,JPNOMA
          STOP 1
        ENDIF
      ENDIF
      KTYPE=ITYPEESB
      KLONG=ILONGESB
      ILNA=ILNAESB
      DO JNAUSR=1,ILNA
        IASCII(JNAUSR)=IASCIIESB(JNAUSR)
      ENDDO
      IF(NMES(KUL).EQ.2) PRINT*,'lfaipos: Article de type ',KTYPE &
      ,' de longueur ',KLONG
      !
      ! -------------------------------------------------
      ! On cree le nom clna de l'article a partir
      ! du tableau d'entiers iascii.
      ! -------------------------------------------------
      !
      INOMA=JPNOMA
      CALL LFAINOMA(IASCII,INOMA,ILNA,CLNA)
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.TRUE.
      NTYPE(KUL)=KTYPE
      NLONG(KUL)=KLONG
      CGNA(KUL)=CLNA
      !
      ! Sortie de renseignements.
      !
      IF(NMES(KUL).EQ.2) PRINT*,'lfaipos: Article teste: ',CLNA(1:ILNA)
      IF(CLNA(1:ILNA).EQ.CDNA(1:ILDNA)) THEN
        !
        ! -------------------------------------------------
        ! Le nom de l'article courant est bien celui
        ! demande.
        ! -------------------------------------------------
        !
        IF(NMES(KUL).EQ.2) PRINT*,'lfaipos: Article trouve.'
        KERR=0
        IF (LHOOK) CALL DR_HOOK('LFAIPOS',1,ZHOOK_HANDLE)
        RETURN
      ELSEIF(CDNA.EQ.' ') THEN
        !
        ! -------------------------------------------------
        ! On cherchait l'article suivant.
        ! -------------------------------------------------
        !
        IF(NMES(KUL).EQ.2) PRINT*,'lfaipos: Article suivant trouve.'
        CDNA=CLNA
        KERR=0
        IF (LHOOK) CALL DR_HOOK('LFAIPOS',1,ZHOOK_HANDLE)
        RETURN
      ELSE
        !
        ! -------------------------------------------------
        ! Le nom de l'article courant n'est pas celui
        ! demande. On saute les donnees associees.
        ! -------------------------------------------------
        !
        READ(KUL)
        !
        ! Position du pointeur.
        !
        LGPOINT(KUL)=.FALSE.
      ENDIF
      GOTO 100
  200 CONTINUE
      !
      ! On est en fin de fichier.
      !
      REWIND(KUL)
      READ(KUL) IVERSIONESB
      IVERSION=IVERSIONESB
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      IF(CDNA.EQ.' ') THEN
        !
        ! On recherchait l'article suivant
        ! mais on est en fait en fin de fichier.
        !
        KERR=-10
        IF (LHOOK) CALL DR_HOOK('LFAIPOS',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      IF(.NOT.LLREW) THEN
        !
        ! Si on est ici, c'est qu'a la fois on n'a pas rencontre
        ! l'article demande et que le fichier n'etait pas
        ! rembobine. Il est desormais rembobine, et on tente
        ! a nouveau la recherche de l'article demande.
        !
        LLREW=.TRUE.
        IF(NMES(KUL).GE.2) &
        PRINT*,'LFAIPOS/RECHERCHE de ',CDNA(1:ILDNA) &
        ,': fin du fichier lu et rebobinage...'
        GOTO 100
      ENDIF
      !
      ! L'article demande n'existe pas dans le fichier.
      ! Testons si le niveau d'erreur en fin de recherche
      ! est compatible avec la poursuite du programme.
      !
      IF(KERR.NE.0.AND.LGERF(KUL)) THEN
        !
        ! Le niveau d'erreur doit provoquer un ABORT.
        !
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAIPOS/ERREUR: article ',CDNA(1:ILDNA) &
          ,' inexistant!...'
          STOP 1
        ELSE
          PRINT*,'LFAIPOS/ERROR: article ',CDNA(1:ILDNA) &
          ,' not found!...'
          STOP 1
        ENDIF
      ELSEIF(KERR.NE.0.AND.NMES(KUL).GE.1) THEN
        !
        ! Le niveau d'erreur doit provoquer un message d'alerte.
        !
        IF(KERR.EQ.-1.AND.NMES(KUL).GE.1) THEN
          CALL LONC(CDNA,ILDNA)
          PRINT*,'LFAIPOS/ATTENTION: article ' &
          ,CDNA(1:ILDNA),' inexistant!...'
        ENDIF
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAIPOS',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAIPOS
!
!
      SUBROUTINE LFAITYPE(KTYPE,CDTYPE)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAITYPE* Type en clair de donnee LFA.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! ktype: entier definissant le type.
      ! En sortie:
      ! cdtype: chaine definissant le type (plus clair pour l'utilisateur).
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KTYPE
      CHARACTER*(*) CDTYPE
      CHARACTER*3 CLLANG
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAITYPE',0,ZHOOK_HANDLE)
      IF(KTYPE.EQ.5) THEN
        ! type=1: reel  128 bits (16 octets).
        CDTYPE='RG'
      ELSEIF(KTYPE.EQ.1) THEN
        ! type=1: reel   64 bits ( 8 octets).
        CDTYPE='R8'
      ELSEIF(KTYPE.EQ.3) THEN
        ! type=3: reel   32 bits ( 4 octets).
        CDTYPE='R4'
      ELSEIF(KTYPE.EQ.4) THEN
        ! type=4: entier 64 bits ( 8 octets).
        CDTYPE='I8'
      ELSEIF(KTYPE.EQ.2) THEN
        ! type=2: entier 32 bits ( 4 octets).
        CDTYPE='I4'
      ELSEIF(KTYPE.LT.0) THEN
        CDTYPE='C '
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAITYPE/ERREUR: type de donnee inconnu: ',KTYPE,'!...'
          STOP 1
        ELSE
          PRINT*,'LFAITYPE/ERROR: unknown data type: ',KTYPE,'!...'
          STOP 1
        ENDIF
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAITYPE',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAITYPE
!
!
      SUBROUTINE LFALAF(KUL,KULOUT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFALAF* Liste des articles d'un fichier LFA.
      ! **** *LFALAF* Article list of a LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul             unite logique du fichier.
      ! kulout          unite logique sur laquelle sortir la liste.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul             logical unit of LFA file.
      ! kulout          logical unit on which print out the list.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KULOUT,IART,ILONG,IERR,ILNA
      INTEGER(KIND=JPINTUSR) IVERSION
      INTEGER(KIND=JPINTESB) IVERSIONESB
      CHARACTER*80 CLNA,CLTYPE
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFALAF',0,ZHOOK_HANDLE)
      IF(LGLANG) THEN
        WRITE(KULOUT,'(3a)') 'LFALAF du fichier de nom ' &
        ,CGNOMF(KUL)(1:INDEX(CGNOMF(KUL),' ')),':'
      ELSE
        WRITE(KULOUT,'(3a)') 'LFALAF from file ' &
        ,CGNOMF(KUL)(1:INDEX(CGNOMF(KUL),' ')),':'
      ENDIF
      !
      ! -------------------------------------------------
      ! On rembobine le fichier.
      ! -------------------------------------------------
      !
      REWIND(KUL)
      READ(KUL) IVERSIONESB
      IVERSION=IVERSIONESB
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      !
      ! -------------------------------------------------
      ! On lit le fichier LFA sequentiellement jusqu'a la fin.
      ! -------------------------------------------------
      !
      IART=0
  100 CONTINUE
      CLNA=' '
      !
      ! -------------------------------------------------
      ! Avancee d'un article dans le fichier LFA.
      ! -------------------------------------------------
      !
      CALL LFACAS(KUL,CLNA,CLTYPE,ILONG,IERR)
      IF(IERR.EQ.0) THEN
        !
        ! On n'est pas en fin de fichier.
        ! On avance d'un article.
        !
        CALL LFAAVAN(KUL)
        IART=IART+1
        !
        ! Sortie des renseignements sur l'article.
        !
        CALL LONC(CLNA,ILNA)
        IF(LGLANG) THEN
          WRITE(KULOUT,'(3a,i9,2a)') &
          'Type |',CLTYPE(1:2) &
          ,'| Longueur ',ILONG &
          ,' | ',CLNA(1:ILNA)
        ELSE
          WRITE(KULOUT,'(3a,i9,2a)') &
          'Type |',CLTYPE(1:2) &
          ,'| Length   ',ILONG &
          ,' | ',CLNA(1:ILNA)
        ENDIF
        !
        ! Lecture de la suite du fichier.
        !
        GOTO 100
      ENDIF
      !
      ! -------------------------------------------------
      ! On remet le fichier au debut.
      ! -------------------------------------------------
      !
      REWIND(KUL)
      READ(KUL) IVERSIONESB
      IVERSION=IVERSIONESB
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      IF (LHOOK) CALL DR_HOOK('LFALAF',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFALAF
!
!
      SUBROUTINE LFALAFT(KUL,CDLIS,KDLIS,KNLIS)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFALAFT* Liste des articles d'un fichier LFA sur tableau de caracteres.
      ! **** *LFALAFT* Article list of a LFA file, on an array.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul             unite logique du fichier.
      ! kdlis           dimension physique du tableau cdlis.
      ! En sortie:
      ! knlis           nombre d'articles du fichier.
      ! Ce nombre est egalement le nombre d'elements ecrits sur cdlis.
      ! cdlis(1, ..., knlis) nom des articles du fichier.
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul            logical unit of LFA file.
      ! kdlis          physical dimension of array cdlis.
      ! Output:
      ! knlis          number of articles on the file.
      ! This number is also the number of elements written on cdlis.
      ! cdlis(1, ..., knlis) article names.
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KDLIS,KNLIS,ILONG,IERR
      INTEGER(KIND=JPINTUSR) IVERSION
      INTEGER(KIND=JPINTESB) IVERSIONESB
      CHARACTER*80 CLNA
      CHARACTER*2 CLTYPE
      CHARACTER*(*) CDLIS(KDLIS)
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFALAFT',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfalaft: entree.'
      ENDIF
      !
      ! -------------------------------------------------
      ! On rembobine le fichier.
      ! -------------------------------------------------
      !
      REWIND(KUL)
      READ(KUL) IVERSIONESB
      IVERSION=IVERSIONESB
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      !
      ! -------------------------------------------------
      ! On lit le fichier LFA sequentiellement jusqu'a la fin.
      ! -------------------------------------------------
      !
      KNLIS=0
  100 CONTINUE
      CLNA=' '
      !
      ! -------------------------------------------------
      ! Avancee d'un article dans le fichier LFA.
      ! -------------------------------------------------
      !
      CALL LFACAS(KUL,CLNA,CLTYPE,ILONG,IERR)
      IF(IERR.EQ.0) THEN
        !
        ! On n'est pas en fin de fichier.
        ! On avance d'un article.
        !
        CALL LFAAVAN(KUL)
        KNLIS=KNLIS+1
        IF(KNLIS.GT.KDLIS) THEN
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFALAFT/ERREUR: trop d''articles dans le fichier!...'
            PRINT*,'Recompiler!...'
            PRINT*,KNLIS,KDLIS
            STOP 1
          ELSE
            PRINT*,'LFALAFT/ERROR: too many articles in file!...'
            PRINT*,'Recompile!...'
            PRINT*,KNLIS,KDLIS
            STOP 1
          ENDIF
        ENDIF
        CDLIS(KNLIS)=CLNA
        !
        ! Lecture de la suite du fichier.
        !
        GOTO 100
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFALAFT',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFALAFT
!
!
      SUBROUTINE LFALECC(KUL,CDNA,KDIMB,CDCAR,KLONG,KERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFALECC* Lecture de caracteres sur fichier LFA.
      ! **** *LFALECC* Read character data on LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul              unite logique du fichier.
      ! cdna             nom de l'article.
      ! kdimb            dimension du tableau cdcar.
      ! En sortie:
      ! klong            nombre de chaines de caracteres lues.
      ! cdcar(1,klong)   chaines lues.
      ! kerr             indicateur d'erreur:
      ! +----------+-----------------------------------------------------+
      ! | Valeur   |             Signification                           |
      ! +----------+-----------------------------------------------------+
      ! | kerr=  0 | Tout est OK!                                        |
      ! | kerr= -1 | Article inexistant                                  |
      ! | kerr= -6 | Article plus long que le tableau devant le recevoir |
      ! | kerr= -8 | Mauvais type de donnees (reelles, entieres, car.)   |
      ! +----------+-----------------------------------------------------+
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul              logical unit of LFA file.
      ! cdna             article name.
      ! kdimb            physical dimension of array cdcar.
      ! Output:
      ! klong            number of character elements read.
      ! cdcar(1,klong)   character elements read.
      ! kerr             error indicator:
      ! +----------+-----------------------------------------------------+
      ! | Value    |             Meaning                                 |
      ! +----------+-----------------------------------------------------+
      ! | kerr=  0 | Everything is OK!                                   |
      ! | kerr= -1 | Article inexistant                                  |
      ! | kerr= -6 | Article bigger than array supposed to receive it    |
      ! | kerr= -8 | Wrong data type (real, integer, char.)              |
      ! +----------+-----------------------------------------------------+
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KDIMB,KLONG,KERR,ILNA,ITYPE &
      ,ITAILTAB,ITAILFIC,IPROD
      CHARACTER*(*) CDNA
      CHARACTER*(*) CDCAR(KDIMB)
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFALECC',0,ZHOOK_HANDLE)
      CALL LONC(CDNA,ILNA)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfalecc: lecture de ',CDNA
      ENDIF
      KLONG=0
      !
      ! Position du pointeur fichier sur l'article dans le LFA.
      !
      CALL LFAIPOS(KUL,CDNA,KERR,ITYPE,KLONG)
      IF(KERR.NE.0) THEN
        !
        ! Erreur retournee par lfaipos.
        ! On la transmet a l'appelant.
        !
        IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
        RETURN
      ELSEIF(ITYPE.GE.0) THEN
        !
        ! L'article cherche n'est pas de type caractere!...
        !
        IF(LGERF(KUL)) THEN
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFALECC/ERREUR: l''article ',CDNA(1:ILNA) &
            ,' n''est pas de type caractere!...'
            STOP 1
          ELSE
            PRINT*,'LFALECC/ERROR: the article ',CDNA(1:ILNA) &
            ,' is not character type!...'
            STOP 1
          ENDIF
        ELSE
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFALECC/ATTENTION: l''article ',CDNA(1:ILNA) &
            ,' n''est pas de type caractere!...'
            KERR=-8
            IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
            RETURN
          ELSE
            PRINT*,'LFALECC/WARNING: the article ',CDNA(1:ILNA) &
            ,' is not character type!...'
            KERR=-8
            IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
        ENDIF
      ELSE
        !
        ! L'article existe et a le bon type.
        ! On verifie que la dimension du tableau
        ! de caracteres permet de recevoir les donnees
        ! du fichier LFA.
        !
        IF(KLONG.GT.KDIMB) THEN
          !
          ! L'article est plus long que la dimension
          ! du tableau de caracteres.
          !
          IF(LGERF(KUL)) THEN
            IF(CLLANG().EQ.'FRA') THEN
              PRINT*,'LFALECC/ERREUR: article ',CDNA(1:ILNA) &
              ,' plus long que le tableau devant le recevoir!...'
              PRINT*,'En effet article de ',KLONG &
              ,' elements, or tableau de ',KDIMB,' elements.'
              STOP 1
            ELSE
              PRINT*,'LFALECC/ERROR: article ',CDNA(1:ILNA) &
              ,' longer than the array supposed to receive it!...'
              PRINT*,'Article with ',KLONG &
              ,' elements, array with ',KDIMB,' elements.'
              STOP 1
            ENDIF
          ELSE
            IF(CLLANG().EQ.'FRA') THEN
              PRINT*,'LFALECC/ATTENTION: article ',CDNA(1:ILNA) &
              ,' plus long que le tableau devant le recevoir!...'
              KERR=-6
              IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
              RETURN
            ELSE
              PRINT*,'LFALECC/WARNING: article ',CDNA(1:ILNA) &
              ,' longer than the array supposed to receive it!...'
              KERR=-6
              IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
              RETURN
            ENDIF
          ENDIF
        ENDIF
        !
        ! On verifie que la taille en caracteres
        ! de chaque element du tableau cdcar
        ! permet de recevoir les donnees du fichier LFA.
        !
        ITAILTAB=LEN(CDCAR(1)) ! taille de chaque element du tableau cdcar.
        ITAILFIC=-ITYPE ! taille de chaque element dans le fichier.
        IF(ITAILFIC.GT.ITAILTAB) THEN
          !
          ! L'article est plus long que la taille
          ! du tableau de caracteres.
          !
          IF(LGERF(KUL)) THEN
            IF(CLLANG().EQ.'FRA') THEN
              PRINT*,'LFALECC/ERREUR: l''article ',CDNA(1:ILNA) &
              ,' est constitue d''elements de taille' &
              ,' plus longue que celle du tableau devant les recevoir!...'
              PRINT*,'En effet, fichier: elements de taille ' &
              ,ITAILFIC,', alors que tableau: ',ITAILTAB
              STOP 1
            ELSE
              PRINT*,'LFALECC/ERROR: the article elements',CDNA(1:ILNA) &
              ,' are longer' &
              ,' than those from array supposed to receive them!...'
              PRINT*,'File: length of elements ' &
              ,ITAILFIC,', array: ',ITAILTAB
              STOP 1
            ENDIF
          ELSE
            IF(CLLANG().EQ.'FRA') THEN
              PRINT*,'LFALECC/ATTENTION: l''article ',CDNA(1:ILNA) &
              ,' est constitue d''elements de taille' &
              ,' plus longue que celle du tableau' &
              ,' devant les recevoir!...'
              KERR=-6
              IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
              RETURN
            ELSE
              PRINT*,'LFALECC/WARNING: the article elements ',CDNA(1:ILNA) &
              ,' are longer' &
              ,' than those from array supposed to receive them!...'
              KERR=-6
              IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
              RETURN
            ENDIF
          ENDIF
        ENDIF
        !
        ! Test de la position du pointeur.
        !
        IF(.NOT.LGPOINT(KUL)) THEN
          !
          ! Le pointeur est avant une autodocumentation.
          ! C'est qu'il y a eu un probleme en amont.
          !
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFAILECC/ERREUR: pointeur non positionne ' &
            ,'avant des donnees!...'
            STOP 1
          ELSE
            PRINT*,'LFAILECC/ERROR: pointer location not before ' &
            ,'data!...'
            STOP 1
          ENDIF
        ENDIF
        !
        ! Les tailles sont OK. On va lire sur le fichier LFA.
        !
        IF(NVERSION(KUL).EQ.8) THEN
          !
          ! Version du logiciel dans laquelle
          ! les caracteres etaient ecrits
          ! sous forme leur suite ASCII d'entiers.
          !
          IPROD=ITAILFIC*KLONG ! nombre total de caracteres sur l'ensemble du tableau.
          CALL LFAILECCLOC8(KUL,CDCAR,KDIMB,ITAILFIC,KLONG,IPROD)
        ELSE
          !
          ! Version du logiciel dans laquelle
          ! les caracteres sont ecrits
          ! sous forme d'un write implicite
          ! du tableau de caracteres.
          !
          CALL LFAILECCLOC(KUL,CDCAR,KDIMB,ITAILFIC,KLONG)
        ENDIF
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFALECC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFALECC
!
!
      SUBROUTINE LFALECI(KUL,CDNA,KDIMB,KENTIER,KLONG,KERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFALECI* Lecture d'entiers sur fichier LFA.
      ! **** *LFALECI* Read integer data on LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul              unite logique du fichier.
      ! cdna             nom de l'article.
      ! kdimb            dimension du tableau kentier.
      ! En sortie:
      ! klong            nombre d'entiers lus.
      ! kentier(1,klong) entiers lus.
      ! kerr             indicateur d'erreur:
      ! +----------+-----------------------------------------------------+
      ! | Valeur   |             Signification                           |
      ! +----------+-----------------------------------------------------+
      ! | kerr=  0 | Tout est OK!                                        |
      ! | kerr= -1 | Article inexistant                                  |
      ! | kerr= -6 | Article plus long que le tableau devant le recevoir |
      ! | kerr= -8 | Mauvais type de donnees (reelles, entieres, car.)   |
      ! +----------+-----------------------------------------------------+
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul              logical unit of LFA file.
      ! cdna             article name.
      ! kdimb            physical dimension of array kentier.
      ! Output:
      ! klong            number of integer elements read.
      ! kentier(1,klong) integer elements read.
      ! kerr             error indicator:
      ! +----------+-----------------------------------------------------+
      ! | Value    |             Meaning                                 |
      ! +----------+-----------------------------------------------------+
      ! | kerr=  0 | Everything is OK!                                   |
      ! | kerr= -1 | Article inexistant                                  |
      ! | kerr= -6 | Article bigger than array supposed to receive it    |
      ! | kerr= -8 | Wrong data type (real, integer, char.)              |
      ! +----------+-----------------------------------------------------+
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KDIMB,KLONG,KERR,ILNA,ITYPE &
      ,JLONG,IPREC
      CHARACTER*2 CLTYPE
      CHARACTER*(*) CDNA
      INTEGER(KIND=JPINTUSR) KENTIER(KDIMB)
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      !
      ! cgdoc: documentation du type d'article LFP ou LFA lu (LIG, LFP, R8-, etc...).
      ! Si cgdoc est écrit par lfalecr ou lfaleci, il peut valoir "R4-", "R8-", "I4-", etc..
      ! Si cgdoc est écrit par lfplecr ou lfpleci, il peut valoir "R4-", "R8-", "I4-", mais aussi "LIG", "LFP", etc..
      !
      CHARACTER*3 CGDOC
      COMMON/LFADOC/CGDOC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFALECI',0,ZHOOK_HANDLE)
      CALL LONC(CDNA,ILNA)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfaleci: lecture de ',CDNA
      ENDIF
      KLONG=0
      !
      ! Position du pointeur fichier sur l'article dans le LFA.
      !
      CALL LFAIPOS(KUL,CDNA,KERR,ITYPE,KLONG)
      CALL LFAITYPE(ITYPE,CLTYPE)
      CGDOC=CLTYPE//'-'
      IF(KERR.NE.0) THEN
        !
        ! Erreur retournee par lfaipos.
        ! On la transmet a l'appelant.
        !
        IF (LHOOK) CALL DR_HOOK('LFALECI',1,ZHOOK_HANDLE)
        RETURN
      ELSEIF(CLTYPE(1:1).NE.'I') THEN
        !
        ! L'article cherche n'est pas de type entier!...
        !
        IF(LGERF(KUL)) THEN
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFALECI/ERREUR: l''article ',CDNA(1:ILNA) &
            ,' n''est pas de type entier!...'
            STOP 1
          ELSE
            PRINT*,'LFALECI/ERROR: the article ',CDNA(1:ILNA) &
            ,' is not integer type!...'
            STOP 1
          ENDIF
        ELSE
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFALECI/ATTENTION: l''article ',CDNA(1:ILNA) &
            ,' n''est pas de type entier!...'
            KERR=-8
            IF (LHOOK) CALL DR_HOOK('LFALECI',1,ZHOOK_HANDLE)
            RETURN
          ELSE
            PRINT*,'LFALECI/WARNING: the article ',CDNA(1:ILNA) &
            ,' is not integer type!...'
            KERR=-8
            IF (LHOOK) CALL DR_HOOK('LFALECI',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
        ENDIF
      ELSE
        !
        ! L'article existe et a le bon type.
        ! On verifie que la dimension du tableau
        ! permet de recevoir les donnees
        ! du fichier LFA.
        !
        IF(KLONG.GT.KDIMB) THEN
          !
          ! L'article est plus long que la dimension
          ! du tableau de caracteres.
          !
          IF(LGERF(KUL)) THEN
            IF(CLLANG().EQ.'FRA') THEN
              PRINT*,'LFALECI/ERREUR: article ',CDNA(1:ILNA) &
              ,' plus long que le tableau devant le recevoir!...'
              PRINT*,'En effet article de ',KLONG &
              ,' elements, or tableau de ' &
              ,KDIMB,' elements.'
              STOP 1
            ELSE
              PRINT*,'LFALECI/ERROR: article ',CDNA(1:ILNA) &
              ,' longer than array supposed to receive it!...'
              PRINT*,'Article ',KLONG &
              ,' elements, array ' &
              ,KDIMB,' elements.'
              STOP 1
            ENDIF
          ELSE
            PRINT*,'LFALECI/ATTENTION: article ',CDNA(1:ILNA) &
            ,' plus long que le tableau devant le recevoir!...'
            KERR=-6
            IF (LHOOK) CALL DR_HOOK('LFALECI',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
        ENDIF
        !
        ! Test de la position du pointeur.
        !
        IF(.NOT.LGPOINT(KUL)) THEN
          !
          ! Le pointeur est avant une autodocumentation.
          ! C'est qu'il y a eu un probleme en amont.
          !
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFAILECI/ERREUR: pointeur non positionne' &
            ,' avant des donnees!...'
            STOP 1
          ELSE
            PRINT*,'LFAILECI/ERROR: pointer location not before' &
            ,' data!...'
            STOP 1
          ENDIF
        ENDIF
        !
        ! La taille est OK. Lecture sur le LFA.
        !
        !
        ! iprec: taille en octets des entiers du fichier.
        !
        IF(ITYPE.EQ.4) THEN
          IPREC=8
        ELSEIF(ITYPE.EQ.2) THEN
          IPREC=4
        ELSE
          PRINT*,'LFALECI/ERREUR: type de donnee non attendu!...'
          PRINT*,ITYPE
          STOP 1
        ENDIF
        IF(IPREC.EQ.JPINTUSR) THEN
          !
          ! La precision des entiers a lire
          ! est celle des entiers passes en argument.
          ! Lecture directe donc.
          !
          READ(KUL) (KENTIER(JLONG),JLONG=1,KLONG)
        ELSEIF(IPREC.EQ.4) THEN
          !
          ! La precision des entiers a lire est de 4 octets.
          !
          CALL LFAILECI4(KUL,KLONG,KENTIER)
        ELSEIF(IPREC.EQ.8) THEN
          !
          ! La precision des entiers a lire est de 8 octets.
          !
          CALL LFAILECI8(KUL,KLONG,KENTIER)
        ELSE
          PRINT*,'LFALECI/ERREUR: precision d''entree impossible: ',IPREC
          STOP 1
        ENDIF
        !
        ! Position du pointeur.
        !
        LGPOINT(KUL)=.FALSE.
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFALECI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFALECI
!
!
      SUBROUTINE LFALECR(KUL,CDNA,KDIMB,PREEL,KLONG,KERR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFALECR* Lecture de reels sur fichier LFA.
      ! **** *LFALECR* Read real data on LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul              unite logique du fichier.
      ! cdna             nom de l'article.
      ! kdimb            dimension du tableau preel.
      ! En sortie:
      ! klong            nombre de reels lus.
      ! preel(1,klong)   reels lus.
      ! kerr             indicateur d'erreur:
      ! +----------+-----------------------------------------------------+
      ! | Valeur   |             Signification                           |
      ! +----------+-----------------------------------------------------+
      ! | kerr=  0 | Tout est OK!                                        |
      ! | kerr= -1 | Article inexistant                                  |
      ! | kerr= -6 | Article plus long que le tableau devant le recevoir |
      ! | kerr= -8 | Mauvais type de donnees (reelles, entieres, car.)   |
      ! +----------+-----------------------------------------------------+
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul              logical unit of LFA file.
      ! cdna             article name.
      ! kdimb            physical dimension of array preel.
      ! Output:
      ! klong            number of real elements read.
      ! preel(1,klong)   real elements read.
      ! kerr             error indicator:
      ! +----------+-----------------------------------------------------+
      ! | Value    |             Meaning                                 |
      ! +----------+-----------------------------------------------------+
      ! | kerr=  0 | Everything is OK!                                   |
      ! | kerr= -1 | Article inexistant                                  |
      ! | kerr= -6 | Article bigger than array supposed to receive it    |
      ! | kerr= -8 | Wrong data type (real, integer, char.)              |
      ! +----------+-----------------------------------------------------+
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KDIMB,KLONG,KERR,ILNA,ITYPE &
      ,JLONG,IPREC
      CHARACTER*2 CLTYPE
      CHARACTER*(*) CDNA
      REAL(KIND=JPREEUSR) PREEL(KDIMB)
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      !
      ! cgdoc: documentation du type d'article LFP ou LFA lu (LIG, LFP, R8-, etc...).
      ! Si cgdoc est écrit par lfalecr ou lfaleci, il peut valoir "R4-", "R8-", "I4-", etc..
      ! Si cgdoc est écrit par lfplecr ou lfpleci, il peut valoir "R4-", "R8-", "I4-", mais aussi "LIG", "LFP", etc..
      !
      CHARACTER*3 CGDOC
      COMMON/LFADOC/CGDOC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFALECR',0,ZHOOK_HANDLE)
      CALL LONC(CDNA,ILNA)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfalecr: lecture de ',CDNA
      ENDIF
      KLONG=0
      !
      ! Position du pointeur fichier sur l'article dans le LFA.
      !
      CALL LFAIPOS(KUL,CDNA,KERR,ITYPE,KLONG)
      CALL LFAITYPE(ITYPE,CLTYPE)
      CGDOC=CLTYPE//'-'
      IF(KERR.NE.0) THEN
        !
        ! Erreur retournee par lfaipos.
        ! On la transmet a l'appelant.
        !
        IF (LHOOK) CALL DR_HOOK('LFALECR',1,ZHOOK_HANDLE)
        RETURN
      ELSEIF(CLTYPE(1:1).NE.'R') THEN
        !
        ! L'article cherche n'est pas de type reel!...
        !
        IF(LGERF(KUL)) THEN
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFALECR/ERREUR: l''article ',CDNA(1:ILNA) &
            ,' n''est pas de type reel!...'
            STOP 1
          ELSE
            PRINT*,'LFALECR/ERROR: the article ',CDNA(1:ILNA) &
            ,' is not real type!...'
            STOP 1
          ENDIF
        ELSE
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFALECR/ATTENTION: l''article ',CDNA(1:ILNA) &
            ,' n''est pas de type reel!...'
            KERR=-8
            IF (LHOOK) CALL DR_HOOK('LFALECR',1,ZHOOK_HANDLE)
            RETURN
          ELSE
            PRINT*,'LFALECR/WARNING: the article ',CDNA(1:ILNA) &
            ,' is not real type!...'
            KERR=-8
            IF (LHOOK) CALL DR_HOOK('LFALECR',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
        ENDIF
      ELSE
        !
        ! L'article existe et a le bon type.
        ! On verifie que la dimension du tableau
        ! permet de recevoir les donnees
        ! du fichier LFA.
        !
        IF(KLONG.GT.KDIMB) THEN
          !
          ! L'article est plus long que la dimension
          ! du tableau de caracteres.
          !
          IF(LGERF(KUL)) THEN
            IF(CLLANG().EQ.'FRA') THEN
              PRINT*,'LFALECR/ERREUR: article ',CDNA(1:ILNA) &
              ,' plus long que le tableau devant le recevoir!...'
              PRINT*,'En effet article de ',KLONG &
              ,' elements, or tableau de ' &
              ,KDIMB,' elements.'
              STOP 1
            ELSE
              PRINT*,'LFALECR/ERROR: article ',CDNA(1:ILNA) &
              ,' longer than array supposed to receive it!...'
              PRINT*,'Article ',KLONG &
              ,' elements, array ' &
              ,KDIMB,' elements.'
              STOP 1
            ENDIF
          ELSE
            PRINT*,'LFALECR/ATTENTION: article ' &
            ,CDNA(1:ILNA),' plus long que le tableau' &
            ,' devant le recevoir!...'
            KERR=-6
            IF (LHOOK) CALL DR_HOOK('LFALECR',1,ZHOOK_HANDLE)
            RETURN
          ENDIF
        ENDIF
        !
        ! Test de la position du pointeur.
        !
        IF(.NOT.LGPOINT(KUL)) THEN
          !
          ! Le pointeur est avant une autodocumentation.
          ! C'est qu'il y a eu un probleme en amont.
          !
          PRINT*,'LFAILECR/ERREUR: pointeur non positionne ' &
          ,'avant des donnees!...'
          STOP 1
        ENDIF
        !
        ! La taille est OK. Lecture sur le LFA.
        !
        !
        ! iprec: taille en octets des reels du fichier.
        !
        IF(ITYPE.EQ.5) THEN
          IPREC=16
        ELSEIF(ITYPE.EQ.1) THEN
          IPREC=8
        ELSEIF(ITYPE.EQ.3) THEN
          IPREC=4
        ELSE
          PRINT*,'LFALECR/ERREUR: type de donnee non attendu!...'
          PRINT*,ITYPE
          STOP 1
        ENDIF
        IF(IPREC.EQ.JPREEUSR) THEN
          !
          ! La precision des reels a lire
          ! est celle des reels passes en argument.
          ! Lecture directe donc.
          !
          READ(KUL) (PREEL(JLONG),JLONG=1,KLONG)
        ELSEIF(IPREC.EQ.4) THEN
          !
          ! La precision des reels a lire est de 4 octets.
          !
          CALL LFAILECR4(KUL,KLONG,PREEL)
        ELSEIF(IPREC.EQ.8) THEN
          !
          ! La precision des reels a lire est de 8 octets.
          !
          CALL LFAILECR8(KUL,KLONG,PREEL)
        ELSEIF(IPREC.EQ.16) THEN
          !
          ! La precision des reels a lire est de 16 octets.
          !
          CALL LFAILECR16(KUL,KLONG,PREEL)
        ELSE
          PRINT*,'LFALECR/ERREUR: precision d''entree impossible: ',IPREC
          STOP 1
        ENDIF
        !
        ! Position du pointeur.
        !
        LGPOINT(KUL)=.FALSE.
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFALECR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFALECR
!
!
      SUBROUTINE LFAMES(KUL,KMES)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAMES* Niveau de messagerie du logiciel LFA.
      ! **** *LFAMES* Print out level of LFA software.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul         unite logique du fichier.
      ! kmes        niveau de messagerie:
      ! si 0 aucun message sorti par le logiciel LFA.
      ! si 1 messages d'ATTENTION et d'ERREUR sorties.
      ! si 2 LFA est bavard (a reserver au debug de LFA...).
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul         logical unit of LFA file.
      ! kmes        print out level:
      ! if 0 no message print out.
      ! if 1 WARNING or ERROR messages print out.
      ! if 2 many comments print out (LFA debug mode only...).
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,KMES
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAMES',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfames: niveau de messagerie de l''unite logique ' &
        ,KUL,' porte a ',KMES
      ENDIF
      NMES(KUL)=KMES
      IF (LHOOK) CALL DR_HOOK('LFAMES',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAMES
!
!
      SUBROUTINE LFAMINM(KUL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAMINM* Extrema de tous les articles d'un fichier LFA.
      ! **** *LFAMINM* Extrema of all articles of a given LFA file.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kul unite logique du fichier LFA d'entree.
      ! En sortie:
      ! Extrema sur output standard.
      ! --------------------------------------------------------------
      ! Input:
      ! kul logical unit of LFA file.
      ! Output:
      ! Extrema on standard output.
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,ILONG,IERR
      INTEGER(KIND=JPINTUSR) IVERSION
      INTEGER(KIND=JPINTESB) IVERSIONESB
      CHARACTER*80 CLNA
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      CHARACTER*2 CLTYPE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAMINM',0,ZHOOK_HANDLE)
      IF(LGLANG) THEN
        WRITE(*,'(3a)') 'LFAMINM du fichier ' &
        ,CGNOMF(KUL)(1:INDEX(CGNOMF(KUL),' ')),':'
      ELSE
        WRITE(*,'(3a)') 'LFAMINM from file ' &
        ,CGNOMF(KUL)(1:INDEX(CGNOMF(KUL),' ')),':'
      ENDIF
      REWIND(KUL)
      READ(KUL) IVERSIONESB
      IVERSION=IVERSIONESB
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
  100 CONTINUE
      !
      ! -------------------------------------------------
      ! Avancee d'un article dans le fichier LFA.
      ! -------------------------------------------------
      !
      CLNA=' '
      CALL LFACAS(KUL,CLNA,CLTYPE,ILONG,IERR)
      IF(IERR.EQ.0) THEN
        !
        ! On n'est pas en fin de fichier.
        !
        IF(CLTYPE(1:1).EQ.'R') THEN
          !
          ! Article de type reel.
          !
          CALL LFAIMINMR(KUL,CLNA,CLTYPE,ILONG)
        ELSEIF(CLTYPE(1:1).EQ.'I') THEN
          !
          ! Article de type entier.
          !
          CALL LFAIMINMI(KUL,CLNA,CLTYPE,ILONG)
        ELSEIF(CLTYPE(1:1).EQ.'C') THEN
          !
          ! Article de type caractere.
          !
          CALL LFAIMINMC(KUL,CLNA,CLTYPE,ILONG)
        ELSE
          PRINT*,'LFAMINM/ATTENTION: type de donnee inconnu!...'
          PRINT*,CLTYPE
        ENDIF
        GOTO 100
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAMINM',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAMINM
!
!
      SUBROUTINE LFAOUV(KUL,CDNOMF,CDTYPO)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAOUV* Ouverture de fichier LFA.
      ! **** *LFAOUV* Open a LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul         unite logique du fichier.
      ! cdnomf      nom du fichier.
      ! cdtypo      type d'ouverture: 'R' READ, 'W' WRITE, 'A' APPEND, 'S' SCRATCH.
      ! En sortie:
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul         logical unit of LFA file.
      ! cdnomf      file name.
      ! cdtypo      opening type: 'R' READ, 'W' WRITE, 'A' APPEND, 'S' SCRATCH.
      ! Output:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL,ILNOMF,ILTYPO
      INTEGER(KIND=JPINTUSR) IVERSION
      INTEGER(KIND=JPINTESB) IVERSIONESB
      LOGICAL LLEX,LLLFA
      CHARACTER*(*) CDNOMF
      CHARACTER*(*) CDTYPO
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAOUV',0,ZHOOK_HANDLE)
      CALL LONC(CDNOMF,ILNOMF)
      CALL LONC(CDTYPO,ILTYPO)
      IVERSION=9
      !
      ! Type d'ouverture du fichier.
      !
      CGTYPO(KUL)=CDTYPO
      !
      ! Niveau de messagerie
      ! (valeur par defaut).
      !
      NMES(KUL)=1
      !
      ! Niveau d'erreur fatale par defaut.
      !
      CALL LFAERF(KUL,.TRUE.)
      !
      ! Sauvegarde du nom du fichier.
      !
      CGNOMF(KUL)=CDNOMF
      !
      ! Precision d'ecriture par defaut.
      !
      NPRECR(KUL)=JPREEDEF
      NPRECI(KUL)=JPINTDEF
      IF(CLLANG().EQ.'FRA') THEN
        LGLANG=.TRUE.
      ELSE
        LGLANG=.FALSE.
      ENDIF
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      IF(CDTYPO(1:ILTYPO).EQ.'R') THEN
        !
        ! READ.
        ! Le fichier existe-t-il?
        !
        INQUIRE(FILE=CDNOMF,EXIST=LLEX)
        IF(.NOT.LLEX) THEN
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFAOUV/ERREUR: fichier d''entree inexistant!...'
          ELSE
            PRINT*,'LFAOUV/ERROR: input file not found!...'
          ENDIF
          PRINT*,CDNOMF(1:ILNOMF)
          STOP 1
        ENDIF
        !
        ! Le fichier est-il bien LFA?
        !
        CALL LFATEST(KUL,CDNOMF,LLLFA)
        IF(.NOT.LLLFA) THEN
          IF(CLLANG().EQ.'FRA') THEN
            PRINT*,'LFAOUV/ERREUR: incompatibilite fichier/logiciel.'
          ELSE
            PRINT*,'LFAOUV/ERROR: file/software inconsistency!...'
          ENDIF
          PRINT*,CDNOMF(1:ILNOMF)
          STOP 1
        ENDIF
        !
        ! Le fichier existe et est bien LFA.
        ! On l'ouvre.
        !
        OPEN(KUL,FILE=CDNOMF,FORM='unformatted',STATUS='old')
        READ(KUL) IVERSIONESB
        IVERSION=IVERSIONESB
        !
        ! On sauvegarde la version qui a produit le fichier lu.
        !
        NVERSION(KUL)=IVERSION
        !
        ! Position du pointeur.
        !
        LGPOINT(KUL)=.FALSE.
      ELSEIF(CDTYPO(1:ILTYPO).EQ.'W') THEN
        !
        ! WRITE.
        !
        OPEN(KUL,FILE=CDNOMF,FORM='unformatted',STATUS='UNKNOWN')
        !
        ! Ecriture de l'en-tete du fichier LFA.
        !
        IVERSIONESB=IVERSION
        WRITE(KUL) IVERSIONESB
      ELSEIF(CDTYPO(1:ILTYPO).EQ.'S') THEN
        !
        ! SCRATCH.
        !
        OPEN(KUL,FORM='unformatted',STATUS='scratch')
        !
        ! Ecriture de l'en-tete du fichier LFA.
        !
        IVERSIONESB=IVERSION
        WRITE(KUL) IVERSIONESB
      ELSEIF(CDTYPO(1:ILTYPO).EQ.'A') THEN
        !
        ! APPEND.
        !
        OPEN(KUL,FILE=CDNOMF,FORM='unformatted' &
        ,STATUS='old',POSITION='append')
      ELSE
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,'LFAOUV/ERREUR: type d''ouverture inconnu!...'
        ELSE
          PRINT*,'LFAOUV/ERROR: unknown open type!...'
        ENDIF
        PRINT*,CDTYPO(1:ILTYPO)
        STOP 1
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAOUV',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAOUV
!
!
      SUBROUTINE LFAPPDEMO
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAPPDEMO* Demonstation du logiciel LFA.
      ! **** *LFAPPDEMO* LFA software demonstration.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      INTEGER(KIND=JPINTUSR) JPR,IR,JR,IUL,ILONG,IERR,ILC,IPREC &
      ,IULOUT
      CHARACTER*80 CLNOMF,CLNARC,CLNARL,CLNAEC,CLNAEL,CLNAC,CLNARA &
      ,CLNARD,CLNAED
      LOGICAL LLLFA,LLENT,LLREE,LLCAR
      PARAMETER(JPR=8) ! dimension physique.
      REAL(KIND=JPREEUSR) ZREEL(JPR),ZINTERM
      INTEGER(KIND=JPINTUSR) IENT(JPR)
      CHARACTER*200 CLC(JPR)
      CHARACTER*3 CLLANG
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPPDEMO',0,ZHOOK_HANDLE)
      IR=JPR/2 ! dimension logique.
      LLCAR=.TRUE.
      LLENT=.TRUE.
      LLREE=.TRUE.
      IF(CLLANG().EQ.'FRA') THEN
        CLNARD='Reels par defaut'
        CLNARC='Reels courts'
        CLNARL='Reels.Longs'
        CLNAED='Entiers defaut'
        CLNAEC='Entiers courts'
        CLNAEL='Entiers LONGS'
        CLNAC='Caracteres'
        CLNARA='Reels_ajoutes'
      ELSE
        CLNARD='Default REAL'
        CLNARC='Short REAL'
        CLNARL='Long REAL'
        CLNAED='Int. by default'
        CLNAEC='Short.integers'
        CLNAEL='Long_integers'
        CLNAC='Some char'
        CLNARA='Appended real'
      ENDIF
      !
      ! -------------------------------------------------
      ! Initialisation des tableaux.
      ! -------------------------------------------------
      !
      DO JR=1,IR
        ZINTERM=REAL(JR)
        ZREEL(JR)=SQRT(ZINTERM)
        IENT(JR)=JR*JR
        WRITE(CLC(JR),FMT='(a,i8)') 'Car ',JR*JR
      ENDDO
      !
      ! -------------------------------------------------
      ! Ouverture du fichier en mode ecriture.
      ! -------------------------------------------------
      !
      IUL=23
      CLNOMF='LFA'
      CALL LFAOUV(IUL,CLNOMF,'W')
      !
      ! -------------------------------------------------
      ! Niveau de messagerie.
      ! -------------------------------------------------
      !
      ! call lfames(iul,2)
      IF(LLREE) THEN
        !
        ! -------------------------------------------------
        ! Ecriture de reels.
        ! -------------------------------------------------
        !
        CALL LFAECRR(IUL,CLNARD,ZREEL,IR)
        !
        ! -------------------------------------------------
        ! Ecriture de reels courts.
        ! -------------------------------------------------
        !
        IPREC=4
        CALL LFAPRECR(IUL,IPREC)
        CALL LFAECRR(IUL,CLNARC,ZREEL,IR)
        !
        ! -------------------------------------------------
        ! Ecriture de reels longs.
        ! -------------------------------------------------
        !
        IPREC=8
        CALL LFAPRECR(IUL,IPREC)
        CALL LFAECRR(IUL,CLNARL,ZREEL,IR)
      ENDIF
      IF(LLCAR) THEN
        !
        ! -------------------------------------------------
        ! Ecriture de caracteres.
        ! -------------------------------------------------
        !
        CALL LFAECRC(IUL,CLNAC,CLC,IR)
      ENDIF
      IF(LLENT) THEN
        !
        ! -------------------------------------------------
        ! Ecriture d'entiers.
        ! -------------------------------------------------
        !
        CALL LFAECRI(IUL,CLNAED,IENT,IR)
        !
        ! -------------------------------------------------
        ! Ecriture d'entiers courts.
        ! -------------------------------------------------
        !
        IPREC=4
        CALL LFAPRECI(IUL,IPREC)
        CALL LFAECRI(IUL,CLNAEC,IENT,IR)
        !
        ! -------------------------------------------------
        ! Ecriture d'entiers longs.
        ! -------------------------------------------------
        !
        IPREC=8
        CALL LFAPRECI(IUL,IPREC)
        CALL LFAECRI(IUL,CLNAEL,IENT,IR)
      ENDIF
      !
      ! -------------------------------------------------
      ! Liste des articles.
      ! -------------------------------------------------
      !
      IULOUT=6
      CALL LFALAF(IUL,IULOUT)
      CALL LFAMINM(IUL)
      !
      ! -------------------------------------------------
      ! Fermeture du fichier.
      ! -------------------------------------------------
      !
      CALL LFAFER(IUL)
      !
      ! -------------------------------------------------
      ! Ouverture du fichier en mode ajout.
      ! -------------------------------------------------
      !
      IUL=23
      CLNOMF='LFA'
      CALL LFAOUV(IUL,CLNOMF,'A')
      IF(LLREE) THEN
        !
        ! -------------------------------------------------
        ! Ecriture de reels.
        ! -------------------------------------------------
        !
        CALL LFAECRR(IUL,CLNARA,ZREEL,IR)
      ENDIF
      !
      ! -------------------------------------------------
      ! Liste des articles.
      ! -------------------------------------------------
      !
      ! call lfaminm(iul)
      !
      ! -------------------------------------------------
      ! Fermeture du fichier.
      ! -------------------------------------------------
      !
      CALL LFAFER(IUL)
      !
      ! -------------------------------------------------
      ! Initialisation des tableaux.
      ! -------------------------------------------------
      !
      DO JR=1,JPR
        ZREEL(JR)=0.
        IENT(JR)=0
        CLC(JR)=' '
      ENDDO
      !
      ! -------------------------------------------------
      ! Ouverture du fichier en mode lecture.
      ! -------------------------------------------------
      !
      IUL=23
      CLNOMF='LFA'
      CALL LFATEST(IUL,CLNOMF,LLLFA)
      IF(LLLFA) THEN
        IF(LGLANG) THEN
          PRINT*,'Le fichier est de type LFA.'
        ELSE
          PRINT*,'The file is a LFA one.'
        ENDIF
      ELSE
        IF(LGLANG) THEN
          PRINT*,'Le fichier n''est pas de type LFA.'
        ELSE
          PRINT*,'The file is not a LFA one.'
        ENDIF
      ENDIF
      CALL LFAOUV(IUL,CLNOMF,'R')
      !
      ! -------------------------------------------------
      ! Niveau de messagerie.
      ! -------------------------------------------------
      !
      ! call lfames(iul,2)
      IF(LLENT) THEN
        !
        ! -------------------------------------------------
        ! Lecture d'entiers.
        ! -------------------------------------------------
        !
        CALL LFALECI(IUL,CLNAED,IR,IENT,ILONG,IERR)
        DO JR=1,IR
          PRINT*,'ient(',JR,')=',IENT(JR)
        ENDDO
      ENDIF
      IF(LLCAR) THEN
        !
        ! -------------------------------------------------
        ! Lecture de caracteres.
        ! -------------------------------------------------
        !
        CALL LFALECC(IUL,CLNAC,IR,CLC,ILONG,IERR)
        DO JR=1,IR
          CALL LONC(CLC(JR),ILC)
          PRINT*,CLC(JR)(1:ILC)
        ENDDO
      ENDIF
      IF(LLREE) THEN
        !
        ! -------------------------------------------------
        ! Lecture de reels.
        ! -------------------------------------------------
        !
        CALL LFALECR(IUL,CLNARD,IR,ZREEL,ILONG,IERR)
        DO JR=1,IR
          PRINT*,'zreel(',JR,')=',ZREEL(JR)
        ENDDO
      ENDIF
      !
      ! -------------------------------------------------
      ! Test de la gestion d'erreur.
      ! -------------------------------------------------
      !
      !
      ! Mise en mode "permissif".
      !
      ! call lfaerf(iul,.false.)
      !
      ! Lecture d'article inexistant.
      !
      ! call lfalecr(iul,'ACHSO',ir,zreel,ilong,ierr)
      ! print*,'Demande d''article inexistant: ierr=',ierr
      !
      ! -------------------------------------------------
      ! Fermeture du fichier.
      ! -------------------------------------------------
      !
      CALL LFAFER(IUL)
      IF (LHOOK) CALL DR_HOOK('LFAPPDEMO',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPPDEMO
!
!
      SUBROUTINE LFAPPLFAC
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAPPLFAC* Extraction sur output standard d'un article de fichier LFA.
      ! **** *LFAPPLFAC* Extract on standard output one LFA article.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   95-03, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      CHARACTER*80 CLNOMF,CLNOMC
      INTEGER(KIND=JPINTUSR) IUL,IULS,ILONG,IERR,ILNOMF
      CHARACTER*2 CLTYPE
      CHARACTER*3 CLLANG
      INTEGER(KIND=4) INARG
      !
      ! Saisie de la ligne de commande.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPPLFAC',0,ZHOOK_HANDLE)
      INARG=1
      CALL GETARGP(INARG,CLNOMF)
      INARG=2
      CALL GETARGP(INARG,CLNOMC)
      IF(CLNOMF.EQ.' '.OR.CLNOMC.EQ.' ') THEN
        IF(CLLANG().EQ.'FRA') THEN
          WRITE(*,*) ' '
          WRITE(*,*) &
          'Extraction sur output standard d''un article de fichier LFA.'
          WRITE(*,*) ' '
          WRITE(*,*) 'Utilisation: lfac nomf nomc '
          WRITE(*,*) 'avec'
          WRITE(*,*) '    nomf nom du fichier LFA d''entree.'
          WRITE(*,*) '    nomc nom du champ a extraire.'
          WRITE(*,*) ' '
          STOP
        ELSE
          WRITE(*,*) ' '
          WRITE(*,*) &
          'Extract on standard output one LFA article.'
          WRITE(*,*) ' '
          WRITE(*,*) 'Usage: lfac FILE ARTICLE'
          WRITE(*,*) 'with'
          WRITE(*,*) '    FILE: LFA file name.'
          WRITE(*,*) '    ARTICLE: article name in the file.'
          WRITE(*,*) ' '
          STOP
        ENDIF
      ENDIF
      !
      ! Ouverture du fichier.
      !
      ILNOMF=80
      IUL=7
      IULS=17
      CALL LFAOUV(IUL,CLNOMF,'R')
      !
      ! Recherche de l'article souhaite.
      !
      CALL LFACAS(IUL,CLNOMC,CLTYPE,ILONG,IERR)
      IF(IERR.EQ.0) THEN
        !
        ! L'article existe dans le fichier.
        !
        IF(CLTYPE(1:1).EQ.'R') THEN
          !
          ! Article de type reel.
          !
          CALL LFAAFFR(ILONG,IUL,CLNOMC)
        ELSEIF(CLTYPE(1:1).EQ.'I') THEN
          !
          ! Article de type entier.
          !
          CALL LFAAFFI(ILONG,IUL,CLNOMC)
        ELSEIF(CLTYPE(1:1).EQ.'C') THEN
          !
          ! Article de type caractere.
          !
          CALL LFAAFFC(ILONG,IUL,CLNOMC)
        ELSE
          PRINT*,'lfac: ERREUR interne: type de champ non prevu!...'
          STOP 1
        ENDIF
        CLOSE(IULS)
      ENDIF
      !
      ! Fermeture du fichier.
      !
      CALL LFAFER(IUL)
      IF (LHOOK) CALL DR_HOOK('LFAPPLFAC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPPLFAC
!
!
      SUBROUTINE LFAPPLFACRE
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAPPLFACRE* Creation d'un fichier LFA a partir de la ligne de commande.
      ! **** *LFAPPLFACRE* Create a LFA file from command line.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) JPARG,IULS,JARG,ILTYPE,IULE &
      ,INOMAL,ILVALFIC,IENTIER(1),IENTIERLOC,IGOL1,IPREC
      REAL(KIND=JPREEUSR) ZREEL(1),ZREELLOC
      LOGICAL LLEX
      CHARACTER*80 CLFS,CLNOMA,CLVALFIC
      CHARACTER*2 CLTYPE
      PARAMETER(JPARG=20)
      INTEGER(KIND=4) INARG
      CHARACTER*3 CLLANG
      !
      ! Saisie de la ligne de commande.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPPLFACRE',0,ZHOOK_HANDLE)
      INARG=1
      CALL GETARGP(INARG,CLFS)
      IF(CLFS.EQ.' ') THEN
        IF(CLLANG().EQ.'FRA') THEN
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') &
          'Creation d''un fichier LFA a partir de la ligne de commande'
          WRITE(*,FMT='(9a)') 'et(ou) de fichier(s) texte.'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') 'Utilisation:'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') &
          'lfacre LFA [nom_1 type_1 valfic_1] ... [nom_n type_n valfic_n]'
          WRITE(*,FMT='(a,i3)') 'n doit valoir au plus ',JPARG
          WRITE(*,FMT='(9a)') &
          'En sortie, le fichier LFA contiendra les n articles'
          WRITE(*,FMT='(9a)') &
          'nom_1 a nom_n, dont le type sera type_1 a type_n' &
          ,' (type: R4, R8, I4, I8 ou C),'
          WRITE(*,FMT='(9a)') &
          'et dont le contenu sera valfic_1 a valfic_n:'
          WRITE(*,FMT='(9a)') &
          '     - Si valfic_i est un fichier, alors le contenu'
          WRITE(*,FMT='(9a)') &
          '       de ce fichier sera le contenu de l''article nom_i.'
          WRITE(*,FMT='(9a)') &
          '     - Si valfic_i n''est pas un fichier, alors c''est la valeur'
          WRITE(*,FMT='(9a)') &
          '       du contenu de l''article de longueur 1 nom_i.'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') 'Exemple:'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') 'cat <<EOF > gol'
          WRITE(*,FMT='(9a)') 'gol1'
          WRITE(*,FMT='(9a)') 'gol2'
          WRITE(*,FMT='(9a)') 'EOF'
          WRITE(*,FMT='(9a)') &
          'lfacre LFA RII0 R8 1370. indice C gol annee I4 2006'
          WRITE(*,FMT='(9a)') &
          'creera le fichier LFA, contenant trois articles' &
          ,', l''article reel RII0'
          WRITE(*,FMT='(9a)') &
          '(longueur 1), l''article caractere indice (longueur 2),'
          WRITE(*,FMT='(9a)') &
          'et l''article entier annee (longueur 1).'
          WRITE(*,FMT='(9a)') ' '
          STOP
        ELSE
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') &
          'Create a LFA file from command line'
          WRITE(*,FMT='(9a)') 'and(or) from ASCII text file(s).'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') 'Usage:'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') &
          'lfacre LFA [article_name_1 type_1 fil_name_1] ' &
          ,'... [article_name_n type_n fil_name_n]'
          WRITE(*,FMT='(a,i3)') 'n has to be less than ',JPARG
          WRITE(*,FMT='(9a)') &
          'In output, the LFA file will contain the n articles'
          WRITE(*,FMT='(9a)') &
          'article_name_1 to article_name_n,' &
          ,' which type will be type_1 to type_n' &
          ,' (type: R4, R8, I4, I8 or C),'
          WRITE(*,FMT='(9a)') &
          'and contents of these articles will be ' &
          ,'fil_name_1 to fil_name_n:'
          WRITE(*,FMT='(9a)') &
          '     - If fil_name_i is a file, then its contents'
          WRITE(*,FMT='(9a)') &
          '       will be put in article_name_i article.'
          WRITE(*,FMT='(9a)') &
          '     - If fil_name_i is not a file, then it is the value'
          WRITE(*,FMT='(9a)') &
          '       of the one-value article article_name_i.'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') 'Example:'
          WRITE(*,FMT='(9a)') ' '
          WRITE(*,FMT='(9a)') 'cat <<EOF > gol'
          WRITE(*,FMT='(9a)') 'gol1'
          WRITE(*,FMT='(9a)') 'gol2'
          WRITE(*,FMT='(9a)') 'EOF'
          WRITE(*,FMT='(9a)') &
          'lfacre LFA RII0 R8 1370. indice C gol year I4 2006'
          WRITE(*,FMT='(9a)') &
          'will create the file LFA, containing tree articles' &
          ,', the real data article RII0'
          WRITE(*,FMT='(9a)') &
          '(length 1), the character data article indice (length 2),'
          WRITE(*,FMT='(9a)') &
          'and the integer data article year (length 1).'
          WRITE(*,FMT='(9a)') ' '
          STOP
        ENDIF
      ENDIF
      !
      ! Ouverture du fichier de sortie.
      !
      IULS=73
      CALL LFAOUV(IULS,CLFS,'W')
      !
      ! Boucle sur la ligne de commande.
      !
      DO JARG=1,JPARG
        INARG=2+3*(JARG-1)
        CALL GETARGP(INARG,CLNOMA)
        INARG=3+3*(JARG-1)
        CALL GETARGP(INARG,CLTYPE)
        INARG=4+3*(JARG-1)
        CALL GETARGP(INARG,CLVALFIC)
        IF(CLNOMA.EQ.' ') THEN
          !
          ! Fin de la ligne de commande.
          ! Fermeture des fichiers.
          !
          CALL LFAFER(IULS)
          STOP
        ELSE
          !
          ! Ligne de commande non vide.
          !
          CALL LONC(CLTYPE,ILTYPE)
          !
          ! -------------------------------------------------
          ! On force la précision d'écriture des entiers et réels.
          ! -------------------------------------------------
          !
          IF(CLTYPE(1:1).EQ.'R') THEN
            !
            ! Article réel.
            !
            READ(CLTYPE(2:2),FMT='(i1)') IPREC
            CALL LFAPRECR(IULS,IPREC)
          ELSEIF(CLTYPE(1:1).EQ.'I') THEN
            !
            ! Article entier.
            !
            READ(CLTYPE(2:2),FMT='(i1)') IPREC
            CALL LFAPRECI(IULS,IPREC)
          ENDIF
          !
          ! -------------------------------------------------
          ! L'utilisateur fournit-il un fichier d'entrée?
          ! -------------------------------------------------
          !
          INQUIRE(FILE=CLVALFIC,EXIST=LLEX)
          IF(LLEX) THEN
            !
            ! Le fichier existe.
            ! On va le lire sur un tableau.
            !
            ! Combien d'articles ce fichier comporte-t-il?
            !
            IULE=1
            OPEN(IULE,FILE=CLVALFIC,FORM='formatted')
            INOMAL=0
  100       CONTINUE
            READ(IULE,FMT='(a)',END=200)
            INOMAL=INOMAL+1
            GOTO 100
  200       CONTINUE
            REWIND(IULE)
            !
            ! Ce fichier comporte inomal articles.
            !
            IF(CLTYPE(1:1).EQ.'R') THEN
              !
              ! Article de reel.
              !
              CALL LFAFORVLR(IULE,IULS,CLTYPE,CLNOMA,INOMAL)
            ELSEIF(CLTYPE(1:1).EQ.'I') THEN
              !
              ! Article d'entier.
              !
              CALL LFAFORVLI(IULE,IULS,CLTYPE,CLNOMA,INOMAL)
            ELSEIF(CLTYPE(1:1).EQ.'C') THEN
              !
              ! Article caractere.
              !
              CALL LFAFORVLC(IULE,IULS,CLTYPE,CLNOMA,INOMAL)
            ELSE
              PRINT*,'LFACRE/ERREUR: type d''article inconnu!...'
              PRINT*,CLTYPE
              STOP 1
            ENDIF
            CLOSE(IULE)
          ELSE
            !
            ! Le fichier n'existe pas.
            ! On ecrit un article de longueur 1.
            !
            IGOL1=1
            IF(CLTYPE(1:1).EQ.'R') THEN
              !
              ! Article de reel.
              !
              CALL LONC(CLVALFIC,ILVALFIC)
              CALL CARREE(CLVALFIC,ILVALFIC,ZREELLOC)
              ZREEL(1)=ZREELLOC
              CALL LFAECRR(IULS,CLNOMA,ZREEL,IGOL1)
            ELSEIF(CLTYPE(1:1).EQ.'I') THEN
              !
              ! Article d'entier.
              !
              CALL LONC(CLVALFIC,ILVALFIC)
              IGOL1=1
              CALL CARINT(CLVALFIC,IGOL1,ILVALFIC,IENTIERLOC)
              IENTIER(1)=IENTIERLOC
              CALL LFAECRI(IULS,CLNOMA,IENTIER,IGOL1)
            ELSEIF(CLTYPE(1:1).EQ.'C') THEN
              !
              ! Article caractere.
              !
              CALL LFAECRC(IULS,CLNOMA,CLVALFIC,IGOL1)
            ELSE
              PRINT*,'LFACRE/ERREUR: type d''article inconnu!...'
              PRINT*,CLTYPE
              STOP 1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAPPLFACRE',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPPLFACRE
!
!
      SUBROUTINE LFAPPLFALAF
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAPPLFALAF* Programme listant les articles de fichiers ASCII autodocumentes.
      ! **** *LFAPPLFALAF* Get the articles list of a LFA file.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   94-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      CHARACTER*80 CLNOMF
      INTEGER(KIND=JPINTUSR) IUL,IULOUT
      INTEGER(KIND=4) INARG
      CHARACTER*3 CLLANG
      !
      ! Ouverture du fichier.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPPLFALAF',0,ZHOOK_HANDLE)
      INARG=1
      CALL GETARGP(INARG,CLNOMF)
      IF(CLNOMF.EQ.' ') THEN
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,' '
          PRINT*,'Sortie sur output standard de la liste des articles' &
          ,' d''un fichier LFA.'
          PRINT*,' '
          PRINT*,'Utilisation: lfalaf nomf'
          PRINT*,' '
          STOP
        ELSE
          PRINT*,' '
          PRINT*,'Get the articles list of a LFA file.'
          PRINT*,' '
          PRINT*,'Usage: lfalaf FILE'
          PRINT*,' '
          STOP
        ENDIF
      ENDIF
      IUL=7
      CALL LFAOUV(IUL,CLNOMF,'R')
      !
      ! Contenu du fichier.
      !
      IULOUT=6
      CALL LFALAF(IUL,IULOUT)
      !
      ! Fermeture du fichier.
      !
      CALL LFAFER(IUL)
      IF (LHOOK) CALL DR_HOOK('LFAPPLFALAF',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPPLFALAF
!
!
      SUBROUTINE LFAPPLFAMINM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAPPLFAMINM* Sortie des extrema des articles de fichiers LFA.
      ! **** *LFAPPLFAMINM* Print out extrema of all articles from some LFA files.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) IARG,IUL,ILF
      INTEGER(KIND=4) JARG,IARG4,IARGC
      CHARACTER*80 CLF
      CHARACTER*3 CLLANG
      !
      ! Saisie de la ligne de commande.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPPLFAMINM',0,ZHOOK_HANDLE)
      IARG4=IARGC()
      IARG=IARG4
      IF(IARG.EQ.0) THEN
        IF(CLLANG().EQ.'FRA') THEN
          WRITE(*,'(a)') ' '
          WRITE(*,'(2a)') 'Sortie des extrema et moyenne des articles' &
          ,' d''un (plusieurs) fichier(s) LFA.'
          WRITE(*,'(a)') ' '
          WRITE(*,'(a)') 'Utilisation: lfaminm LFA1 [LFA2 ... LFAn]'
          WRITE(*,'(a)') ' '
          STOP
        ELSE
          WRITE(*,'(a)') ' '
          WRITE(*,'(2a)') 'Prints out extrema, mean and rms ' &
          ,'of all articles from one (or more) LFA file(s).'
          WRITE(*,'(a)') ' '
          WRITE(*,'(a)') 'Usage: lfaminm LFA1 [LFA2 ... LFAn]'
          WRITE(*,'(a)') ' '
          STOP
        ENDIF
      ENDIF
      DO JARG=1,IARG
        CALL GETARGP(JARG,CLF)
        !
        ! Ouverture du fichier.
        !
        IUL=72
        CALL LFAOUV(IUL,CLF,'R')
        CALL LONC(CLF,ILF)
        !
        ! Determination des extrema, avec sortie sur output standard.
        !
        CALL LFAMINM(IUL)
        !
        ! Fermeture du fichier.
        !
        CALL LFAFER(IUL)
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LFAPPLFAMINM',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPPLFAMINM
!
!
      SUBROUTINE LFAPPLFAREU
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAPPLFAREU* Reunion de deux fichiers de LFA.
      ! **** *LFAPPLFAREU* Fuse two LFA files.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) JPLIS,IULE1,IULE2,IULS,ILIS1 &
      ,ILIS2,JLIS1,ILLIS1,JLIS2,ILLIS2
      LOGICAL LLINTER
      CHARACTER*80 CLFE1,CLFE2,CLFS
      PARAMETER(JPLIS=600) ! nombre maxi d'articles dans les fichiers LFA d'entree.
      CHARACTER*80 CLLIS1(JPLIS) ! liste des articles du fichier d'entree 1.
      CHARACTER*80 CLLIS2(JPLIS) ! liste des articles du fichier d'entree 2.
      INTEGER(KIND=4) INARG
      CHARACTER*3 CLLANG
      !
      ! Saisie de la ligne de commande.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPPLFAREU',0,ZHOOK_HANDLE)
      INARG=1
      CALL GETARGP(INARG,CLFE1)
      INARG=2
      CALL GETARGP(INARG,CLFE2)
      INARG=3
      CALL GETARGP(INARG,CLFS)
      IF(CLFS.EQ.' ') THEN
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,' '
          PRINT*,'Reunion de deux fichiers LFA.'
          PRINT*,' '
          PRINT*,'Utilisation: lfareu F1 F2 Fres'
          PRINT*,'      En entree: F1 et F2, en sortie: Fres.'
          PRINT*,'      F2 est prioritaire sur F1, i.e. si un article'
          PRINT*,'      est present dans F1 et F2,' &
          ,' c''est celui de F2 qui sera copie.'
          PRINT*,' '
          STOP
        ELSE
          PRINT*,' '
          PRINT*,'Fuse two LFA files.'
          PRINT*,' '
          PRINT*,'Usage: lfareu F1 F2 Fres'
          PRINT*,'      In input: F1 and F2, in output: Fres.'
          PRINT*,'      F2 has higher priority than F1, i.e. if an article'
          PRINT*,'      is present in both F1 and F2,' &
          ,' the article from F2 will be copied.'
          PRINT*,' '
          STOP
        ENDIF
      ENDIF
      !
      ! Ouverture des fichiers.
      !
      IULE1=72
      IULE2=74
      IULS=73
      CALL LFAOUV(IULE1,CLFE1,'R')
      CALL LFAOUV(IULE2,CLFE2,'R')
      CALL LFAOUV(IULS,CLFS,'W')
      !
      ! Liste des articles.
      !
      CALL LFALAFT(IULE1,CLLIS1,JPLIS,ILIS1)
      CALL LFALAFT(IULE2,CLLIS2,JPLIS,ILIS2)
      !
      ! On boucle sur les articles du fichier 1.
      ! On copie en sortie chaque article,
      ! en le lisant sur 1 s'il n'est que dans 1,
      ! et en le lisant sur 2 s'il est dans les deux.
      !
      DO JLIS1=1,ILIS1
        LLINTER=.FALSE. ! vrai si l'article est present dans l'intersection des deux fichiers.
        CALL LONC(CLLIS1(JLIS1),ILLIS1)
        DO JLIS2=1,ILIS2
          CALL LONC(CLLIS2(JLIS2),ILLIS2)
          IF(CLLIS2(JLIS2)(1:ILLIS2).EQ.CLLIS1(JLIS1)(1:ILLIS1)) &
          LLINTER=.TRUE.
        ENDDO
        IF(LLINTER) THEN
          !
          ! On copie l'article du fichier 2 au fichier de sortie.
          !
          CALL LFACOP(IULE2,CLLIS1(JLIS1),' ',IULS)
        ELSE
          !
          ! On copie l'article du fichier 1 au fichier de sortie.
          !
          CALL LFACOP(IULE1,CLLIS1(JLIS1),' ',IULS)
        ENDIF
      ENDDO
      !
      ! On boucle sur les articles du fichier 2.
      ! On copie en sortie les articles,
      ! de 2 qui ne sont pas dans 1.
      !
      DO JLIS2=1,ILIS2
        LLINTER=.FALSE. ! vrai si l'article est present dans l'intersection des deux fichiers.
        CALL LONC(CLLIS2(JLIS2),ILLIS2)
        DO JLIS1=1,ILIS1
          CALL LONC(CLLIS1(JLIS1),ILLIS1)
          IF(CLLIS2(JLIS2)(1:ILLIS2).EQ.CLLIS1(JLIS1)(1:ILLIS1)) &
          LLINTER=.TRUE.
        ENDDO
        IF(.NOT.LLINTER) THEN
          !
          ! L'article est seulement dans le 2.
          !
          CALL LFACOP(IULE2,CLLIS2(JLIS2),' ',IULS)
        ENDIF
      ENDDO
      !
      ! Fermeture des fichiers.
      !
      CALL LFAFER(IULE1)
      CALL LFAFER(IULE2)
      CALL LFAFER(IULS)
      IF (LHOOK) CALL DR_HOOK('LFAPPLFAREU',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPPLFAREU
!
!
      SUBROUTINE LFAPPLFATEST
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFAPPLFATEST* Programme testant si un fichier est LFA.
      ! **** *LFAPPLFATEST* Prints out if a given file is a LFA one.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   98-01, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      CHARACTER*80 CLNOMF
      INTEGER(KIND=JPINTUSR) IUL
      LOGICAL LLLFA
      INTEGER(KIND=4) INARG
      CHARACTER*3 CLLANG
      !
      ! Ouverture du fichier.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPPLFATEST',0,ZHOOK_HANDLE)
      INARG=1
      CALL GETARGP(INARG,CLNOMF)
      IF(CLNOMF.EQ.' ') THEN
        IF(CLLANG().EQ.'FRA') THEN
          PRINT*,' '
          PRINT*,'Programme testant si un fichier est de type LFA.'
          PRINT*,' '
          PRINT*,'Utilisation: lfatest FICHIER'
          PRINT*,' '
          PRINT*,'La reponse, sur output standard, sera'
          PRINT*,'      - "lfa"     si le fichier est LFA,'
          PRINT*,'      - "non-lfa" sinon.'
          PRINT*,' '
          STOP
        ELSE
          PRINT*,' '
          PRINT*,'Prints out if a given file is a LFA one.'
          PRINT*,' '
          PRINT*,'Usage: lfatest FILE'
          PRINT*,' '
          PRINT*,'The result, on standard output, will be'
          PRINT*,'      - "lfa"     if the file is a LFA one,'
          PRINT*,'      - "non-lfa" else case.'
          PRINT*,' '
          STOP
        ENDIF
      ENDIF
      IUL=7
      CALL LFATEST(IUL,CLNOMF,LLLFA)
      IF(LLLFA) THEN
        WRITE(*,'(a)') 'lfa'
      ELSE
        WRITE(*,'(a)') 'non-lfa'
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAPPLFATEST',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPPLFATEST
!
!
      SUBROUTINE LFAPRECI(KUL,KPREC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAPRECI* Forcage de la precision d'ecriture des entiers.
      ! **** *LFAPRECI* Force integer data writing precision.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-02, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kul      unite logique du fichier LFA.
      ! kprec    precision des entiers a ecrire ulterieurement, en octets.
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! kul      logical unit of LFA file.
      ! kprec    precision of integer data to write, in bytes.
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      INTEGER(KIND=JPINTUSR) KUL,KPREC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPRECI',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfapreci: forcage de la precision ',KPREC
      ENDIF
      IF(KPREC.EQ.4.OR.KPREC.EQ.8) THEN
        NPRECI(KUL)=KPREC
      ELSE
        PRINT*,'LFAPRECI/ERREUR: precision en octets' &
        ,' non recevable!...',KPREC
        STOP 1
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAPRECI',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPRECI
!
!
      SUBROUTINE LFAPRECR(KUL,KPREC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAPRECR* Forcage de la precision d'ecriture des reels.
      ! **** *LFAPRECR* Force real data writing precision.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-02, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kul      unite logique du fichier LFA.
      ! kprec    precision des reels a ecrire ulterieurement, en octets.
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! kul      logical unit of LFA file.
      ! kprec    precision of real data to write, in bytes.
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      INTEGER(KIND=JPINTUSR) KUL,KPREC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAPRECR',0,ZHOOK_HANDLE)
      IF(NMES(KUL).EQ.2) THEN
        !
        ! Messagerie bavarde.
        !
        PRINT*,'++ lfaprecr: forcage de la precision ',KPREC
      ENDIF
      IF(KPREC.EQ.4.OR.KPREC.EQ.8) THEN
        NPRECR(KUL)=KPREC
      ELSE
        PRINT*,'LFAPRECR/ERREUR: precision en octets' &
        ,' non recevable!...',KPREC
        STOP 1
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFAPRECR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAPRECR
!
!
      SUBROUTINE LFAREW(KUL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *LFAREW* Rebobinage d'un fichier LFA.
      ! **** *LFAREW* Rewind a LFA file.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   97-10, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! kul: unite logique du fichier LFA.
      ! En sortie:
      ! --------------------------------------------------------------
      ! Input:
      ! kul: logical unit of LFA file.
      ! En sortie:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL
      INTEGER(KIND=JPINTUSR) IVERSION
      INTEGER(KIND=JPINTESB) IVERSIONESB
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      !
      ! Rebobinage.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFAREW',0,ZHOOK_HANDLE)
      REWIND(KUL)
      !
      ! Lecture de l'en-tete..
      !
      READ(KUL) IVERSIONESB
      IVERSION=IVERSIONESB
      !
      ! Position du pointeur.
      !
      LGPOINT(KUL)=.FALSE.
      IF (LHOOK) CALL DR_HOOK('LFAREW',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFAREW
!
!
      SUBROUTINE LFATEST(KUL,CDNOMF,LDLFA)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LFATEST* Teste si un fichier est bien de type LFA.
      ! **** *LFATEST* Test if a file is a LFA one.
      ! --------------------------------------------------------------------------
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
      ! Auteur:   97-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En entree:
      ! kul         unite logique du fichier;
      ! .           ce doit etre une unite disponible:
      ! .           le fichier va etre ouvert sous cette unite logique.
      ! cdnomf      nom du fichier.
      ! En sortie:
      ! ldlfa=.true. si le fichier est de type LFA, .false. sinon.
      ! --------------------------------------------------------------------------
      ! Input:
      ! kul         logical unit of file.
      ! .           this unit has to be free:
      ! .           the file will be opened with this logical unit.
      ! cdnomf      file name.
      ! Output:
      ! ldlfa=.true. if the file is a LFA one, .false. else case.
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KUL
      INTEGER(KIND=JPINTUSR) IVERSION
      INTEGER(KIND=JPINTESB) IVERSIONESB
      LOGICAL LDLFA,LLEX
      !
      ! -------------------------------------------------
      ! Parametres.
      ! -------------------------------------------------
      !
      INTEGER JPNOMA
      PARAMETER(JPNOMA=80) ! nombre maximal de caractères des noms d'article.
      INTEGER JPERF
      PARAMETER(JPERF=100) ! nombre maximal de fichiers ouvrables.
      !
      ! -------------------------------------------------
      ! Tableaux globaux renseignant la position
      ! du pointeur dans le fichier:
      ! lgpoint est
      ! .    - faux si le pointeur est avant une autodocumentation
      ! .    - vrai s'il est avant un article de donnees.
      ! Dans le cas ou lgpoint est vrai:
      ! .    - ntype est le type de l'article.
      ! .    - nlong est sa longueur.
      ! .    - cgna est son nom.
      ! Dans le cas ou lgpoint est faux ces 3 tableaux
      ! .    ne sont pas utilises.
      ! -------------------------------------------------
      !
      !
      ! -------------------------------------------------
      ! Logiques.
      ! -------------------------------------------------
      !
      !
      ! lgerf(kul): .true. si toute erreur sur le fichier kul
      ! doit etre fatale.
      LOGICAL LGERF(JPERF)
      LOGICAL LGPOINT(JPERF)
      LOGICAL LGLANG ! .true. if french, .false. if english.
      COMMON/LFACOML/LGERF,LGPOINT,LGLANG
      !
      ! -------------------------------------------------
      ! Entiers.
      ! -------------------------------------------------
      !
      !
      ! nmes(kul): niveau de messagerie associe au fichier
      ! d'unite logique kul:
      ! 0 rien sauf erreurs fatales.
      ! 1 erreurs fatales + messages d'attention.
      ! 2 bavard (mode utile en développement du logiciel uniquement).
      !
      INTEGER NMES(JPERF)
      !
      ! nversion: version du logiciel qui a produit le fichier lu.
      !
      INTEGER NVERSION(JPERF)
      INTEGER NTYPE(JPERF)
      INTEGER NLONG(JPERF)
      !
      ! Precision des reels et entiers en octets.
      !
      INTEGER NPRECR(JPERF)
      INTEGER NPRECI(JPERF)
      COMMON/LFACOMI/NMES,NTYPE,NLONG,NVERSION,NPRECR,NPRECI
      !
      ! -------------------------------------------------
      ! Caracteres.
      ! -------------------------------------------------
      !
      ! cgnomf(kul): nom en clair du fichier d'unite logique kul.
      CHARACTER*80 CGNOMF(JPERF)
      !
      ! cgtypo(kul): type d'ouverture du fichier: 'R' READ, 'W' WRITE, 'S' scratch.
      CHARACTER*1 CGTYPO(JPERF)
      CHARACTER*(JPNOMA) CGNA(JPERF)
      COMMON/LFACOMC/CGNOMF,CGTYPO,CGNA
      CHARACTER*(*) CDNOMF
      !
      ! Le fichier existe-t-il?
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFATEST',0,ZHOOK_HANDLE)
      INQUIRE(FILE=CDNOMF,EXIST=LLEX)
      IF(.NOT.LLEX) THEN
        !
        ! Le fichier n'existe pas.
        ! On se contente de retourner qu'il est non-LFA.
        !
        LDLFA=.FALSE.
      ELSE
        !
        ! Le fichier existe.
        ! On lit les premiers caracteres afin de voir s'ils
        ! sont conformes a ce qu'on est en droit d'attendre
        ! d'un LFA!...
        !
        OPEN(KUL,FILE=CDNOMF,FORM='unformatted',STATUS='old',ERR=100)
        READ(KUL,ERR=100,END=100) IVERSIONESB
        IVERSION=IVERSIONESB
        CLOSE(KUL)
        IF(IVERSION.EQ.8.OR.IVERSION.EQ.9) THEN
          !
          ! L'en-tete lue est bien celle
          ! attendue.
          !
          LDLFA=.TRUE.
        ELSE
          !
          ! L'en-tete lue n'est pas celle
          ! attendue.
          !
          LDLFA=.FALSE.
        ENDIF
        IF (LHOOK) CALL DR_HOOK('LFATEST',1,ZHOOK_HANDLE)
        RETURN
  100   CONTINUE
        !
        ! Si on est ici, c'est qu'il y a eu erreur
        ! lors de la lecture de l'en-tete sur le fichier.
        ! Le fichier n'est donc pas un LFA.
        !
        CLOSE(KUL)
        LDLFA=.FALSE.
      ENDIF
      IF (LHOOK) CALL DR_HOOK('LFATEST',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LFATEST
!
!
      SUBROUTINE LFAIUNASSIGN(KUL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
      IMPLICIT NONE
!
      INTEGER :: KUL
!RJ: where is assign?
      ENDSUBROUTINE LFAIUNASSIGN
!
!
      SUBROUTINE GETARGP(KARG,CDARG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *getargp* GET Command Line Arguments.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Cette routine est creee pour contourner une bug du 1 HP,
      ! qui decale de un par rapport au standard l'indice de l'argument
      ! a getarg!...
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! Input:
      ! Output:
      ! --------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 KARG,IARG
      CHARACTER*(*) CDARG
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('GETARGP',0,ZHOOK_HANDLE)
      CALL GETARG(KARG,CDARG)
      IF (LHOOK) CALL DR_HOOK('GETARGP',1,ZHOOK_HANDLE)
      ENDSUBROUTINE GETARGP
!
!
      SUBROUTINE LONC(CDCHA,KREP)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *LONC*  - Reperage de la longueur significative de la chaine.
      ! --------------------------------------------------------------------------
      ! Sujet:
      ! ------
      ! Arguments explicites: /
      ! ---------------------
      ! Arguments implicites: /
      ! ---------------------
      ! Methode:
      ! --------
      ! Externes: /
      ! ---------
      ! Auteur:  J.M.Piriou
      ! -------
      ! Original : 91-11
      ! ----------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! En Entree       :
      ! cdCha       : chaine a estimer
      ! En Sortie       :
      ! krep        : position du premier non blanc a partir de la droite.
      ! --------------------------------------------------------------------------
      ! Exemple:
      ! cdCha       : 'CHAINE    '
      ! =====>
      ! krep        : 6
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KREP,JTEST,LEN
      CHARACTER *(*) CDCHA
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LONC',0,ZHOOK_HANDLE)
      KREP=0
      DO JTEST=LEN(CDCHA),1,-1
        IF(CDCHA(JTEST:JTEST).NE.' ') THEN
          KREP=JTEST
          IF (LHOOK) CALL DR_HOOK('LONC',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
      ENDDO
      IF (LHOOK) CALL DR_HOOK('LONC',1,ZHOOK_HANDLE)
      ENDSUBROUTINE LONC
!
!
      SUBROUTINE CARINT(CDENT,KPOS1,KPOS2,KRES)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *CARINT *  - Conversion chaine de caractere > entier
      ! --------------------------------------------------------------------------
      ! Auteur:  J.M.Piriou
      ! -------
      ! Modifications.
      ! --------------
      ! Original : 91-01
      ! --------------------------------------------------------------------------
      ! En entree: cdent     chaine de caracteres (l.max. 80 caracteres)
      ! kpos1,kpos2 abscisses de debut de fin d'extraction dans cdent:
      ! on extrait du caractere kpos1 au caractere kpos2
      ! En sortie: kres  resultat entier
      ! Remarque:  les seuls caracteres numeriques sont traites
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KPOS1,KPOS2,KRES,ISIGN,JPOS,IASCII
      CHARACTER *80 CDENT
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CARINT',0,ZHOOK_HANDLE)
      KRES=0
      ISIGN=1
      DO JPOS=KPOS1,KPOS2
        IASCII=ICHAR(CDENT(JPOS:JPOS))
        IF ((IASCII.GE.48).AND.(IASCII.LE.57)) THEN
          !
          ! Chiffre entre 0 et 9.
          !
          KRES=10*KRES+IASCII-48
        ELSEIF(IASCII.EQ.45) THEN
          !
          ! Signe -.
          !
          ISIGN=-ISIGN
        ENDIF
      ENDDO
      KRES=ISIGN*KRES
      IF (LHOOK) CALL DR_HOOK('CARINT',1,ZHOOK_HANDLE)
      ENDSUBROUTINE CARINT
!
!
      SUBROUTINE CARREE(CDENT,KLENT,PREEL)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *CARREE *  Extraction d'un reel d'une chaine de caracteres.
      ! --------------------------------------------------------------------------
      ! Sujet:
      ! ------
      ! Arguments explicites: /
      ! ---------------------
      ! Arguments implicites: /
      ! ---------------------
      ! Methode:
      ! --------
      ! Externes: /
      ! ---------
      ! Auteur:           92-10, J.M. Piriou.
      ! -------
      ! Modifications:
      ! --------------------------------------------------------------------------
      ! en entree: cdent chaine contenant le reel
      ! klent nombre de caracteres significatifs ecrits sur cdent
      ! en sortie: preel reel lu
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KLENT,IPOSPNB,JENT
      REAL(KIND=JPREEUSR) PREEL
      CHARACTER*(*) CDENT
      CHARACTER*01 CL1
      !
      ! Test de validite de l'entree.
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CARREE',0,ZHOOK_HANDLE)
      IPOSPNB=KLENT+1 ! position du premier caractère non blanc.
      DO JENT=1,KLENT
        CL1=CDENT(JENT:JENT)
        IF(CL1.EQ.' ') THEN
          ! Blanc.
        ELSEIF((CL1.LE.'9'.AND.CL1.GE.'0').OR.CL1.EQ.'-' &
        .OR.CL1.EQ.'+'.OR.CL1.EQ.'E'.OR.CL1.EQ.'e' &
        .OR.CL1.EQ.'.') THEN
          ! Partie de nombre.
          IPOSPNB=MIN(IPOSPNB,JENT)
        ELSE
          PRINT*,'CARREE/ERREUR: extraction de reel' &
          ,' impossible sur la chaine ',CDENT(1:KLENT)
          PRINT*,'Caractere incrimine: ''',CL1,''''
          STOP 1
        ENDIF
      ENDDO
  200 CONTINUE
      !
      ! Extraction du reel.
      !
      IF(CDENT.NE.' ') THEN
        READ(CDENT(IPOSPNB:KLENT),FMT=*) PREEL
      ELSE
        PREEL=0.
      ENDIF
      IF (LHOOK) CALL DR_HOOK('CARREE',1,ZHOOK_HANDLE)
      ENDSUBROUTINE CARREE
!
!
      SUBROUTINE REECAR(PX,KOPT,KNC,CDSOR,KSOR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------------------
      ! **** *REECAR *  - Conversion reel > chaine de caracteres.
      ! --------------------------------------------------------------------------
      !
      ! Sujet:
      ! ------
      !
      ! Arguments explicites: /
      ! ---------------------
      !
      ! Arguments implicites: /
      ! ---------------------
      !
      ! Methode:
      ! --------
      !
      ! Externes: /
      ! ---------
      !
      ! Original:
      ! ---------
      ! 92-05-19, J.M. Piriou
      !
      ! Modifications:
      ! --------------
      !
      ! --------------------------------------------------------------------------
      ! En Entree       :
      ! px          : reel a convertir
      ! kopt: option:
      ! si kopt=-1 on veut knc chiffres significatifs.
      ! si kopt=-2 on veut un affichage fixe ou scientifique
      ! ...........suivant la valeur absolue du réel
      ! ...........(fixe pour les valeurs absolues "proches" de 1).
      ! ...........Si knc=3, 3 chiffres significatifs garantis, 7 sinon.
      ! si kopt>=0 et knc=0 affichage a kopt chiffres
      ! ........... apres la virgule.
      ! si kopt>=0 et knc>0 affichage a kopt chiffres
      ! ...........apres la virgule, et si de plus px est voisin d'un entier
      ! ...........a mieux que 10.**(-knc), on l'affiche
      ! ...........au format entier.
      ! En Sortie       :
      ! cdsor         : chaine de caracteres resultante
      ! ksor          : nombre de caracteres ecrits sur cdsor
      ! --------------------------------------------------------------------------
      ! Exemples:
      ! px          kopt        knc         cdsor       ksor
      ! 6.22        -1          2           '6.2'       3
      ! 1.E-6       -1          3           '1.00E-6'   8
      ! 1.E-6       -1          indifférent '1.00E-6'   8
      ! -6.22       1           0           '-6.2'      4
      ! 6.00012     2           0           '6.00'      4
      ! 6.00012     2           2           '6'         1
      ! 1.E-6       3           0           '0.000'     5
      ! 1.E-6       3           5           '0'         1
      ! .1          2           4           '0.10'      4
      ! --------------------------------------------------------------------------
      IMPLICIT NONE
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention de passer
      ! en argument au logiciel LFA.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEUSR=JPRB ! reels
      INTEGER, PARAMETER :: JPINTUSR=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers et reels
      ! que l'utilisateur a l'intention d'ecrire par defaut
      ! sur les fichiers.
      ! Cette taille peut etre modifiee
      ! lors de l'execution, par appel a lfaprec.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPREEDEF=4 ! reels
      INTEGER, PARAMETER :: JPINTDEF=4 ! entiers
      !
      ! -------------------------------------------------
      ! Taille en octets des entiers a usage interne LFA,
      ! et qui renseignent dans le fichier
      ! la version du logiciel, la taille des articles, etc...
      ! Afin que les fichiers soient relisibles
      ! entre les differentes versions du logiciel LFA
      ! et les differentes precisions possibles,
      ! NE PAS CHANGER la valeur 4 ci-dessous.
      ! -------------------------------------------------
      !
      INTEGER, PARAMETER :: JPINTESB=4 ! entiers
      INTEGER(KIND=JPINTUSR) KOPT,KNC,KSOR,ISIGN,ILOG &
      ,ILONLOG,IX,IY,IAFF,ILOG2,IAFFTMP,JCH
      REAL(KIND=JPREEUSR) PX,ZX,ZLOG,ZMANT,Z10,ZVALA,ZXNOR,ZINT
      CHARACTER *(*) CDSOR
      CHARACTER *144  CLXPROV
      CHARACTER *11 CLMANT
      CHARACTER *10 CLMANT2
      CHARACTER *04 CLLOG
      CHARACTER *22 CLFOR
      CHARACTER *08 CLFORF
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('REECAR',0,ZHOOK_HANDLE)
      IF(KOPT.EQ.-1) THEN
        ! On veut un certain nombre de chiffres SIGNIFICATIFS.
        IF(PX.EQ.0.) THEN
          ! Cas nombre nul
          CDSOR='0.'
          KSOR=2
        ELSE
          ! Cas general
          ZX=PX
          IF(ZX.LT.0.) THEN
            ZX=-ZX
            ISIGN=-1
          ELSE
            ISIGN=1
          ENDIF
          ZLOG=LOG(ZX)/LOG(10.)
          IF(ZLOG.LT.0.) THEN
            ILOG=INT(ZLOG)-1
          ELSE
            ILOG=INT(ZLOG)
          ENDIF
          ! On calcule la mantisse de zx, soit un nombre
          ! zmant tel que 1.<=zmant<10.
          ZMANT=ZX/10.**REAL(ILOG)
          ! On ne conserve que knc chiffres significatifs
          Z10=10.**REAL(KNC-1)
          ZMANT=NINT(ZMANT*Z10)/Z10
          WRITE(CLMANT,FMT=1000) ZMANT
          IF(CLMANT(1:1).EQ.'*') THEN
            ! Il y a eu un probleme lors de l'arrondi
            ! permettant l'obtention de la mantisse de px,
            ! ce pour les valeurs de mantisse voisines de 1. ou 10.
            ! On reecrit ici le resultat.
            CLMANT='1.000000000'
            IF(ZMANT.LT.5.) THEN
              ILOG=ILOG-1
            ELSE
              ILOG=ILOG+1
            ENDIF
          ENDIF
          IF(ILOG.LE.-100) THEN
            WRITE(CLLOG,FMT='(i4)') ILOG
            ILONLOG=4
          ELSEIF(ILOG.LE.-10) THEN
            WRITE(CLLOG,FMT='(i3)') ILOG
            ILONLOG=3
          ELSEIF(ILOG.LT.0) THEN
            WRITE(CLLOG,FMT='(i2)') ILOG
            ILONLOG=2
          ELSEIF(ILOG.LT.10) THEN
            WRITE(CLLOG,FMT='(i1)') ILOG
            ILONLOG=1
          ELSEIF(ILOG.LT.100.) THEN
            WRITE(CLLOG,FMT='(i2)') ILOG
            ILONLOG=2
          ELSE
            WRITE(CLLOG,FMT='(i3)') ILOG
            ILONLOG=3
          ENDIF
          CLMANT2=CLMANT(1:1)//CLMANT(3:11)
          IF(ILOG.GE.5) THEN
            !
            ! Cas des nombres x <= 10**5.
            !
            CLXPROV=CLMANT(1:KNC+1)//'E'//CLLOG
            KSOR=KNC+1+1+ILONLOG
          ELSEIF(ILOG.GE.0) THEN
            !
            ! Cas des nombres 1.<= x < 10**5.
            !
            IF(KNC.GE.ILOG+2) THEN
              CLXPROV=CLMANT2(1:ILOG+1)//'.' &
              //CLMANT2(ILOG+2:KNC)
              KSOR=KNC+1
            ELSE
              CLXPROV=CLMANT2(1:ILOG+1)//'.'
              KSOR=ILOG+2
            ENDIF
          ELSEIF(ILOG.GE.-1) THEN
            CLXPROV='0.'//CLMANT2(1:KNC)
            KSOR=KNC+2
          ELSEIF(ILOG.GE.-2) THEN
            CLXPROV='0.0'//CLMANT2(1:KNC)
            KSOR=KNC+3
          ELSE
            CLXPROV=CLMANT(1:KNC+1)//'E'//CLLOG
            KSOR=KNC+1+1+ILONLOG
          ENDIF
          IF(ISIGN.EQ.-1) THEN
            CDSOR='-'//CLXPROV(1:KSOR)
            KSOR=KSOR+1
          ELSE
            CDSOR=CLXPROV(1:KSOR)
          ENDIF
        ENDIF
      ELSEIF(KOPT.EQ.-2) THEN
        !
        ! On veut un affichage fixe ou scientifique
        ! suivant la valeur absolue du réel.
        !
        IF(KNC.EQ.3) THEN
          !
          ! On affiche 3 chiffres significatifs.
          !
          ZVALA=ABS(PX)
          ZX=PX
          IF(ZVALA.EQ.0.) THEN
            !
            ! Nombre nul.
            ! Choix du format fixe.
            !
            CLFOR='(F11.4)'
            KSOR=11
          ELSEIF(ZVALA.GE.100000..OR.ZVALA.LT.0.01) THEN
            !
            ! Grandes ou petites valeurs absolues.
            ! Choix du format scientifique.
            !
            CLFOR='(E11.5)'
            KSOR=11
          ELSE
            !
            ! Valeurs absolues moyennes.
            ! Choix du format fixe.
            !
            CLFOR='(F11.4)'
            KSOR=11
          ENDIF
        ELSE
          !
          ! On affiche 7 chiffres significatifs.
          !
          ZVALA=ABS(PX)
          ZX=PX
          IF(ZVALA.EQ.0.) THEN
            !
            ! Nombre nul.
            ! Choix du format fixe.
            !
            CLFOR='(F15.8)'
            KSOR=15
          ELSEIF(ZVALA.GE.10000.) THEN
            !
            ! Grandes valeurs absolues.
            ! Choix du format scientifique.
            !
            CLFOR='(E15.7)'
            KSOR=15
          ELSEIF(ZVALA.GE.1.E-03) THEN
            !
            ! Valeurs absolues moyennes.
            ! Choix du format fixe.
            !
            CLFOR='(F15.9)'
            KSOR=15
          ELSE
            !
            ! Petites valeurs absolues.
            ! Choix du format scientifique.
            !
            CLFOR='(E15.7)'
            KSOR=15
          ENDIF
        ENDIF
        WRITE(CDSOR,FMT=CLFOR) ZX
      ELSE
        ! On veut un affichage a kopt chiffres apres la virgule.
        ZX=PX
        IF(ZX.NE.0.) THEN
          IF(ZX.GT.0.) THEN
            ISIGN=1
          ELSE
            ISIGN=-1
            ZX=-ZX
          ENDIF
          ILOG=NINT(LOG(ZX)/LOG(10.)-.5)
          ! Pour les nombres dont la mantisse est voisine de 1., il y a ris
          ! d'erreur d'arrondi. On teste si ce cas d'erreur est arrive,
          ! et si oui on modifie la valeur de ilog.
          ZXNOR=ZX/10.**REAL(ILOG)
          IF(ZXNOR.LT.1.) THEN
            ILOG=ILOG-1
          ELSEIF(ZXNOR.GE.10.) THEN
            ILOG=ILOG+1
          ENDIF
        ELSE
          ISIGN=1
          ILOG=0
        ENDIF
        ZINT=ZX*10.**REAL(KNC)
        ! Arrondi de zint a l'entier le plus proche.
        IF(ZINT.EQ.0.) THEN
          IX=0
        ELSE
          IX=NINT(ZINT-.5)
        ENDIF
        IY=NINT(10.**REAL(KNC))
        IF(MOD(IX,IY).EQ.0.AND.KNC.NE.0) THEN
          ! Le reel est voisin d'un entier
          IAFF=IX/NINT(10.**REAL(KNC))
          IF(IAFF.EQ.0) THEN
            CDSOR='0'
            KSOR=1
          ELSE
            ! Obtention du nombre de chiffres ilog2 de l'entier iaff
            ILOG2=1
            IAFFTMP=IAFF
            DO JCH=0,12
              IAFFTMP=IAFFTMP/10
              IF(IAFFTMP.NE.0) ILOG2=ILOG2+1
            ENDDO
            WRITE(CLFOR,FMT='(A2,I3.3,A1)')'(I',ILOG2,')'
            WRITE(CDSOR,FMT=CLFOR) IAFF
            KSOR=ILOG2
            IF(ISIGN.EQ.-1) THEN
              CDSOR='-'//CDSOR(1:KSOR)
              KSOR=KSOR+1
            ENDIF
          ENDIF
        ELSE
          ! Le reel n'est pas voisin d'un entier.
          ! On l'affiche au format reel + kopt ch. apres la virgule.
          IF(ILOG.LT.0) ILOG=0
          KSOR=ILOG+2+KOPT
          WRITE(CLFORF,FMT='(A2,I2.2,A1,I2.2,A1)') &
          '(F',KSOR,'.',KOPT,')'
          WRITE(CDSOR,FMT=CLFORF) ZX
          IF(CDSOR(1:1).EQ.' ') THEN
            CDSOR(1:1)='0'
          ENDIF
          IF(ISIGN.EQ.-1) THEN
            CDSOR='-'//CDSOR(1:KSOR)
            KSOR=KSOR+1
          ENDIF
        ENDIF
      ENDIF
 1000 FORMAT(F11.09)
      IF (LHOOK) CALL DR_HOOK('REECAR',1,ZHOOK_HANDLE)
      ENDSUBROUTINE REECAR
!
!
      FUNCTION CLLANG()
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
      ! --------------------------------------------------------------
      ! **** *CLLANG* Returns current language choice, depending on LANG environment variable.
      ! --------------------------------------------------------------
      ! Sujet:
      ! Arguments explicites:
      ! Arguments implicites:
      ! Methode:
      ! Externes:
      ! Auteur:   98-03, J.M. Piriou.
      ! Modifications:
      ! --------------------------------------------------------------
      ! En entree:
      ! En sortie:
      ! cllang='FRA', 'ENG', ...
      ! --------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*80 CLVAR
      CHARACTER*3 CLLANG
      !
      ! -------------------------------------------------
      ! Get $HOME current value.
      ! -------------------------------------------------
      !
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CLLANG',0,ZHOOK_HANDLE)
      CALL GET_ENVIRONMENT_VARIABLE('LANG',CLVAR)
      IF(CLVAR(1:2).EQ.'fr') THEN
        CLLANG='FRA'
      ELSE
        CLLANG='ENG'
      ENDIF
      IF (LHOOK) CALL DR_HOOK('CLLANG',1,ZHOOK_HANDLE)
      ENDFUNCTION CLLANG
