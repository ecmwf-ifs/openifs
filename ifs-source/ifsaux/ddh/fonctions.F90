SUBROUTINE CIN_CAPE(KLEV,PT,PQV,PP,PCIN,PCAPE,KLC,KLFC,KLNB,KCLOUD,PTHETAE_TOP,PPREC)
! --------------------------------------------------------------
! **** *cin_cape* Compute CIN, CAPE, and vertical location of clouds.
! --------------------------------------------------------------
! Subject:
! Method:
!   Ascents will be computed starting from each level of the input profile.
!   The parcel is raised from its level of origin LO to the level of condensation LC,
!   then to its level of free convection LFC (the parcel becomes buoyant),
!   then to the level of neutral buoyancy LNB, where the the parcel becomes unbuoyant.
!   All computations are done without any entrainment of environmental air.
! Externals:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! In input:
!   klev: number of levels.
!   pt: temperature (K).
!   pqv: specific humidity (no dim).
!   pp: pressure (Pa).
! WARNING: pressure values are supposed to icrease from pp(1) to pp(klev).
! In output:
!   pcin (J/kg): CIN, Convection INhibition,
!          massic energy to raise the parcel from
!          LO to LC then to LFC (see above).
!          Only negative terms are cumulated.
!   pcape (J/kg): CAPE, Convection Available Potential Energy,
!          massic energy  provided by the raise of the parcel
!          from LFC to LNB (see above).
!          Only positive terms are cumulated.
!   klc(jlev) LC of the parcel raised from jlev.
!   klfc(jlev) LFC of the parcel raised from jlev.
!   klnb(jlev) LNB of the parcel raised from jlev.
!   kcloud: 1 if convective cloud at this level,
!           2 if "near saturation" cloud at this level,
!           0 if no cloud at this level.
!           If both convective and "near saturation" are present, the ouput is 1.
!   pthetae_top (K): array receiving the potential temperature of the parcel raised up to the LNB.
!   pprec (kg/m2): precipitations cumulated over all ascents.
! --------------------------------------------------------------
!
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY: RG,RTT
IMPLICIT NONE
INTEGER(KIND=JPIM) :: IPOS
INTEGER(KIND=JPIM) :: JLEV
INTEGER(KIND=JPIM) :: JLEV1
LOGICAL :: LLDEBUG
LOGICAL :: LLBUOY
LOGICAL :: LLSAT
REAL(KIND=JPRB) :: ZDLOG
REAL(KIND=JPRB) :: ZBUOY
REAL(KIND=JPRB) :: ZBUOY_PREC
REAL(KIND=JPRB) :: ZPREC
REAL(KIND=JPRB) :: ZPRESS
REAL(KIND=JPRB) :: ZQS
REAL(KIND=JPRB) :: ZQV
REAL(KIND=JPRB) :: ZQV1,ZQV2
REAL(KIND=JPRB) :: ZRT
REAL(KIND=JPRB) :: ZT,ZRETAINED,ZADD
REAL(KIND=JPRB) :: ZT1,ZT2,ZQL,ZQI,ZQC,TV,THETAV,THETAVL,THETAL_KE,TRHO
INTEGER(KIND=JPIM) :: KLEV
REAL(KIND=JPRB) :: PPREC
REAL(KIND=JPRB) PT(KLEV)
REAL(KIND=JPRB) PQV(KLEV)
REAL(KIND=JPRB) PP(KLEV)
REAL(KIND=JPRB) PCIN(KLEV)
REAL(KIND=JPRB) PCAPE(KLEV)
INTEGER(KIND=JPIM) KLC(KLEV)
INTEGER(KIND=JPIM) KLFC(KLEV)
INTEGER(KIND=JPIM) KLNB(KLEV)
INTEGER(KIND=JPIM) KCLOUD(KLEV)
REAL(KIND=JPRB) PTHETAE_TOP(KLEV)
REAL(KIND=JPRB) QS,R_HUM,T_EVOL_AD_HUM,THETA,T_EVOL_AD_SECHE
!
!-------------------------------------------------
! Default initializations.
!-------------------------------------------------
!
LLDEBUG=.FALSE.
!
!-------------------------------------------------
! zretained: fraction of the condensates which is retained, i.e. which does not precipitate.
!       if zretained=1. ==> reversible moist ascent.
!                      it is assumed that all the parcel's condensed
!                      water is retained, thus liquid and ice sustents reduce the buoyancy.
!        if zretained=0. ==> "irreversible" (pseudo-adiabatic) moist ascent.
!                       liquid and ice sustents precipitate
!                       instantaneously and thus do not affect the buoyancy.
!       zretained can be used with values between 0. and 1..
!-------------------------------------------------
!
ZRETAINED=0._JPRB
!
!-------------------------------------------------
! Initialize to zero.
!-------------------------------------------------
!
PCIN=0._JPRB
PCAPE=0._JPRB
KLC=0
KLFC=0
KLNB=0
KCLOUD=0
PTHETAE_TOP=0._JPRB
PPREC=0._JPRB
!
!-------------------------------------------------
! Loop over origin of ascents.
! This loop can be done indifferently upwards or downwards.
!-------------------------------------------------
!
DO JLEV1=KLEV,1,-1
    !
    !-------------------------------------------------
    ! ipos:
    !   0 if parcel between LO and LC.
    !   1 if parcel between LC and LFC.
    !   2 if parcel between LFC and LNB.
    !   3 if parcel between LNB and top of profile.
    !-------------------------------------------------
    !
    IPOS=0
    !
    !-------------------------------------------------
    ! For each LO, one will raise the parcel up to the top.
    !-------------------------------------------------
    !
    ZT=PT(JLEV1)
    ZQV=PQV(JLEV1)
    IF(ZQV/QS(ZT,PP(JLEV1)) > 0.9) THEN
        !
        !-------------------------------------------------
        ! "Near saturation" cloud.
        !-------------------------------------------------
        !
        IF(KCLOUD(JLEV1) == 0) KCLOUD(JLEV1)=2
    ENDIF
    ZPREC=0._JPRB
    ZQL=0._JPRB
    ZQI=0._JPRB
    DO JLEV=JLEV1,2,-1
        !
        !-------------------------------------------------
        ! Saturation specific humidity.
        !-------------------------------------------------
        !
        ZQS=QS(ZT,PP(JLEV))
        !
        !-------------------------------------------------
        ! Pressure and saturation.
        !-------------------------------------------------
        !
        ZPRESS=PP(JLEV)
        LLSAT=ZQV > 0.99*ZQS ! true if saturated parcel.
        !
        !-------------------------------------------------
        ! The parcel buoyancy is computed from the ratio
        ! between density of the parcel (which is a mixture of dry air, water vapour zqv
        ! and condensates zqc) and density of the environmental air.
        ! Note that zqc is zql+zqi and thus will be 0. if zretained=0.,
        ! i.e. if liquid and ice sustents are supposed to precipitate instantaneously.
        !-------------------------------------------------
        !
        ZQC=ZQL+ZQI
        ZBUOY=(TRHO(ZT,ZQV,ZQC)/TV(PT(JLEV),PQV(JLEV))-1._JPRB)*R_HUM(PQV(JLEV))*PT(JLEV)
        LLBUOY=ZBUOY >= 0._JPRB ! true if buoyant parcel.
        !
        !-------------------------------------------------
        ! CIN and CAPE integrals.
        !-------------------------------------------------
        !
        IF(JLEV < JLEV1) THEN
            ZDLOG=LOG(PP(JLEV+1)/PP(JLEV))
            ZRT=0.5_JPRB*(ZBUOY+ZBUOY_PREC)*ZDLOG
            !zrt=0.5_JPRB*(zbuoy+zbuoy_prec)*(pp(jlev+1)-pp(jlev))/(0.5_JPRB*(pp(jlev)+pp(jlev+1)))
            IF(ZRT > 0._JPRB) THEN
                !
                !-------------------------------------------------
                ! Cumulate CAPE if positive contribution.
                !-------------------------------------------------
                !
                PCAPE(JLEV1)=PCAPE(JLEV1)+ZRT
            ELSE
                !
                !-------------------------------------------------
                ! Cumulate CIN if negative contribution and below LFC.
                !-------------------------------------------------
                !
                IF(IPOS <= 1) PCIN(JLEV1)=PCIN(JLEV1)+ZRT
            ENDIF
        ENDIF
        ZBUOY_PREC=ZBUOY
        IF(LLDEBUG .AND. JLEV1 == KLEV) THEN
            !
            !-------------------------------------------------
            ! Saturation and CAPE profiles.
            !-------------------------------------------------
            !
            WRITE(77,FMT='(3i3,5(a,g9.3),a,i3,2(a,l3),2(a,g9.3))') &
            & JLEV1,JLEV,KLEV-JLEV+1,' zqv=',ZQV,' zqs=',ZQS,' zt=',ZT-RTT &
            & ,' pt=',PT(JLEV)-RTT,' pp=',PP(JLEV)/100.,' ipos=',IPOS,' sat=',LLSAT,' buoy=',LLBUOY &
            &,' pcin=',PCIN(JLEV1),' pcape=',PCAPE(JLEV1)
            !
            !-------------------------------------------------
            ! Theta, thetaV and thetaVL profiles.
            !-------------------------------------------------
            !
            WRITE(78,FMT='(1(a,i2.2),5(a,g10.4),2(a,l3),9(a,g10.4))') &
            & 'niv=',JLEV,' theta = ',THETA(PP(JLEV),ZT),' thetaV = ',THETAV(PP(JLEV),ZT,ZQV) &
            & ,' thetaVL = ',THETAVL(PP(JLEV),ZT,ZQV,ZQC)
            !
            !-------------------------------------------------
            ! Theta, thetaV and thetaVL profiles.
            !-------------------------------------------------
            !
            WRITE(79,FMT='(9g16.7)') THETA(PP(JLEV),ZT),-REAL(JLEV)
            WRITE(80,FMT='(9g16.7)') THETAV(PP(JLEV),ZT,ZQV),-REAL(JLEV)
            WRITE(81,FMT='(9g16.7)') THETAL_KE(PP(JLEV),ZT,ZQV,ZQC),-REAL(JLEV)
            WRITE(82,FMT='(9g16.7)') THETAVL(PP(JLEV),ZT,ZQV,ZQC),-REAL(JLEV)
            !
            !-------------------------------------------------
            ! qv and qc=ql+qi profiles.
            !-------------------------------------------------
            !
            WRITE(83,FMT='(9g16.7)') ZQV,-REAL(JLEV)
            WRITE(84,FMT='(9g16.7)') ZQC,-REAL(JLEV)
        ENDIF
        !
        !-------------------------------------------------
        ! Check-up transitions between LC, LFC and LNB.
        !-------------------------------------------------
        !
        IF(LLSAT) THEN
            !
            !-------------------------------------------------
            ! Saturated parcel.
            !-------------------------------------------------
            !
            IF(LLBUOY) THEN
                !
                !-------------------------------------------------
                ! Buoyant parcel.
                !-------------------------------------------------
                !
                IF(IPOS == 0) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LC, one has found both LC and LFC!...
                    !-------------------------------------------------
                    !
                    IPOS=2
                    KLC(JLEV1)=JLEV
                    KLFC(JLEV1)=JLEV
                ELSEIF(IPOS == 1) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LFC, one has found LFC.
                    !-------------------------------------------------
                    !
                    IPOS=2
                    KLFC(JLEV1)=JLEV
                ELSEIF(IPOS == 2) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LNB, one has to go on raising.
                    !-------------------------------------------------
                    !
                ELSE
                    PRINT*,'cin_cape/ERROR: unexpected ipos!...'
                    STOP 'call abort'
                ENDIF
            ELSE
                !
                !-------------------------------------------------
                ! Unbuoyant parcel.
                !-------------------------------------------------
                !
                IF(IPOS == 0) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LC, one has found LC.
                    !-------------------------------------------------
                    !
                    IPOS=1
                    KLC(JLEV1)=JLEV
                ELSEIF(IPOS == 1) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LFC, one has to go on raising.
                    !-------------------------------------------------
                    !
                ELSEIF(IPOS == 2) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LNB, one has found LNB.
                    ! Raising from jlev1 can be stopped here.
                    !-------------------------------------------------
                    !
                    IPOS=3
                    KLNB(JLEV1)=JLEV
                    EXIT
                ELSE
                    PRINT*,'cin_cape/ERROR: unexpected ipos!...'
                    STOP 'call abort'
                ENDIF
            ENDIF
        ELSE
            !
            !-------------------------------------------------
            ! Unsaturated parcel.
            !-------------------------------------------------
            !
            IF(LLBUOY) THEN
                !
                !-------------------------------------------------
                ! Buoyant parcel.
                !-------------------------------------------------
                !
                IF(IPOS == 0) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LC, one has to go on raising.
                    !-------------------------------------------------
                    !
                ELSEIF(IPOS == 1) THEN
                    !
                    !-------------------------------------------------
                    ! The parcel is unsaturated and buoyant, above LC.
                    ! Thus another LC point exists above.
                    ! One restarts the search for LC.
                    !-------------------------------------------------
                    !
                    IPOS=0
                ELSEIF(IPOS == 2) THEN
                    !
                    !-------------------------------------------------
                    ! Go on raising to LNB (ipos=2).
                    !-------------------------------------------------
                    !
                ELSE
                    PRINT*,'cin_cape/ERROR: unexpected ipos!...'
                    STOP 'call abort'
                ENDIF
            ELSE
                !
                !-------------------------------------------------
                ! Unbuoyant parcel.
                !-------------------------------------------------
                !
                IF(IPOS == 0) THEN
                    !
                    !-------------------------------------------------
                    ! Go on raising to LC (ipos=0).
                    !-------------------------------------------------
                    !
                ELSEIF(IPOS == 1) THEN
                    !
                    !-------------------------------------------------
                    ! Go on raising to LFC (ipos=1).
                    !-------------------------------------------------
                    !
                ELSEIF(IPOS == 2) THEN
                    !
                    !-------------------------------------------------
                    ! While raising to LNB, one has found LNB.
                    ! Raising from jlev1 can be stopped here.
                    !-------------------------------------------------
                    !
                    IPOS=3
                    KLNB(JLEV1)=JLEV
                    EXIT
                ELSE
                    PRINT*,'cin_cape/ERROR: unexpected ipos!...'
                    STOP 'call abort'
                ENDIF
            ENDIF
        ENDIF
        IF(LLSAT) THEN
            !
            !-------------------------------------------------
            ! If the parcel is supersaturated, supersaturation is removed.
            ! This is done trough an isobaric transformation from (zt,zqv) to (zt1,zqv1).
            !-------------------------------------------------
            !
            CALL EVOL_AD_HUM_ISOBARE(ZT,ZQV,PP(JLEV),ZT1,ZQV1)
            ZADD=ZRETAINED*(ZQV-ZQV1)
            IF(ZT1 >= RTT) THEN
                ZQL=ZQL+ZADD
            ELSE
                ZQI=ZQI+ZADD
            ENDIF
            !
            !-------------------------------------------------
            ! Moist adiabatic ascent.
            ! Transformation from (zt1,zqv1) to (zt2,zqv2).
            !-------------------------------------------------
            !
            ZT2=T_EVOL_AD_HUM(ZT1,PP(JLEV),PP(JLEV-1))
            ZQV2=QS(ZT2,PP(JLEV-1))
            IF(ZQV1 > ZQV2) THEN
                ZADD=ZRETAINED*(ZQV1-ZQV2)
                IF(ZT2 >= RTT) THEN
                    ZQL=ZQL+ZADD
                ELSE
                    ZQI=ZQI+ZADD
                ENDIF
            ENDIF
            !
            !-------------------------------------------------
            ! Cumulate precipitations from ascent started at jlev1.
            !-------------------------------------------------
            !
            ZPREC=ZPREC+(ZQV-ZQV2)*(PP(JLEV)-PP(JLEV-1))/RG
            IF(ZQV < ZQV2) THEN
                PRINT*,'cin_cape/ERROR: negative precipitations!... ' &
                & ,ZPREC,ZQV,ZQV2,PP(JLEV),PP(JLEV-1)
                STOP 'call abort'
            ENDIF
            !
            !-------------------------------------------------
            ! Update parcel state.
            !-------------------------------------------------
            !
            ZT=ZT2
            ZQV=ZQV2
        ELSE
            !
            !-------------------------------------------------
            ! Dry adiabatic ascent.
            !-------------------------------------------------
            !
            ZT=T_EVOL_AD_SECHE(ZT,ZQV,PP(JLEV),PP(JLEV-1))
        ENDIF
    ENDDO
    IF(LLDEBUG .AND. JLEV1 == KLEV) PRINT*,'Final ql=',ZQL,', final qi=',ZQI
    !
    !-------------------------------------------------
    ! If CAPE is 0, CIN is put also to 0, in order
    ! not to saturate the graphics with very negatice CINs
    ! where anyway no convection is to be expected.
    !-------------------------------------------------
    !
    !if(pcape(jlev1) == 0._JPRB) pcin(jlev1) = 0._JPRB
    !
    !-------------------------------------------------
    ! If some CAPE exists, al levels from LC to LNB are cloudy ones.
    !-------------------------------------------------
    !
    IF(PCAPE(JLEV1) /= 0._JPRB) THEN
        IF(KLC(JLEV1) == 0 .OR. KLNB(JLEV1) == 0) THEN
            !
            !-------------------------------------------------
            ! Buoyancy in dry air. No cloud.
            !-------------------------------------------------
            !
        ELSE
            !
            !-------------------------------------------------
            ! Usual case, for which buoyancy appeared
            ! at LFC, i.e. AFTER saturation.
            !-------------------------------------------------
            !
            IF(LLDEBUG) PRINT*,'Parcel raised from level ',JLEV1,': LNB at ',KLNB(JLEV1),', LC at ',KLC(JLEV1)
            DO JLEV=KLNB(JLEV1),KLC(JLEV1)
                KCLOUD(JLEV)=1
            ENDDO
            PPREC=PPREC+ZPREC
        ENDIF
    ENDIF
    !
    !-------------------------------------------------
    ! Compute potential temperature at LNB.
    !-------------------------------------------------
    !
    IF(IPOS == 3) THEN
        PTHETAE_TOP(JLEV1)=THETA(ZPRESS,ZT)
    ELSE
        PTHETAE_TOP(JLEV1)=THETA(PP(JLEV1),PT(JLEV1))
    ENDIF
ENDDO
END SUBROUTINE CIN_CAPE
FUNCTION HCLA(KLEV,PTHETA,PZ,CDMETHODE)
! --------------------------------------------------------------
! **** *hcla* Hauteur de la couche limite atmosphérique.
! --------------------------------------------------------------
! Subject:
! Explicit arguments:
! Implicit arguments:
! Method:
! Externals:
! Auteur/author:   2002-03, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! Input:
!   klev: nombre de niveaux du profil de theta.
!   ptheta: profil vertical de theta (K).
!      Remarque: s'il est disponible, il est préférable de fournir en entrée
!      theta_vl plutôt que theta: le diagnostic de PBL sera plus précis.
!   pz: profil vertical de l'élévation (m).
!      ATTENTION: pz doit croître de klev à 1.
!      pz doit être une élévation au dessus de la surface,
!      i.e. ne doit pas être une altitude. pz(klev) peut être
!      nul ou non, i.e. les tableaux d'entrée
!      peuvent contenir ou non la surface.
!   cdmethode: méthode de calcul de la hauteur de CLA:
!      plusieurs méthodes ont en effet été codées dans la présente fonction.
!      Voir les commentaires plus bas pour détail des méthodes.
! Output:
!   hcla: hauteur de la couche limite atmosphérique (m).
! --------------------------------------------------------------
!
!-------------------------------------------------
! Types implicites.
!-------------------------------------------------
!
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
CHARACTER*(*) CDMETHODE
INTEGER(KIND=JPIM) KLEV,JLEV
REAL(KIND=JPRB) PTHETA(KLEV), PZ(KLEV)
REAL(KIND=JPRB) ZTHETA(KLEV+1), ZZ(KLEV+1)
REAL(KIND=JPRB) HCLA,ZKHI0,ZKHI1,ZINT,ZMOY,ZX,ZPSI
INTEGER(KIND=JPIM) ILEV
!
!-------------------------------------------------
! Si les tableaux d'entrée ne contiennent pas la surface, on l'y ajoute.
!-------------------------------------------------
!
IF(PZ(KLEV) /= 0._JPRB) THEN
    !
    !-------------------------------------------------
    ! Il faut ajouter la surface.
    !-------------------------------------------------
    !
    ILEV=KLEV+1
    DO JLEV=1,KLEV
        ZZ(JLEV)=PZ(JLEV)
        ZTHETA(JLEV)=PTHETA(JLEV)
    ENDDO
    ZZ(KLEV+1)=0._JPRB
    ZTHETA(KLEV+1)=PTHETA(KLEV)
ELSE
    !
    !-------------------------------------------------
    ! La surface est déjà présente. Rien à ajouter.
    !-------------------------------------------------
    !
    ILEV=KLEV
    ZZ=PZ
    ZTHETA=PTHETA
ENDIF
!
!-------------------------------------------------
! Test sur les méthodes.
!-------------------------------------------------
!
IF(CDMETHODE(1:LEN_TRIM(CDMETHODE)) == 'DINT') THEN
    !
    !-------------------------------------------------
    ! Méthode par double intégrale.
    ! Auteur: J.M. Piriou 2002.
    ! Avantages:
    !   - robustesse numérique
    ! Défauts:
    !-------------------------------------------------
    !
    ZKHI0=0.25
    ZKHI1=0.12
    !
    !-------------------------------------------------
    ! Intégrale ascendante.
    !-------------------------------------------------
    !
    HCLA=0._JPRB
    ZINT=0._JPRB
    DO JLEV=ILEV-1,1,-1
        !
        !-------------------------------------------------
        ! Intégrale de theta.
        !-------------------------------------------------
        !
        ZINT=ZINT+(ZZ(JLEV)-ZZ(JLEV+1))*0.5_JPRB*(ZTHETA(JLEV)+ZTHETA(JLEV+1))
        !
        !-------------------------------------------------
        ! Ecart entre theta du niveau courant
        ! et la moyenne de theta entre la surface et le niveau courant.
        !-------------------------------------------------
        !
        ZMOY=1._JPRB/ZZ(JLEV)*ZINT
        ZX=ZTHETA(JLEV)-ZMOY
        !
        !-------------------------------------------------
        ! Fonction poids psi.
        !-------------------------------------------------
        !
        !zpsi=1._JPRB/(1._JPRB+exp(max(-10.,min(10.,(zx-zkhi0)/(zkhi1/6.)))))
        ZPSI=MAX(0._JPRB,MIN(1._JPRB,(ZKHI0-ZX)/ZKHI1+0.5_JPRB))
        !
        !-------------------------------------------------
        ! Intégrale de dz, avec pour poids psi.
        !-------------------------------------------------
        !
        HCLA=HCLA+(ZZ(JLEV)-ZZ(JLEV+1))*ZPSI
    ENDDO
ELSEIF(CDMETHODE(1:LEN_TRIM(CDMETHODE)) == 'AY1996') THEN
    !
    !-------------------------------------------------
    ! Méthode par simple intégrale, ref Ayotte and al., 1996, in BLM.
    ! Avantages:
    !   - simplicité
    !   - robustesse numérique
    ! Défauts:
    !   - sensible aux "Diracs", i.e. aux niveaux où une forte
    !     fluctuation de theta existe, même si le niveau n'est pas d'épaisseur significative;
    !     exemples: présence d'instabilités numériques, forte inversion
    !     entre le niveau le plus bas et la surface,
    !     conduisant alors à une hauteur estimée de CLA nulle.
    !   - l'épaisseur de CLA est discrète, i.e. passe brusquement
    !     de la hauteur d'un niveau à celle du suivant, ce même si
    !     le profil de theta évolue continûment.
    !-------------------------------------------------
    !
    !
    !-------------------------------------------------
    ! Intégrale ascendante.
    !-------------------------------------------------
    !
    ZINT=0._JPRB
    DO JLEV=ILEV-1,1,-1
        !
        !-------------------------------------------------
        ! Intégrale de theta.
        !-------------------------------------------------
        !
        ZINT=ZINT+(ZZ(JLEV)-ZZ(JLEV+1))*0.5_JPRB*(ZTHETA(JLEV)+ZTHETA(JLEV+1))
        !
        !-------------------------------------------------
        ! Ecart entre theta du niveau courant
        ! et la moyenne de theta entre la surface et le niveau courant.
        !-------------------------------------------------
        !
        ZMOY=1._JPRB/ZZ(JLEV)*ZINT
        ZX=ZTHETA(JLEV)-ZMOY
        !
        !-------------------------------------------------
        ! .
        !-------------------------------------------------
        !
        IF(ZX > 0.25) THEN
            HCLA=ZZ(JLEV+1)
            RETURN
        ENDIF
    ENDDO
ELSE
    WRITE(*,FMT=*) 'hcla/ERREUR: méthode de calcul de la CLA non attendue!...'
    WRITE(*,FMT=*) CDMETHODE(1:LEN_TRIM(CDMETHODE))
    STOP 'call abort'
ENDIF
END FUNCTION HCLA
FUNCTION CP(PQV,PQL,PQI)
! --------------------------------------------------------------
! **** *cp* Chaleur massique de l'air humide.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! pqv  humidité spécifique vapeur (sans dimension).
! pql  humidité spécifique liquide (sans dimension).
! pqi  humidité spécifique glace (sans dimension).
! En sortie: chaleur massique de l'air humide (J/kg/K).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RCPV
USE CONST_THER, ONLY : RCW
USE CONST_THER, ONLY : RCS
IMPLICIT NONE
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PQL
REAL(KIND=JPRB) :: PQI
REAL(KIND=JPRB) CP
CP=RCPD*(1._JPRB-PQV-PQL-PQI)+RCPV*PQV+RCW*PQL+RCS*PQI
END FUNCTION CP
FUNCTION E(PQ,PP)
! --------------------------------------------------------------
! **** *e* Tension de vapeur en fonction de l'humidité spécifique et de la pression.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree: pq humidité spécifique (sans dimension).
! pp pression en Pa.
! En sortie: e en Pa.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RV

IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PQ
REAL(KIND=JPRB) E
E=PP*PQ/((1._JPRB-PQ)*RD/RV+PQ)
END FUNCTION E
FUNCTION ES(PT)
! --------------------------------------------------------------
! **** *es* Fonction es(T) par rapport à l'eau ou la glace.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.F. Geleyn.
! Modifications:
! --------------------------------------------------------------
! En entree: température en K.
! En sortie: es en Pa.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RALPD
USE CONST_THER, ONLY : RALPW
USE CONST_THER, ONLY : RBETD
USE CONST_THER, ONLY : RBETW
USE CONST_THER, ONLY : RGAMD
USE CONST_THER, ONLY : RGAMW
USE CONST_THER, ONLY : RTT
IMPLICIT NONE
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) ES
ES=EXP( &
& (RALPW+RALPD*MAX(0._JPRB,SIGN(1._JPRB,RTT-PT))) &
& -(RBETW+RBETD*MAX(0._JPRB,SIGN(1._JPRB,RTT-PT)))/PT &
& -(RGAMW+RGAMD*MAX(0._JPRB,SIGN(1._JPRB,RTT-PT)))*LOG(PT))
END FUNCTION ES
FUNCTION ESS(PT)
! --------------------------------------------------------------
! **** *ess* Fonction es(T) par rapport à la glace.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.F. Geleyn.
! Modifications:
! --------------------------------------------------------------
! En entree: température en K.
! En sortie: es en Pa.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RALPS
USE CONST_THER, ONLY : RBETS
USE CONST_THER, ONLY : RGAMS
IMPLICIT NONE
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) ESS
ESS=EXP(RALPS-RBETS/PT-RGAMS*LOG(PT))
END FUNCTION ESS
FUNCTION ESW(PT)
! --------------------------------------------------------------
! **** *esw* Fonction es(T) par rapport à l'eau.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.F. Geleyn.
! Modifications:
! --------------------------------------------------------------
! En entree: température en K.
! En sortie: es en Pa.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RALPW
USE CONST_THER, ONLY : RBETW
USE CONST_THER, ONLY : RGAMW
IMPLICIT NONE
REAL(KIND=JPRB) :: PT,ESW
ESW=EXP(RALPW-RBETW/PT-RGAMW*LOG(PT))
END FUNCTION ESW
SUBROUTINE EVOL_AD_HUM_ISOBARE(PT0,PQV0,PP0,PT,PQV)
! --------------------------------------------------------------
! **** *evol_ad_hum_isobare* Calcul de l'état final d'une condensation isobare.
! --------------------------------------------------------------
! Sujet:
! On passe d'un état (T0,qv0,p0) avec qv0 > qs(T0,p0)
! à un état (T,qv,p0), en vérifiant cp*dt+L*dq=0, et qv=qs(T,p0), et à pression constante.
! Arguments explicites:
! Arguments implicites:
! Methode:
!   On résout en T
!      f(T)=cp*(T-T0)+L*(q-q0)=0
!   avec la contrainte q=qs(T,p0),
!   On résout par la méthode de Newton, en itérant
!   T --> T-f(T)/f'(T), avec pour point de départ T0.
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pt0: température de départ (K).
!   pqv0: humidité spécifique de départ (kg/kg).
!   pp0: pression de départ et d'arrivée (Pa).
! En sortie:
!   pt: température d'arrivée (K).
!   pqv: humidité spécifique d'arrivée (kg/kg).
!       Elle est égale à qs(pt,pp0).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RCPV
USE CONST_THER, ONLY : RTT
IMPLICIT NONE
INTEGER(KIND=JPIM) :: JIT
REAL(KIND=JPRB) :: PP0
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PQV0
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: PT0
REAL(KIND=JPRB) :: ZCP
REAL(KIND=JPRB) :: ZDELARG
REAL(KIND=JPRB) :: ZL
REAL(KIND=JPRB) :: ZQI
REAL(KIND=JPRB) :: ZQL
REAL(KIND=JPRB) :: ZT_DEPART,QS,FOLH,CP,FDERQS,FDERFOLH

PT=PT0
PQV=QS(PT,PP0)
ZQL=0._JPRB
ZQI=0._JPRB
DO JIT=1,10
    ZT_DEPART=PT
    !
    !-------------------------------------------------
    ! Chaleur latente.
    !-------------------------------------------------
    !
    ZDELARG=MAX(0._JPRB,SIGN(1._JPRB,RTT-PT))
    ZL=FOLH(PT,ZDELARG)
    !
    !-------------------------------------------------
    ! Itération de Newton.
    !-------------------------------------------------
    !
    ZCP=CP(PQV,ZQL,ZQI)
    PT=PT-(ZCP*(PT-PT0)+ZL*(PQV-PQV0)) &
    & /(ZCP+(RCPV-RCPD+ZL)*FDERQS(PT,PP0)+PQV*FDERFOLH(ZDELARG))
    !
    !-------------------------------------------------
    ! Vapeur saturante.
    !-------------------------------------------------
    !
    PQV=QS(PT,PP0)
    !
    !-------------------------------------------------
    ! On sort de la boucle de Newton
    ! si on a la solution à epsilon près.
    !-------------------------------------------------
    !
    IF(ABS(PT-ZT_DEPART) < 0.01) EXIT
    ZT_DEPART=PT
ENDDO
END SUBROUTINE EVOL_AD_HUM_ISOBARE
FUNCTION FDERQS(PT,PP)
! --------------------------------------------------------------
! **** *fderqs* Fonction dérivée partielle par rapport à la température de l'humidité spécifique saturante.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:  2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree: pt température en K.
! pp pression en Pa.
! En sortie: humidité spécifique saturante (sans dimension).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RETV
USE CONST_THER, ONLY : RV
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: ZES
REAL(KIND=JPRB) :: ZRAPP,FDERQS,ES,FODLES

ZES=ES(PT)
ZRAPP=RV/RD*PP/ZES
FDERQS=ZRAPP/(ZRAPP-RETV)**2*FODLES(PT)
END FUNCTION FDERQS
FUNCTION FODLES(PT)
! --------------------------------------------------------------
! **** *fodles* Fonction d(ln(es(T)))/dT.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.F. Geleyn.
! Modifications:
! --------------------------------------------------------------
! En entree: température en K.
! En sortie: es en Pa.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RBETD
USE CONST_THER, ONLY : RBETW
USE CONST_THER, ONLY : RGAMD
USE CONST_THER, ONLY : RGAMW
USE CONST_THER, ONLY : RTT
IMPLICIT NONE
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: ZDELARG,FODLES

!
! FONCTION DERIVEE DU LOGARITHME NEPERIEN DE LA PRECEDENTE (FOEW) .
! INPUT : pt = TEMPERATURE
! PDELARG = 0 SI EAU (QUELQUE SOIT pt)
! 1 SI GLACE (QUELQUE SOIT pt).
ZDELARG=MAX(0._JPRB,SIGN(1._JPRB,RTT-PT))
FODLES = ( &
& ( RBETW+ZDELARG*RBETD ) &
& - ( RGAMW+ZDELARG*RGAMD ) * PT ) &
& / ( PT*PT )
END FUNCTION FODLES
FUNCTION FOLH(PT,PDELARG)
! --------------------------------------------------------------
! **** *folh* Fonction chaleur latente vapeur/eau ou vapeur/glace.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-09, J.F. Geleyn.
! Modifications:
! --------------------------------------------------------------
! FONCTION CHALEUR LATENTE .
! Entrée: pt = TEMPERATURE
! PDELARG = 0 SI EAU (QUELQUE SOIT pt)
! 1 SI GLACE (QUELQUE SOIT pt).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RBETD
USE CONST_THER, ONLY : RBETW
USE CONST_THER, ONLY : RGAMD
USE CONST_THER, ONLY : RGAMW
USE CONST_THER, ONLY : RV
USE CONST_THER, ONLY : RTT
IMPLICIT NONE
REAL(KIND=JPRB) :: PDELARG
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) FOLH

FOLH=RV*((RBETW+PDELARG*RBETD)-(RGAMW+PDELARG*RGAMD)*PT) ! L complet (eau liquide/glace).
!folh=(2.501-0.00237*(pt-rtt))*1.e6 ! L Bolton (i.e. L vapeur/eau liquide).
!folh=2.5e6 ! L constant.

END FUNCTION FOLH
FUNCTION FDERFOLH(PDELARG)
! --------------------------------------------------------------
! **** *fderfolh* Fonction dérivée partielle par rapport à la température de la fonction chaleur latente vapeur/eau ou vapeur/glace.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:  2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! FONCTION CHALEUR LATENTE .
! Entrée:
! PDELARG = 0 SI EAU (QUELQUE SOIT pt)
! 1 SI GLACE (QUELQUE SOIT pt).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RGAMD
USE CONST_THER, ONLY : RGAMW
USE CONST_THER, ONLY : RV
IMPLICIT NONE
REAL(KIND=JPRB) :: PDELARG
REAL(KIND=JPRB) FDERFOLH

FDERFOLH=-RV*(RGAMW+PDELARG*RGAMD) ! L complet (eau liquide/glace).
!fderfolh=-0.00237e6 ! L Bolton (i.e. L vapeur/eau liquide).
!fderfolh=0._JPRB ! L constant.

END FUNCTION FDERFOLH
FUNCTION HR(PP,PT,PQV)
! --------------------------------------------------------------
! **** *hr* Fonction humidité relative.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-10, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
! En sortie:
!   hr (sans dimension).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PT,HR,E,ES

HR=E(MAX(0._JPRB,PQV),PP)/ES(PT)
END FUNCTION HR
SUBROUTINE POINT_CONDENS(PT0,PQV0,PP0,PTCOND,PPCOND)
! --------------------------------------------------------------
! **** *point_condens* Calcul du point de condensation d'une particule donnée par (T0, qv0, p0).
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
!   On résout en p
!      qs(T,p)=qv0
!      avec T=T0*(p0/p)**(R/cp)
! On résout par la méthode de Newton, en itérant
!   p --> p-f(p)/f'(p), avec pour point de départ p0.
! Externes:
! Auteur/author:   2000-10, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pt0: température de départ (K).
!   pqv0: humidité spécifique de départ (Pa).
!   pp0: pression de départ (Pa).
! En sortie:
!   ptcond température du point de condensation (K).
!   ppcond pression du point de condensation (Pa).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM) :: JIT
REAL(KIND=JPRB) :: PP0
REAL(KIND=JPRB) :: PPCOND
REAL(KIND=JPRB) :: PQV0
REAL(KIND=JPRB) :: PT0
REAL(KIND=JPRB) :: PTCOND
REAL(KIND=JPRB) :: ZDERI
REAL(KIND=JPRB) :: ZDP
REAL(KIND=JPRB) :: ZFONC
REAL(KIND=JPRB) :: ZPPREC
REAL(KIND=JPRB) :: ZQV
REAL(KIND=JPRB) :: ZTPLUS,QS,T_EVOL_AD_SECHE

PTCOND=PT0
PPCOND=PP0
ZDP=5.
ZQV=0._JPRB
DO JIT=1,10
    ZPPREC=PPCOND
    ZFONC=QS(PTCOND,PPCOND) ! valeur en p.
    ZTPLUS=T_EVOL_AD_SECHE(PTCOND,ZQV,PPCOND,PPCOND+ZDP) ! T(p+dp).
    ZDERI=(QS(ZTPLUS,PPCOND+ZDP)-ZFONC)/ZDP ! dérivée [qs(T(p+dp),p+dp)-qs(T,p)]/dp.
    PPCOND=PPCOND-(ZFONC-PQV0)/ZDERI
    PTCOND=T_EVOL_AD_SECHE(PTCOND,ZQV,ZPPREC,PPCOND)
    IF(ABS(PPCOND-ZPPREC) < 50.) EXIT
ENDDO
END SUBROUTINE POINT_CONDENS
FUNCTION QS(PT,PP)
! --------------------------------------------------------------
! **** *qs* humidité spécifique saturante en fonction de T et p.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree: pt température en K.
! pp pression en Pa.
! En sortie: humidité spécifique saturante (sans dimension).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: ZES,QS,ES,QV

ZES=ES(PT)
QS=QV(ZES,PP)
END FUNCTION QS
FUNCTION QV(PE,PP)
! --------------------------------------------------------------
! **** *qv* qv en fonction de la tension de vapeur et de la pression.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree: pe tension de vapeur en Pa.
! pp pression en Pa.
! En sortie: q (sans dimension).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RV
IMPLICIT NONE
REAL(KIND=JPRB) :: PE
REAL(KIND=JPRB) :: PP,QV

QV=PE/(RV/RD*(PP-PE)+PE)
END FUNCTION QV
FUNCTION R_HUM(PQV)
! --------------------------------------------------------------
! **** *r_hum* Constante spécifique de l'air humide.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-10, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
! pqv  humidité spécifique vapeur (sans dimension).
! En sortie: constante spécifique de l'air humide (J/kg/K).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RV
IMPLICIT NONE
REAL(KIND=JPRB) :: PQV,R_HUM

R_HUM=RD+(RV-RD)*PQV
END FUNCTION R_HUM
FUNCTION T_EVOL_AD_HUM(PT0,PP0,PP)
! --------------------------------------------------------------
! **** *t_evol_ad_hum* Calcul de l'état final d'une évolution pseudo-adiabatique humide.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
!   On résout en T
!      f(T)=cp*(T-T0)+L*(q-q0)+phi-phi0=0
!   soit f(T)=cp*(T-T0)+L*(q-q0)-R*T*log(p/p0)=0
!   avec la contrainte q=qs(T,p),
!   et sachant que q0=qs(T0,p0).
!   On résout par la méthode de Newton, en itérant
!   T --> T-f(T)/f'(T), avec pour point de départ T0.
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pt0: température de départ (K).
!   pp0: pression de départ (Pa).
!   pp: pression d'arrivée (Pa).
! En sortie:
!   température d'arrivée (K).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RCPV
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RTT
USE CONST_THER, ONLY : RV
IMPLICIT NONE
INTEGER(KIND=JPIM) :: IETAPES
INTEGER(KIND=JPIM) :: JETAPES
INTEGER(KIND=JPIM) :: JIT
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PP0
REAL(KIND=JPRB) :: PT0
REAL(KIND=JPRB) :: ZCP
REAL(KIND=JPRB) :: ZDELARG
REAL(KIND=JPRB) :: ZDERL
REAL(KIND=JPRB) :: ZETAPE
REAL(KIND=JPRB) :: ZL
REAL(KIND=JPRB) :: ZLOG
REAL(KIND=JPRB) :: ZP_ARRIVEE
REAL(KIND=JPRB) :: ZP_DEPART
REAL(KIND=JPRB) :: ZQI
REAL(KIND=JPRB) :: ZQL
REAL(KIND=JPRB) :: ZQV
REAL(KIND=JPRB) :: ZQV_DEPART
REAL(KIND=JPRB) :: ZRDLOG
REAL(KIND=JPRB) :: ZT
REAL(KIND=JPRB) :: ZT_DEPART,T_EVOL_AD_HUM,QS,FOLH,FDERFOLH,R_HUM,CP,FDERQS
REAL(KIND=JPRB) :: ZT_PREC

!
!-------------------------------------------------
! L'écart de pression à effectuer doit être
! cassé en suffisamment d'étapes pour que le calcul
! discret soit suffisamment précis.
!-------------------------------------------------
!
ZETAPE=5000. ! pas de pression conseillé (Pa).
IETAPES=NINT(ABS(PP-PP0)/ZETAPE)+1
ZT=PT0
DO JETAPES=1,IETAPES
    ZP_DEPART=PP0+(PP-PP0)*REAL(JETAPES-1)/REAL(IETAPES)
    ZP_ARRIVEE=PP0+(PP-PP0)*REAL(JETAPES)/REAL(IETAPES)
    ZT_DEPART=ZT
    ZQV_DEPART=QS(ZT_DEPART,ZP_DEPART)
    ZQL=0._JPRB
    ZQI=0._JPRB
    ZLOG=LOG(ZP_ARRIVEE/ZP_DEPART)
    DO JIT=1,10
        !
        !-------------------------------------------------
        ! Vapeur saturante.
        !-------------------------------------------------
        !
        ZQV=QS(ZT,ZP_ARRIVEE)
        !
        !-------------------------------------------------
        ! Chaleur latente.
        !-------------------------------------------------
        !
        ZDELARG=MAX(0._JPRB,SIGN(1._JPRB,RTT-ZT_DEPART))
        !
        ! Calcul où l'on considère que toute la chaleur latente libérée
        ! est récupérée par la vapeur d'eau pour se réchauffer.
        !
        ZL=FOLH(ZT,ZDELARG) ; ZDERL=FDERFOLH(ZDELARG)
;        !
        ! Calcul où l'on considère qu'une part seulement de la chaleur latente libérée
        ! est récupérée par la vapeur d'eau pour se réchauffer:
        ! le reste reste dans l'eau liquide ou glace.
        ! Ce calcul fournit des CAPE plus réalistes,
        ! mais n'est pas conservatif d'une énergie totale
        ! qui ne prendrait pas en compte l'énergie interne des particules
        ! d'eau et de glace.
        !
        !zratio=0.92 ; zl=zratio*folh(zt,zdelarg) ; zderl=zratio*fderfolh(zdelarg)
        !
        !-------------------------------------------------
        ! Itération de Newton.
        !-------------------------------------------------
        !
        ZRDLOG=R_HUM(ZQV)*ZLOG
        ZCP=CP(ZQV,ZQL,ZQI)
        ZT_PREC=ZT
        ZT=ZT-(ZCP*(ZT-ZT_DEPART)+ZL*(ZQV-ZQV_DEPART)-ZT*ZRDLOG) &
        & /(ZCP+(RCPV-RCPD+ZL-ZLOG*ZT*(RV-RD))*FDERQS(ZT,ZP_ARRIVEE)+ZQV*ZDERL-ZRDLOG)
        !
        !-------------------------------------------------
        ! On sort de la boucle de Newton
        ! si on a la solution à epsilon près.
        !-------------------------------------------------
        !
        IF(ABS(ZT-ZT_PREC) < 0.01) EXIT
    ENDDO
ENDDO
T_EVOL_AD_HUM=ZT
END FUNCTION T_EVOL_AD_HUM
FUNCTION T_EVOL_AD_SECHE(PT0,PQV0,PP0,PP)
! --------------------------------------------------------------
! **** *t_evol_ad_seche* Calcul de l'état final d'une évolution adiabatique sèche.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pt0: température de départ (K).
!   pqv0: humidité spécifique de départ (et d'arrivée, puisqu'aucune condensation ici) (kg/kg).
!         l'humidité spécifique sert simplement à calculer R et cp.
!   pp0: pression de départ (Pa).
!   pp: pression d'arrivée (Pa).
! En sortie:
!   t_evol_ad_seche: température d'arrivée (K).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PP0
REAL(KIND=JPRB) :: PQV0
REAL(KIND=JPRB) :: PT0
REAL(KIND=JPRB) :: ZQI
REAL(KIND=JPRB) :: ZQL,T_EVOL_AD_SECHE,R_HUM,CP

ZQL=0._JPRB
ZQI=0._JPRB
T_EVOL_AD_SECHE=PT0*(PP/PP0)**(R_HUM(PQV0)/CP(PQV0,ZQL,ZQI))
END FUNCTION T_EVOL_AD_SECHE
FUNCTION TD(PE)
! --------------------------------------------------------------
! **** *td* Calcul de Td en fonction de la tension de vapeur.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   96-04, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree: pe tension de vapeur en Pa.
! En sortie: Td en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM) :: JIT
REAL(KIND=JPRB) :: PE
REAL(KIND=JPRB) :: ZDERI
REAL(KIND=JPRB) :: ZDT
REAL(KIND=JPRB) :: ZES
REAL(KIND=JPRB) :: ZT
REAL(KIND=JPRB) :: ZTPREC,ES,TD

ZT=280.
!
! On itère une boucle de Newton
! pour annuler la fonction es(t)-e.
!
ZDT=0.1
DO JIT=1,10
    ZES=ES(ZT)
    ZDERI=(ES(ZT+ZDT)-ZES)/ZDT
    ZTPREC=ZT
    ZT=ZT-(ZES-PE)/ZDERI
    IF(ABS(ZT-ZTPREC) < 0.01) EXIT
ENDDO
TD=ZT
END FUNCTION TD
FUNCTION THETA(PP,PT)
! --------------------------------------------------------------
! **** *theta* Fonction température potentielle.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   98-01, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pt température en K.
!   pp pression en Pa.
! En sortie:
!   theta en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RD
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) THETA

THETA=PT*(RATM/PP)**(RD/RCPD)
END FUNCTION THETA
FUNCTION T(PP,PTHETA)
! --------------------------------------------------------------
! **** *t* Compute temperature from pressure and potential temperature.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2001-01, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   ptheta potential temperature in K.
!   pp pressure in Pa.
! En sortie:
!   t in K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RD
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PTHETA
REAL(KIND=JPRB) T

T=PTHETA*(PP/RATM)**(RD/RCPD)
END FUNCTION T
FUNCTION THETAD(PP,PQV)
! --------------------------------------------------------------
! **** *thetad* Fonction température de rosée potentielle.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
! En sortie:
!   thetad en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RD
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PQV,THETAD,E,TD

THETAD=TD(E(PQV,PP))*(RATM/PP)**(RD/RCPD)
END FUNCTION THETAD
FUNCTION THETAE(PP,PT,PQV)
! --------------------------------------------------------------
! **** *thetae* Fonction température potentielle équivalente (calcul discret précis).
! --------------------------------------------------------------
! Sujet:
! Le thetae d'une particule est la température qu'elle aurait
! si on la montait selon une adiabatique sèche jusqu'en son point
! de condensation, puis selon une adiabatique humide jusqu'à
! épuiser son humidité spécifique, puis on la redescendait
! selon une adiabatique sèche jusqu'au niveau de pression standard (ratm dans le code).
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-10, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
! En sortie:
!   thetae en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
INTEGER(KIND=JPIM) :: IETAPES
INTEGER(KIND=JPIM) :: JETAPES
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: ZDELTAP
REAL(KIND=JPRB) :: ZP
REAL(KIND=JPRB) :: ZPLIM
REAL(KIND=JPRB) :: ZPPREC
REAL(KIND=JPRB) :: ZQV
REAL(KIND=JPRB) :: ZQV1
REAL(KIND=JPRB) :: ZQVPREC
REAL(KIND=JPRB) :: ZT
REAL(KIND=JPRB) :: ZT1
REAL(KIND=JPRB) :: ZTPREC,THETAE,QS,T_EVOL_AD_HUM,T_EVOL_AD_SECHE,THETA

LOGICAL LLSAT
!
!-------------------------------------------------
! Ascendance jusqu'à épuisement de l'humidité spécifique,
! ce qu'on suppose arriver au niveau p=zplim
!-------------------------------------------------
!
LLSAT=PQV >= QS(PT,PP)
IF(LLSAT) THEN
    !
    !-------------------------------------------------
    ! La particule est sursaturée dès le départ.
    ! On élimine cette sursaturation.
    !-------------------------------------------------
    !
    CALL EVOL_AD_HUM_ISOBARE(PT,PQV,PP,ZT,ZQV)
ELSE
    ZT=PT
    ZQV=PQV
ENDIF
ZPLIM=10000. ! pression d'arrêt d'ascendance (Pa).
ZDELTAP=500. ! pas de pression pour la discrétisation (Pa).
IETAPES=NINT(ABS(PP-ZPLIM)/ZDELTAP)+1
ZPPREC=PP
ZTPREC=ZT
ZQVPREC=ZQV
DO JETAPES=1,IETAPES
    ZP=PP+(ZPLIM-PP)*REAL(JETAPES)/REAL(IETAPES)
    IF(LLSAT) THEN
        !
        !-------------------------------------------------
        ! Particule saturée. Ascendance adiabatique humide.
        !-------------------------------------------------
        !
        ZT=T_EVOL_AD_HUM(ZTPREC,ZPPREC,ZP)
        ZQV=QS(ZT,ZP)
    ELSE
        !
        !-------------------------------------------------
        ! Particule insaturée. Ascendance adiabatique sèche.
        !-------------------------------------------------
        !
        ZT=T_EVOL_AD_SECHE(ZTPREC,ZQVPREC,ZPPREC,ZP)
        IF(ZQVPREC > QS(ZT,ZP)) THEN
            !
            !-------------------------------------------------
            ! On atteint le point de condensation.
            ! On élimine cette sursaturation.
            !-------------------------------------------------
            !
            CALL EVOL_AD_HUM_ISOBARE(ZT,ZQVPREC,ZP,ZT1,ZQV1)
            LLSAT=.TRUE.
            ZT=ZT1
            ZQV=ZQV1
        ELSE
            !
            !-------------------------------------------------
            ! On n'atteint pas le point de condensation.
            ! qv est reconduit égal à lui-même.
            !-------------------------------------------------
            !
            ZQV=ZQVPREC
        ENDIF
    ENDIF
    !
    !-------------------------------------------------
    ! Le niveau courant devient le précédent.
    !-------------------------------------------------
    !
    ZTPREC=ZT
    ZQVPREC=ZQV
    ZPPREC=ZP
ENDDO
!
!-------------------------------------------------
! thetae n'est autre que le theta de la particule
! parvenue au sommet.
!-------------------------------------------------
!
THETAE=THETA(ZPPREC,ZT)
END FUNCTION THETAE
FUNCTION THETAE_BOLTON(PP,PT,PQV)
! --------------------------------------------------------------
! **** *thetae_bolton* Fonction température potentielle équivalente (calcul direct approché).
! --------------------------------------------------------------
! Sujet:
! Le thetae d'une particule est la température qu'elle aurait
! si on la montait selon une adiabatique sèche jusqu'en son point
! de condensation, puis selon une adiabatique humide jusqu'à
! épuiser son humidité spécifique, puis on la redescendait
! selon une adiabatique sèche jusqu'au niveau de pression standard.
! Arguments explicites:
! Arguments implicites:
! Methode: David Bolton, MWR 1980.
! La formule se veut rapide à calculer, et donc
! fait des hypothèses, telle une dépendance affine de L en T
! sur toute la plage de températures, etc...
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
! En sortie:
!   thetae_bolton en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: ZE
REAL(KIND=JPRB) :: ZR
REAL(KIND=JPRB) :: ZTCOND,THETAE_BOLTON,E

ZE=E(MAX(0._JPRB,PQV),PP)
ZR=PQV/(1._JPRB-PQV)*1000. ! r en g/kg.
IF(ZE /= 0._JPRB) THEN
    ZTCOND=2840./(3.5*LOG(PT)-LOG(ZE/100.)-4.805)+55. ! equ. (21) de Bolton 1980.
    !zhr=ze/es(pt) ; ztcond=1._JPRB/(1._JPRB/(pt-55.)-log(zhr)/2840.)+55. ! equ. (22) de Bolton 1980.
ELSE
    ZTCOND=PT
ENDIF
THETAE_BOLTON=PT*(RATM/PP)**(0.2854*(1._JPRB-0.28E-3*ZR)) &
&*EXP((3.376/ZTCOND-0.00254)*ZR*(1+0.81E-3*ZR)) ! equ. (43) de Bolton 1980.
END FUNCTION THETAE_BOLTON
FUNCTION THETAES(PP,PT)
! --------------------------------------------------------------
! **** *thetae* Fonction température potentielle équivalente de saturation (calcul discret précis).
! --------------------------------------------------------------
! Sujet:
! Le thetae d'une particule est la température qu'elle aurait
! si elle était saturée au départ (qv=qs(T,p)),
! si on la montait selon une adiabatique sèche jusqu'en son point
! de condensation, puis selon une adiabatique humide jusqu'à
! épuiser son humidité spécifique, puis on la redescendait
! selon une adiabatique sèche jusqu'au niveau de pression standard (ratm dans le code).
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-10, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
! En sortie:
!   thetaes en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: THETAES,THETAE,QS

THETAES=THETAE(PP,PT,QS(PT,PP))
END FUNCTION THETAES
FUNCTION THETAES_BOLTON(PP,PT)
! --------------------------------------------------------------
! **** *thetaes_bolton* Fonction température potentielle équivalente de saturation (calcul direct approché).
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
! En sortie:
!   thetaes_bolton en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PT,THETAE_BOLTON,QS,THETAES_BOLTON
THETAES_BOLTON=THETAE_BOLTON(PP,PT,QS(PT,PP))
END FUNCTION THETAES_BOLTON
FUNCTION THETAPW(PP,PT,PQV,PQL,PQI)
! --------------------------------------------------------------
! **** *thetapw* Fonction theta'w.
! --------------------------------------------------------------
! Sujet:
! CALCUL DE LA TEMPERATURE POTENTIELLE PSEUDO-ADIABATIQUE
! DU THERMOMETRE MOUILLE.
! Arguments explicites:
! Arguments implicites:
!
! Methode.
! --------
!
! CETTE ROUTINE N'A QU'UNE DIMENSION POUR SES VARIABLES D'ENTREE
! AFIN D'ETRE LA PLUS GENERALE POSSIBLE (UTILISATION SUR LES NIVEAUX
! DU MODELE COMME DU POST-PROCESSING, PAR EXEMPLE). TOUT ETAT DE
! L'AIR REALISTE EST ADMIS EN ENTREE ET L'ALGORITHME PREND EN COMPTE
! AUTOMATIQUEMENT UNE POSSIBLE TRANSITION DE PHASE LIQUIDE/GLACE.
! TROIS EQUATIONS IMPLICITES SONT RESOLUES PAR METHODE DE NEWTON:
! - RECHERCHE DU POINT DE SATURATION D'ENTROPIE EGALE PAR
! TRANSFORMATION REVERSIBLE ;
! - RECHERCHE DU POINT DE TEMPERATURE EGALE A CELLE DU POINT
! TRIPLE LE LONG DE L'ADIABATIQUE SATUREE IRREVERSIBLE ;
! - RECHERCHE DU POINT DE PRESSION EGALE A LA REFERENCE
! ATMOSPHERIQUE LE LONG D'UNE AUTRE (PARFOIS LA MEME) ADIABATIQUE
! IRREVERSIBLE.
! REMARQUES :
! - POUR LA PREMIERE ETAPE LA FORME SYMETRIQUE DE L'ENTROPIE
! HUMIDE PROPOSEE PAR P. MARQUET EST UTILISEE AFIN DE PERMETTRE UN
! MELANGE DE PHASES LIQUIDE ET GLACE DANS L'ETAT DE L'AIR ;
! - POUR LES DEUX DERNIERES ETAPES, PLUTOT QUE DE NEGLIGER
! COMME DE COUTUME LE TERME CONTENANT LE CP DU CONDENSAT, L'AUTEUR
! DE LA ROUTINE EN A DERIVE UNE APPROXIMATION QUASI-EXACTE ET PLUTOT
! BON MARCHE ;
! - POUR CES DEUX MEMES ETAPES, LES EBAUCHES DES BOUCLES DE
! NEWTON SONT OBTENUES PAR EXTRAPOLATION D'UNE LINEARISATION LOCALE
! APPROCHEE DES EQUATIONS ADIABATIQUE SATUREES.
!
! THIS ROUTINE HAS ONLY ONE DIMENSIONAL INPUT/OUTPUT ARRAYS IN
! ORDER TO BE THE MOST GENERAL POSSIBLE (USE ON MODEL OR ON POST-
! PROCESSING LEVELS, FOR EXAMPLE). ALL POSSIBLE REALISTIC INPUT
! STATES ARE ALLOWED AND THE ALGORITHM AUTOMATICALLY TAKES INTO
! ACCOUNT THE POTENTIAL LIQUID/ICE WATER TRANSITION.
! THREE IMPLICIT EQUATIONS ARE SOLVED BY NEWTON METHODS :
! - SEARCH OF THE SATURATION POINT OF EQUAL ENTROPY UNDER A
! REVERSIBLE TRANSFORM ;
! - SEARCH OF THE POINT OF TEMPERATURE EQUAL TO THAT OF THE
! TRIPLE POINT ALONG THE IRREVERSIBLE MOIST ADIABAT ;
! - SEARCH OF THE POINT OF REFERENCE ATMOSPHERIC PRESSURE
! ALONG ANOTHER (SOMETIMES IDENTICAL) IRREVERSIBLE MOIST ADIABAT.
! REMARKS :
! - FOR THE FIRST STEP THE SYMETRIC FORM OF THE MOIST ENTROPY
! PROPOSED BY P. MARQUET IS USED IN ORDER TO ALLOW A MIX OF LIQUID
! AND ICE WATER IN THE ATMOSPHERIC STATE ;
! - FOR THE TWO LAST STEPS, RATHER THAN THE USUAL NEGLECTION
! OF THE TERM MULTIPLIED BY CP OF THE CONDENSATE, THE ROUTINE'S
! AUTHOR DERIVED A QUASI EXACT AND NOT TOO EXPENSIVE ANALYTICAL
! APPROXIMATION FOR IT ;
! - FOR THE SAME STEPS, THE GUESSES OF THE NEWTON LOOP ARE
! OBTAINED BY VERTICAL EXTRAPOLATION OF A LINEAR LOCAL APPROXIMATION
! OF THE MOIST ADIABATS.
!
! Auteur/author: 92-09, J.F. Geleyn.
!
! Modifications.
! --------------
! 96-04, J. Calvo: Introduced a minimun in RH instead of a mini-
! mun in PQV. Added a security threshold in the
! calculation of the triple  point pressure
! first guess.
! --------------------------------------------------------------
! En entree:
! En sortie:
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RCPV
USE CONST_THER, ONLY : RCS
USE CONST_THER, ONLY : RCW
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RESTT
USE CONST_THER, ONLY : RKAPPA
USE CONST_THER, ONLY : RTT
USE CONST_THER, ONLY : RV, NBITER
IMPLICIT NONE
INTEGER(KIND=JPIM) :: JIT
REAL(KIND=JPRB) :: ZCPLNTT
REAL(KIND=JPRB) :: ZCWS
REAL(KIND=JPRB) :: ZDELTA
REAL(KIND=JPRB) :: ZDF
REAL(KIND=JPRB) :: ZDLEW
REAL(KIND=JPRB) :: ZDLSTT
REAL(KIND=JPRB) :: ZE
REAL(KIND=JPRB) :: ZEPSP
REAL(KIND=JPRB) :: ZEPSRH
REAL(KIND=JPRB) :: ZEW
REAL(KIND=JPRB) :: ZF
REAL(KIND=JPRB) :: ZKAPPA
REAL(KIND=JPRB) :: ZKDI
REAL(KIND=JPRB) :: ZKDW
REAL(KIND=JPRB) :: ZKNI
REAL(KIND=JPRB) :: ZKNW
REAL(KIND=JPRB) :: ZLH
REAL(KIND=JPRB) :: ZLHZ
REAL(KIND=JPRB) :: ZLHZI
REAL(KIND=JPRB) :: ZLHZW
REAL(KIND=JPRB) :: ZLSTC
REAL(KIND=JPRB) :: ZLSTT
REAL(KIND=JPRB) :: ZLSTTI
REAL(KIND=JPRB) :: ZLSTTW
REAL(KIND=JPRB) :: ZPRS
REAL(KIND=JPRB) :: ZQI
REAL(KIND=JPRB) :: ZQL
REAL(KIND=JPRB) :: ZQV
REAL(KIND=JPRB) :: ZQVSAT
REAL(KIND=JPRB) :: ZRDSV
REAL(KIND=JPRB) :: ZRI
REAL(KIND=JPRB) :: ZRL
REAL(KIND=JPRB) :: ZRS
REAL(KIND=JPRB) :: ZRV
REAL(KIND=JPRB) :: ZTIN
REAL(KIND=JPRB) :: ZTMAX
REAL(KIND=JPRB) :: ZTMIN,FOLH,QS,ES,FODLES,THETAPW

REAL(KIND=JPRB) PP,PT,PQV,PQL,PQI
REAL(KIND=JPRB) ZRT,ZCPT,ZRRT,ZS1,ZCONS,ZTITER,ZDEL,ZS2,ZFUNCT,ZPITER
!
! *
! ------------------------------------------------------------------
! I - CALCUL DES const_ther DERIVEES DE CELLES ISSUES DE YOMCST ET
! const_ther DE SECURITE.
!
! COMPUTATION OF CONSTANTS DERIVED FROM THOSE CONTAINED IN
! YOMCST AND SECURITY CONSTANTS.
!
ZRDSV=RD/RV
ZLSTTW= FOLH (RTT,0._JPRB)/RTT+RCW*RV/( FOLH (RTT,0._JPRB)/RTT)
ZLSTTI= FOLH (RTT,1._JPRB)/RTT+RCS*RV/( FOLH (RTT,1._JPRB)/RTT)
ZLHZW= FOLH (0._JPRB,0._JPRB)
ZLHZI= FOLH (0._JPRB,1._JPRB)
ZKNW=RESTT* FOLH (RTT,0._JPRB)/(RV*RTT)
ZKNI=RESTT* FOLH (RTT,1._JPRB)/(RV*RTT)
ZKDW=RKAPPA*RESTT*( FOLH (RTT,0._JPRB)/(RV*RTT))**2
ZKDI=RKAPPA*RESTT*( FOLH (RTT,1._JPRB)/(RV*RTT))**2
ZCPLNTT=RCPD*LOG(RTT)
!
ZEPSP=10.
ZEPSRH=0.001
ZTMIN=155.
ZTMAX=355.
!
!
! *
! ------------------------------------------------------------------
! II - CALCUL DE LA TEMPERATURE DE SATURATION (EN GENERAL MAIS PAS
! FORCEMENT TEMPERATURE DU POINT DE CONDENSATION). LE RAPPORT DE
! MELANGE TOTAL ET LA TEMPERATURE POTENTIELLE HUMIDE -REVERSIBLE-
! SONT GARDES CONSTANTS DURANT L'OPERATION.
!
! COMPUTATION OF THE SATURATION TEMPERATURE (IN GENERAL BUT NOT
! SYSTEMATICALLY LIFTING CONDENSATION TEMPERATURE). THE TOTAL MIXING
! RATIO AND THE MOIST -REVERSIBLE- POTENTIAL TEMPERATURE ARE KEPT
! CONSTANT DURING THE PROCESS.
!
! - TEMPORAIRES .
!
! ZRT        : RAPPORT DE MELANGE TOTAL DE L'EAU.
! : TOTAL WATER MIXING RATIO.
! ZCPT       : PARTIE "CP" DU "KAPPA" IMPLICITE DE LA TEMP. CONSERVEE.
! : "CP" PART OF THE IMPLICIT "KAPPA" OF THE CONSERVED TEMP..
! ZRRT       : PARTIE "R" DU "KAPPA" IMPLICITE DE LA TEMP. CONSERVEE.
! : "R" PART OF THE IMPLICIT "KAPPA" OF THE CONSERVED TEMP..
! ZS1        : EXPRESSION DE L'ENTROPIE ASSOCIE A LA TEMP. CONSERVEE.
! : ENTROPY'S EXPRESSION LINKED TO THE CONSERVED TEMPERATURE.
! ZCONS      : CONSTANTE AUXILIAIRE POUR LA BOUCLE DE NEWTON.
! : AUXILIARY CONSTANT FOR THE NEWTON LOOP.
! ZTITER     : EBAUCHE POUR LA SOLUTION DE LA BOUCLE DE NEWTON EN TEMP..
! : FIRST GUESS FOR THE SOLUTION OF THE NEWTON LOOP ON TEMP..
! ZQVSAT     :
! : SATURATED SPECIFIC HUMIDITY
!
! CALCULS PRELIMINAIRES.
! PRELIMINARY COMPUTATIONS.
!
!
! SECURITES.
! SECURITIES.
!
! QVSAT CALCULATION DEPENING ON THE SNOW OPTION.
!
ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,RTT-PT))
!
ZQL=MAX(0._JPRB,PQL)
ZQI=MAX(0._JPRB,PQI)
ZTIN=MAX(ZTMIN,MIN(ZTMAX,PT))
ZPRS=MAX(ZEPSP,PP)
ZQVSAT=QS(ZTIN,ZPRS)
ZQV=MAX(ZEPSRH*ZQVSAT,PQV)
!
ZRV=ZQV/(1._JPRB-ZQV-ZQL-ZQI)
ZRL=ZQL/(1._JPRB-ZQV-ZQL-ZQI)
ZRI=ZQI/(1._JPRB-ZQV-ZQL-ZQI)
ZRT=ZRV+ZRL+ZRI
ZCPT=RCPD+RCPV*ZRT
ZRRT=RD+RV*ZRT
ZE=(ZRV*ZPRS)/(ZRV+ZRDSV)
ZS1=ZCPT*LOG(ZTIN)-RD*LOG(ZPRS-ZE)-RV*ZRT &
& *LOG(ZE)-( FOLH (ZTIN,0._JPRB)*ZRL+ FOLH (ZTIN,1._JPRB)*ZRI)/ZTIN
ZCONS=ZS1+RD*LOG(ZRDSV/ZRT)
ZTITER=ZTIN
!
! BOUCLE DE NEWTON.
! NEWTON LOOP.
!
DO JIT=1,NBITER
    !
    ! CALCULS DEPENDANT DE L'OPTION NEIGE.
    ! SNOW OPTION DEPENDENT CALCULATIONS.
    !
    ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,RTT-ZTITER))
    !
    ZEW= ES (ZTITER)
    ZDLEW= FODLES (ZTITER)
    ZF=ZCONS+ZRRT*LOG(ZEW)-ZCPT*LOG(ZTITER)
    ZDF=ZRRT*ZDLEW-ZCPT/ZTITER
    ZTITER=ZTITER-ZF/ZDF
ENDDO
!
! *
! ------------------------------------------------------------------
! III - CALCUL DE LA PRESSION CORRESPONDANT AU POINT TRIPLE LE LONG
! DE L'ADIABATIQUE SATUREE IRREVERSIBLE PASSANT PAR LE POINT CALCULE
! PRECEDEMMENT. DANS LE CAS "LNEIGE=.T." LA TEMPERATURE DU POINT EN
! QUESTION DETERMINE LE CHOIX DES PARAMETRES LIES A "L" ET "CP".
!
! COMPUTATION OF PRESSURE CORRESPONDING TO THE TRIPLE POINT ON
! THE IRREVERSIBLE SATURATED ADIABAT PASSING THROUGH THE PREVIOUSLY
! OBTAINED POINT. IN THE "LNEIGE=.T." CASE THE LATTER'S TEMPERATURE
! DETERMINES THE CHOICE OF THE PARAMETERS LINKED TO "L" AND "CP".
!
!
! - TEMPORAIRES .
!
! ZDEL       : "MEMOIRE" DE LA VALEUR ZDELTA (EAU 0 / GLACE 1).
! : "MEMORY" OF THE ZDELTA (WATER 0 / ICE 1) VALUE.
! ZFUNCT     : EXPRESSION UTILISEE DANS LA BOUCLE DE NEWTON.
! : FUNCTIONAL EXPRESSION USED IN THE NEWTON LOOP.
! ZS2        : EXPRESSION DE LA PSEUDO ENTROPIE DE L'AD. IRREVERSIBLE.
! : PSEUDO ENTROPY'S EXPRESSION FOR THE IRREVERSIBLE AD..
! ZPITER     : EBAUCHE POUR LA SOLUTION DE LA BOUCLE DE NEWTON EN PRES..
! : FIRST GUESS FOR THE SOLUTION OF THE NEWTON LOOP ON PRES..
!
! CALCULS PRELIMINAIRES.
! PRELIMINARY COMPUTATIONS.
!
!
! CALCULS DEPENDANT DE L'OPTION NEIGE.
! SNOW OPTION DEPENDENT CALCULATIONS.
!
ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,RTT-ZTITER))
!
ZDEL=ZDELTA
ZEW= ES (ZTITER)
ZLSTT=ZLSTTW+ZDELTA*(ZLSTTI-ZLSTTW)
ZCWS=RCW+ZDELTA*(RCS-RCW)
ZLSTC= FOLH (ZTITER,ZDELTA)/ZTITER+ZCWS*RV &
& /( FOLH (ZTITER,ZDELTA)/ZTITER)
ZFUNCT=ZRDSV*RESTT*ZLSTT
ZS2=ZS1+ZRT*(ZLSTC+RV*LOG(ZEW)-RCPV &
& *LOG(ZTITER))
ZCONS=ZS2-ZCPLNTT
ZKAPPA=RKAPPA*(1._JPRB+ZRT* FOLH (ZTITER,ZDELTA)/(RD &
& *ZTITER))/(1._JPRB+RKAPPA*ZRT &
& * FOLH (ZTITER,ZDELTA)**2/(RD*RV*ZTITER**2))
ZPITER=(ZRDSV*ZEW/ZRT)*(RTT/ZTITER)**(1._JPRB/ZKAPPA) &
& -RESTT
ZPITER=MAX(ZPITER,ZEPSP)
!
! BOUCLE DE NEWTON (UNE ITERATION DE PLUS POUR P QUE POUR T).
! NEWTON LOOP (ONE MORE ITERATION FOR P THAN FOR T).
!
DO JIT=1,NBITER+1
    ZF=ZCONS+RD*LOG(ZPITER)-ZFUNCT/ZPITER
    ZDF=(RD*ZPITER+ZFUNCT)/ZPITER**2
    ZPITER=ZPITER-ZF/ZDF
ENDDO
!
! RETOUR A LA PRESSION REELLE.
! RETURN TO THE REAL PRESSURE.
!
ZPITER=ZPITER+RESTT
!
! *
! ------------------------------------------------------------------
! IV - CALCUL DE LA TEMPERATURE CORRESPONDANT A P STANDARD LE LONG
! DE L'ADIABATIQUE SATUREE IRREVERSIBLE PASSANT PAR LE POINT CALCULE
! PRECEDEMMENT. DANS LE CAS "LNEIGE=.T." LA PRESSION DU POINT EN
! QUESTION DETERMINE LE CHOIX DES PARAMETRES LIES A "L" ET "CP".
!
! COMPUTATION OF THE TEMPERATURE CORRESPONDING TO THE STD. P ON
! THE IRREVERSIBLE SATURATED ADIABAT PASSING THROUGH THE PREVIOUSLY
! OBTAINED POINT. IN THE "LNEIGE=.T." CASE THE LATTER'S PRESSURE
! DETERMINES THE CHOICE OF THE PARAMETERS LINKED TO "L" AND "CP".
!
! CALCULS PRELIMINAIRES.
! PRELIMINARY COMPUTATIONS.
!
!
! CALCULS DEPENDANT DE L'OPTION NEIGE.
! SNOW OPTION DEPENDENT CALCULATIONS.
!
ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,ZPITER-RATM))
!
ZDLSTT=(ZDELTA-ZDEL)*(ZLSTTI-ZLSTTW)
ZDEL=ZDELTA
ZS2=ZS2+ZDLSTT*ZRDSV*RESTT/(ZPITER-RESTT)
ZKAPPA=RKAPPA*(1._JPRB+(ZKNW+ZDELTA*(ZKNI-ZKNW))/ZPITER)/(1._JPRB &
& +(ZKDW+ZDELTA*(ZKDI-ZKDW))/ZPITER)
ZTITER=RTT*(RATM/ZPITER)**ZKAPPA
!
! BOUCLE DE NEWTON.
! NEWTON LOOP.
!
DO JIT=1,NBITER
    ZEW= ES(ZTITER)
    ZCWS=RCW+ZDEL*(RCS-RCW)
    ZLHZ=ZLHZW+ZDEL*(ZLHZI-ZLHZW)
    ZLH= FOLH (ZTITER,ZDEL)
    ZLSTC=ZLH/ZTITER+ZCWS*RV/(ZLH/ZTITER)
    ZRS=ZRDSV*ZEW/(RATM-ZEW)
    ZF=ZS2-RCPD*LOG(ZTITER)+RD*LOG(RATM-ZEW)-ZRS*ZLSTC
    ZDF=-RCPD/ZTITER-ZRS*((RATM/(RATM-ZEW))*ZLSTC*ZLH/(RV &
&     *ZTITER**2)+ZCWS*RV*ZLHZ/ZLH**2+(RCPV-ZCWS) &
&     /ZTITER)
    ZTITER=ZTITER-ZF/ZDF
ENDDO
!
! STOCKAGE DU RESULTAT.
! RESULT'S STORAGE.
!
THETAPW=ZTITER
END FUNCTION THETAPW
FUNCTION THETAV(PP,PT,PQV)
! --------------------------------------------------------------
! **** *thetav* Fonction température virtuelle potentielle.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
! En sortie:
!   thetav en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RD
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: THETAV,TV

THETAV=TV(PT,PQV)*(RATM/PP)**(RD/RCPD)
END FUNCTION THETAV
FUNCTION THETAL_KE(PP,PT,PQV,PQC)
! --------------------------------------------------------------
! **** *thetal* Fonction température potentielle avec prise en compte de la flottabilité des condensats.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
!   Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press, 200 Madison Avenue, New York.
!   thetal calculé ici est celui associé à une transformation adiabatique humide réversible.
!   Si pqc est non nul, il est donc plus consistant d'avoir pqv=qsat.
! Externes:
! Auteur/author:   2001-07, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
!   pqc humidité spécifique de la somme de tous les condensats ayant atteint leur vitesse limite (qliq+qice+qrain+...) (sans dimension).
! En sortie:
!   thetal en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RCPV
USE CONST_THER, ONLY : RV
USE CONST_THER, ONLY : RTT
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PQC
REAL(KIND=JPRB) :: THETAL_KE,FOLH
REAL(KIND=JPRB) :: ZRT,ZRL,ZKAPPA,ZGAMMA,ZDELARG,ZLV,ZMULT1,ZMULT2,ZMULT3
!
!-------------------------------------------------
! zrt: rapport de mélange total: vapeur + liq + glace + pluie + ...
!-------------------------------------------------
!
ZRT=(PQV+PQC)/(1._JPRB-PQV-PQC)
!
!-------------------------------------------------
! zrl: rapport de mélange des condensats.
!-------------------------------------------------
!
ZRL=PQC/(1._JPRB-PQV-PQC)
!
!-------------------------------------------------
! Kappa equ. (4.5.16) p122 Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press.
!-------------------------------------------------
!
ZKAPPA=(RD+ZRT*RV)/(RCPD+ZRT*RCPV)
!
!-------------------------------------------------
! Gamma equ. (4.5.16) p122 Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press.
!-------------------------------------------------
!
ZGAMMA=ZRT*RV/(RCPD+ZRT*RCPV)
!
!-------------------------------------------------
! zlv: chaleur latente.
!-------------------------------------------------
!
ZDELARG=MAX(0._JPRB,SIGN(1._JPRB,RTT-PT))
ZLV=FOLH(PT,ZDELARG)
!
!-------------------------------------------------
! thetal equ. (4.5.15) p121 Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press.
!-------------------------------------------------
!
ZMULT1=(1._JPRB-ZRL/(RD/RV+ZRT))**ZKAPPA
IF(ZRL < ZRT) THEN
    ZMULT2=(1._JPRB-ZRL/ZRT)**(-ZGAMMA)
ELSE
    ZMULT2=1._JPRB
ENDIF
ZMULT3=EXP(-ZLV*ZRL/((RCPD+ZRT*RCPV)*PT))
THETAL_KE=PT*(RATM/PP)**ZKAPPA &
& *ZMULT1*ZMULT2*ZMULT3
END FUNCTION THETAL_KE
FUNCTION THETAVL(PP,PT,PQV,PQC)
! --------------------------------------------------------------
! **** *thetavl* Fonction température virtuelle potentielle avec prise en compte de la flottabilité des condensats.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
!   Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press, 200 Madison Avenue, New York.
!   Thetavl calculé ici est celui associé à une transformation adiabatique humide réversible.
!   Si pqc est non nul, il est donc plus consistant d'avoir pqv=qsat.
! Externes:
! Auteur/author:   2001-07, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
!   pqc humidité spécifique de la somme de tous les condensats ayant atteint leur vitesse limite (qliq+qice+qrain+...) (sans dimension).
! En sortie:
!   thetavl en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RATM
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RCPV
USE CONST_THER, ONLY : RV
USE CONST_THER, ONLY : RTT
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PQC
REAL(KIND=JPRB) :: THETAVL,TV,FOLH
REAL(KIND=JPRB) :: ZRT,ZRL,ZKAPPA,ZGAMMA,ZDELARG,ZLV,ZMULT1,ZMULT2,ZMULT3,ZMULT4
!
!-------------------------------------------------
! zrt: rapport de mélange total: vapeur + liq + glace + pluie + ...
!-------------------------------------------------
!
ZRT=(PQV+PQC)/(1._JPRB-PQV-PQC)
!
!-------------------------------------------------
! zrl: rapport de mélange des condensats.
!-------------------------------------------------
!
ZRL=PQC/(1._JPRB-PQV-PQC)
!
!-------------------------------------------------
! Kappa equ. (4.5.16) p122 Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press.
!-------------------------------------------------
!
ZKAPPA=(RD+ZRT*RV)/(RCPD+ZRT*RCPV)
!
!-------------------------------------------------
! Gamma equ. (4.5.16) p122 Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press.
!-------------------------------------------------
!
ZGAMMA=ZRT*RV/(RCPD+ZRT*RCPV)
!
!-------------------------------------------------
! zlv: chaleur latente.
!-------------------------------------------------
!
ZDELARG=MAX(0._JPRB,SIGN(1._JPRB,RTT-PT))
ZLV=FOLH(PT,ZDELARG)
!
!-------------------------------------------------
! Thetavl equ. (4.5.18) p122 Kerry A. Emanuel, "Atmospheric Convection", 1994 Oxford University Press.
!-------------------------------------------------
!
ZMULT1=(1._JPRB-ZRL/(1._JPRB+ZRT))
ZMULT2=(1._JPRB-ZRL/(RD/RV+ZRT))**(ZKAPPA-1._JPRB)
IF(ZRL < ZRT) THEN
    ZMULT3=(1._JPRB-ZRL/ZRT)**(-ZGAMMA)
ELSE
    ZMULT3=1._JPRB
ENDIF
ZMULT4=EXP(-ZLV*ZRL/((RCPD+ZRT*RCPV)*PT))
THETAVL=TV(PT,PQV)*(RATM/PP)**ZKAPPA &
& *ZMULT1*ZMULT2*ZMULT3*ZMULT4
END FUNCTION THETAVL
FUNCTION TV(PT,PQV)
! --------------------------------------------------------------
! **** *tv* Fonction température virtuelle, aussi appelée température de densité.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-09, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
! En sortie:
!   tv en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RV
IMPLICIT NONE
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: ZRV
REAL(KIND=JPRB) :: PT,TV
!
!-------------------------------------------------
! Calcul des rapports de mélange.
!-------------------------------------------------
!
ZRV=PQV/(1._JPRB-PQV)
!
!-------------------------------------------------
! Expression de Tv.
!-------------------------------------------------
!
TV=PT*(1._JPRB+ZRV*RV/RD)/(1._JPRB+ZRV)
END FUNCTION TV
FUNCTION WETPOINT(PP,PT,PQV)
! --------------------------------------------------------------
! **** *wetpoint* TEMPERATURE PSEUDO-ADIABATIQUE DU THERMOMETRE MOUILLE (point "bleu").
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur/author:   2000-10, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pp pression en Pa.
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
! En sortie:
!   wetpoint (K).
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
IMPLICIT NONE
REAL(KIND=JPRB) :: PP
REAL(KIND=JPRB) :: PQV
REAL(KIND=JPRB) :: PT
REAL(KIND=JPRB) :: ZPCOND
REAL(KIND=JPRB) :: ZTCOND,WETPOINT,T_EVOL_AD_HUM

!
!-------------------------------------------------
! On cherche le point de condensation.
!-------------------------------------------------
!
CALL POINT_CONDENS(PT,PQV,PP,ZTCOND,ZPCOND)
!
!-------------------------------------------------
! On revient au niveau courant
! selon une adiabatique humide.
!-------------------------------------------------
!
WETPOINT=T_EVOL_AD_HUM(ZTCOND,ZPCOND,PP)
END FUNCTION WETPOINT
FUNCTION WETPOINT_ARP(PP,PT,PQV,PQL,PQI)
! --------------------------------------------------------------
! **** *wetpoint_arp* TEMPERATURE PSEUDO-ADIABATIQUE DU THERMOMETRE MOUILLE (point "bleu").
! --------------------------------------------------------------
! Sujet:
! CALCUL DE LA TEMPERATURE PSEUDO-ADIABATIQUE DU THERMOMETRE MOUILLE.
! Arguments explicites:
! Arguments implicites:
!
! Methode.
! --------
!
! CETTE ROUTINE N'A QU'UNE DIMENSION POUR SES VARIABLES D'ENTREE
! AFIN D'ETRE LA PLUS GENERALE POSSIBLE (UTILISATION SUR LES NIVEAUX
! DU MODELE COMME DU POST-PROCESSING, PAR EXEMPLE). TOUT ETAT DE
! L'AIR REALISTE EST ADMIS EN ENTREE ET L'ALGORITHME PREND EN COMPTE
! AUTOMATIQUEMENT UNE POSSIBLE TRANSITION DE PHASE LIQUIDE/GLACE.
! TROIS EQUATIONS IMPLICITES SONT RESOLUES PAR METHODE DE NEWTON:
! - RECHERCHE DU POINT DE SATURATION D'ENTROPIE EGALE PAR
! TRANSFORMATION REVERSIBLE ;
! - RECHERCHE DU POINT DE TEMPERATURE EGALE A CELLE DU POINT
! TRIPLE LE LONG DE L'ADIABATIQUE SATUREE IRREVERSIBLE ;
! - RECHERCHE DU POINT DE PRESSION EGALE A LA REFERENCE
! ATMOSPHERIQUE LE LONG D'UNE AUTRE (PARFOIS LA MEME) ADIABATIQUE
! IRREVERSIBLE.
! REMARQUES :
! - POUR LA PREMIERE ETAPE LA FORME SYMETRIQUE DE L'ENTROPIE
! HUMIDE PROPOSEE PAR P. MARQUET EST UTILISEE AFIN DE PERMETTRE UN
! MELANGE DE PHASES LIQUIDE ET GLACE DANS L'ETAT DE L'AIR ;
! - POUR LES DEUX DERNIERES ETAPES, PLUTOT QUE DE NEGLIGER
! COMME DE COUTUME LE TERME CONTENANT LE CP DU CONDENSAT, L'AUTEUR
! DE LA ROUTINE EN A DERIVE UNE APPROXIMATION QUASI-EXACTE ET PLUTOT
! BON MARCHE ;
! - POUR CES DEUX MEMES ETAPES, LES EBAUCHES DES BOUCLES DE
! NEWTON SONT OBTENUES PAR EXTRAPOLATION D'UNE LINEARISATION LOCALE
! APPROCHEE DES EQUATIONS ADIABATIQUE SATUREES.
!
! THIS ROUTINE HAS ONLY ONE DIMENSIONAL INPUT/OUTPUT ARRAYS IN
! ORDER TO BE THE MOST GENERAL POSSIBLE (USE ON MODEL OR ON POST-
! PROCESSING LEVELS, FOR EXAMPLE). ALL POSSIBLE REALISTIC INPUT
! STATES ARE ALLOWED AND THE ALGORITHM AUTOMATICALLY TAKES INTO
! ACCOUNT THE POTENTIAL LIQUID/ICE WATER TRANSITION.
! THREE IMPLICIT EQUATIONS ARE SOLVED BY NEWTON METHODS :
! - SEARCH OF THE SATURATION POINT OF EQUAL ENTROPY UNDER A
! REVERSIBLE TRANSFORM ;
! - SEARCH OF THE POINT OF TEMPERATURE EQUAL TO THAT OF THE
! TRIPLE POINT ALONG THE IRREVERSIBLE MOIST ADIABAT ;
! - SEARCH OF THE POINT OF REFERENCE ATMOSPHERIC PRESSURE
! ALONG ANOTHER (SOMETIMES IDENTICAL) IRREVERSIBLE MOIST ADIABAT.
! REMARKS :
! - FOR THE FIRST STEP THE SYMETRIC FORM OF THE MOIST ENTROPY
! PROPOSED BY P. MARQUET IS USED IN ORDER TO ALLOW A MIX OF LIQUID
! AND ICE WATER IN THE ATMOSPHERIC STATE ;
! - FOR THE TWO LAST STEPS, RATHER THAN THE USUAL NEGLECTION
! OF THE TERM MULTIPLIED BY CP OF THE CONDENSATE, THE ROUTINE'S
! AUTHOR DERIVED A QUASI EXACT AND NOT TOO EXPENSIVE ANALYTICAL
! APPROXIMATION FOR IT ;
! - FOR THE SAME STEPS, THE GUESSES OF THE NEWTON LOOP ARE
! OBTAINED BY VERTICAL EXTRAPOLATION OF A LINEAR LOCAL APPROXIMATION
! OF THE MOIST ADIABATS.
!
! Auteur/author: 92-09, J.F. Geleyn.
!
! Modifications.
! --------------
! 96-04, J. Calvo: Introduced a minimun in RH instead of a mini-
! mun in PQV. Added a security threshold in the
! calculation of the triple  point pressure
! first guess.
! --------------------------------------------------------------
! En entree:
! En sortie:
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RCPD
USE CONST_THER, ONLY : RCPV
USE CONST_THER, ONLY : RCS
USE CONST_THER, ONLY : RCW
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RESTT
USE CONST_THER, ONLY : RKAPPA
USE CONST_THER, ONLY : RTT
USE CONST_THER, ONLY : RV,NBITER
IMPLICIT NONE
INTEGER(KIND=JPIM) :: JIT
REAL(KIND=JPRB) :: ZCPLNTT
REAL(KIND=JPRB) :: ZCWS
REAL(KIND=JPRB) :: ZDELTA
REAL(KIND=JPRB) :: ZDF
REAL(KIND=JPRB) :: ZDLEW
REAL(KIND=JPRB) :: ZDLSTT
REAL(KIND=JPRB) :: ZE
REAL(KIND=JPRB) :: ZEPSP
REAL(KIND=JPRB) :: ZEPSRH
REAL(KIND=JPRB) :: ZEW
REAL(KIND=JPRB) :: ZF
REAL(KIND=JPRB) :: ZKAPPA
REAL(KIND=JPRB) :: ZKDI
REAL(KIND=JPRB) :: ZKDW
REAL(KIND=JPRB) :: ZKNI
REAL(KIND=JPRB) :: ZKNW
REAL(KIND=JPRB) :: ZLH
REAL(KIND=JPRB) :: ZLHZ
REAL(KIND=JPRB) :: ZLHZI
REAL(KIND=JPRB) :: ZLHZW
REAL(KIND=JPRB) :: ZLSTC
REAL(KIND=JPRB) :: ZLSTT
REAL(KIND=JPRB) :: ZLSTTI
REAL(KIND=JPRB) :: ZLSTTW
REAL(KIND=JPRB) :: ZPRS
REAL(KIND=JPRB) :: ZQI
REAL(KIND=JPRB) :: ZQL
REAL(KIND=JPRB) :: ZQV
REAL(KIND=JPRB) :: ZQVSAT
REAL(KIND=JPRB) :: ZRDSV
REAL(KIND=JPRB) :: ZRI
REAL(KIND=JPRB) :: ZRL
REAL(KIND=JPRB) :: ZRS
REAL(KIND=JPRB) :: ZRV
REAL(KIND=JPRB) :: ZTIN
REAL(KIND=JPRB) :: ZTMAX
REAL(KIND=JPRB) :: ZTMIN,WETPOINT_ARP

REAL(KIND=JPRB) PP,PT,PQV,PQL,PQI,FOLH,QS,ES,FODLES
REAL(KIND=JPRB) ZRT,ZCPT,ZRRT,ZS1,ZCONS,ZTITER,ZDEL,ZS2,ZFUNCT,ZPITER
!
! *
! ------------------------------------------------------------------
! I - CALCUL DES const_ther DERIVEES DE CELLES ISSUES DE YOMCST ET
! const_ther DE SECURITE.
!
! COMPUTATION OF CONSTANTS DERIVED FROM THOSE CONTAINED IN
! YOMCST AND SECURITY CONSTANTS.
!
ZRDSV=RD/RV
ZLSTTW= FOLH (RTT,0._JPRB)/RTT+RCW*RV/( FOLH (RTT,0._JPRB)/RTT)
ZLSTTI= FOLH (RTT,1._JPRB)/RTT+RCS*RV/( FOLH (RTT,1._JPRB)/RTT)
ZLHZW= FOLH (0._JPRB,0._JPRB)
ZLHZI= FOLH (0._JPRB,1._JPRB)
ZKNW=RESTT* FOLH (RTT,0._JPRB)/(RV*RTT)
ZKNI=RESTT* FOLH (RTT,1._JPRB)/(RV*RTT)
ZKDW=RKAPPA*RESTT*( FOLH (RTT,0._JPRB)/(RV*RTT))**2
ZKDI=RKAPPA*RESTT*( FOLH (RTT,1._JPRB)/(RV*RTT))**2
ZCPLNTT=RCPD*LOG(RTT)
!
ZEPSP=10.
ZEPSRH=0.001
ZTMIN=155.
ZTMAX=355.
!
!
! *
! ------------------------------------------------------------------
! II - CALCUL DE LA TEMPERATURE DE SATURATION (EN GENERAL MAIS PAS
! FORCEMENT TEMPERATURE DU POINT DE CONDENSATION). LE RAPPORT DE
! MELANGE TOTAL ET LA TEMPERATURE POTENTIELLE HUMIDE -REVERSIBLE-
! SONT GARDES CONSTANTS DURANT L'OPERATION.
!
! COMPUTATION OF THE SATURATION TEMPERATURE (IN GENERAL BUT NOT
! SYSTEMATICALLY LIFTING CONDENSATION TEMPERATURE). THE TOTAL MIXING
! RATIO AND THE MOIST -REVERSIBLE- POTENTIAL TEMPERATURE ARE KEPT
! CONSTANT DURING THE PROCESS.
!
! - TEMPORAIRES .
!
! ZRT        : RAPPORT DE MELANGE TOTAL DE L'EAU.
! : TOTAL WATER MIXING RATIO.
! ZCPT       : PARTIE "CP" DU "KAPPA" IMPLICITE DE LA TEMP. CONSERVEE.
! : "CP" PART OF THE IMPLICIT "KAPPA" OF THE CONSERVED TEMP..
! ZRRT       : PARTIE "R" DU "KAPPA" IMPLICITE DE LA TEMP. CONSERVEE.
! : "R" PART OF THE IMPLICIT "KAPPA" OF THE CONSERVED TEMP..
! ZS1        : EXPRESSION DE L'ENTROPIE ASSOCIE A LA TEMP. CONSERVEE.
! : ENTROPY'S EXPRESSION LINKED TO THE CONSERVED TEMPERATURE.
! ZCONS      : CONSTANTE AUXILIAIRE POUR LA BOUCLE DE NEWTON.
! : AUXILIARY CONSTANT FOR THE NEWTON LOOP.
! ZTITER     : EBAUCHE POUR LA SOLUTION DE LA BOUCLE DE NEWTON EN TEMP..
! : FIRST GUESS FOR THE SOLUTION OF THE NEWTON LOOP ON TEMP..
! ZQVSAT     :
! : SATURATED SPECIFIC HUMIDITY
!
! CALCULS PRELIMINAIRES.
! PRELIMINARY COMPUTATIONS.
!
!
! SECURITES.
! SECURITIES.
!
! QVSAT CALCULATION DEPENING ON THE SNOW OPTION.
!
ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,RTT-PT))
!
ZQL=MAX(0._JPRB,PQL)
ZQI=MAX(0._JPRB,PQI)
ZTIN=MAX(ZTMIN,MIN(ZTMAX,PT))
ZPRS=MAX(ZEPSP,PP)
ZQVSAT=QS(ZTIN,ZPRS)
ZQV=MAX(ZEPSRH*ZQVSAT,PQV)
!
ZRV=ZQV/(1._JPRB-ZQV-ZQL-ZQI)
ZRL=ZQL/(1._JPRB-ZQV-ZQL-ZQI)
ZRI=ZQI/(1._JPRB-ZQV-ZQL-ZQI)
ZRT=ZRV+ZRL+ZRI
ZCPT=RCPD+RCPV*ZRT
ZRRT=RD+RV*ZRT
ZE=(ZRV*ZPRS)/(ZRV+ZRDSV)
ZS1=ZCPT*LOG(ZTIN)-RD*LOG(ZPRS-ZE)-RV*ZRT &
& *LOG(ZE)-( FOLH (ZTIN,0._JPRB)*ZRL+ FOLH (ZTIN,1._JPRB)*ZRI)/ZTIN
ZCONS=ZS1+RD*LOG(ZRDSV/ZRT)
ZTITER=ZTIN
!
! BOUCLE DE NEWTON.
! NEWTON LOOP.
!
DO JIT=1,NBITER
    !
    ! CALCULS DEPENDANT DE L'OPTION NEIGE.
    ! SNOW OPTION DEPENDENT CALCULATIONS.
    !
    ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,RTT-ZTITER))
    !
    ZEW= ES (ZTITER)
    ZDLEW= FODLES (ZTITER)
    ZF=ZCONS+ZRRT*LOG(ZEW)-ZCPT*LOG(ZTITER)
    ZDF=ZRRT*ZDLEW-ZCPT/ZTITER
    ZTITER=ZTITER-ZF/ZDF
ENDDO
!
! *
! ------------------------------------------------------------------
! III - CALCUL DE LA PRESSION CORRESPONDANT AU POINT TRIPLE LE LONG
! DE L'ADIABATIQUE SATUREE IRREVERSIBLE PASSANT PAR LE POINT CALCULE
! PRECEDEMMENT. DANS LE CAS "LNEIGE=.T." LA TEMPERATURE DU POINT EN
! QUESTION DETERMINE LE CHOIX DES PARAMETRES LIES A "L" ET "CP".
!
! COMPUTATION OF PRESSURE CORRESPONDING TO THE TRIPLE POINT ON
! THE IRREVERSIBLE SATURATED ADIABAT PASSING THROUGH THE PREVIOUSLY
! OBTAINED POINT. IN THE "LNEIGE=.T." CASE THE LATTER'S TEMPERATURE
! DETERMINES THE CHOICE OF THE PARAMETERS LINKED TO "L" AND "CP".
!
! - TEMPORAIRES .
!
! ZDEL       : "MEMOIRE" DE LA VALEUR ZDELTA (EAU 0 / GLACE 1).
! : "MEMORY" OF THE ZDELTA (WATER 0 / ICE 1) VALUE.
! ZFUNCT     : EXPRESSION UTILISEE DANS LA BOUCLE DE NEWTON.
! : FUNCTIONAL EXPRESSION USED IN THE NEWTON LOOP.
! ZS2        : EXPRESSION DE LA PSEUDO ENTROPIE DE L'AD. IRREVERSIBLE.
! : PSEUDO ENTROPY'S EXPRESSION FOR THE IRREVERSIBLE AD..
! ZPITER     : EBAUCHE POUR LA SOLUTION DE LA BOUCLE DE NEWTON EN PRES..
! : FIRST GUESS FOR THE SOLUTION OF THE NEWTON LOOP ON PRES..
!
! CALCULS PRELIMINAIRES.
! PRELIMINARY COMPUTATIONS.
!
!
! CALCULS DEPENDANT DE L'OPTION NEIGE.
! SNOW OPTION DEPENDENT CALCULATIONS.
!
ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,RTT-ZTITER))
!
ZDEL=ZDELTA
ZEW= ES (ZTITER)
ZLSTT=ZLSTTW+ZDELTA*(ZLSTTI-ZLSTTW)
ZCWS=RCW+ZDELTA*(RCS-RCW)
ZLSTC= FOLH (ZTITER,ZDELTA)/ZTITER+ZCWS*RV &
& /( FOLH (ZTITER,ZDELTA)/ZTITER)
ZFUNCT=ZRDSV*RESTT*ZLSTT
ZS2=ZS1+ZRT*(ZLSTC+RV*LOG(ZEW)-RCPV &
& *LOG(ZTITER))
ZCONS=ZS2-ZCPLNTT
ZKAPPA=RKAPPA*(1._JPRB+ZRT* FOLH (ZTITER,ZDELTA)/(RD &
& *ZTITER))/(1._JPRB+RKAPPA*ZRT &
& * FOLH (ZTITER,ZDELTA)**2/(RD*RV*ZTITER**2))
ZPITER=(ZRDSV*ZEW/ZRT)*(RTT/ZTITER)**(1._JPRB/ZKAPPA) &
& -RESTT
ZPITER=MAX(ZPITER,ZEPSP)
!
! BOUCLE DE NEWTON (UNE ITERATION DE PLUS POUR P QUE POUR T).
! NEWTON LOOP (ONE MORE ITERATION FOR P THAN FOR T).
!
DO JIT=1,NBITER+1
    ZF=ZCONS+RD*LOG(ZPITER)-ZFUNCT/ZPITER
    ZDF=(RD*ZPITER+ZFUNCT)/ZPITER**2
    ZPITER=ZPITER-ZF/ZDF
ENDDO
!
! RETOUR A LA PRESSION REELLE.
! RETURN TO THE REAL PRESSURE.
!
ZPITER=ZPITER+RESTT
!
! *
! ------------------------------------------------------------------
! IV - CALCUL DE LA TEMPERATURE CORRESPONDANT A P STANDARD LE LONG
! DE L'ADIABATIQUE SATUREE IRREVERSIBLE PASSANT PAR LE POINT CALCULE
! PRECEDEMMENT. DANS LE CAS "LNEIGE=.T." LA PRESSION DU POINT EN
! QUESTION DETERMINE LE CHOIX DES PARAMETRES LIES A "L" ET "CP".
!
! COMPUTATION OF THE TEMPERATURE CORRESPONDING TO THE STD. P ON
! THE IRREVERSIBLE SATURATED ADIABAT PASSING THROUGH THE PREVIOUSLY
! OBTAINED POINT. IN THE "LNEIGE=.T." CASE THE LATTER'S PRESSURE
! DETERMINES THE CHOICE OF THE PARAMETERS LINKED TO "L" AND "CP".
!
! CALCULS PRELIMINAIRES.
! PRELIMINARY COMPUTATIONS.
!
!
! CALCULS DEPENDANT DE L'OPTION NEIGE.
! SNOW OPTION DEPENDENT CALCULATIONS.
!
ZDELTA=MAX(0._JPRB,SIGN(1._JPRB,ZPITER-PP))
!
ZDLSTT=(ZDELTA-ZDEL)*(ZLSTTI-ZLSTTW)
ZDEL=ZDELTA
ZS2=ZS2+ZDLSTT*ZRDSV*RESTT/(ZPITER-RESTT)
ZKAPPA=RKAPPA*(1._JPRB+(ZKNW+ZDELTA*(ZKNI-ZKNW))/ZPITER)/(1._JPRB &
& +(ZKDW+ZDELTA*(ZKDI-ZKDW))/ZPITER)
ZTITER=RTT*(PP/ZPITER)**ZKAPPA
!
! BOUCLE DE NEWTON.
! NEWTON LOOP.
!
DO JIT=1,NBITER
    ZEW= ES(ZTITER)
    ZCWS=RCW+ZDEL*(RCS-RCW)
    ZLHZ=ZLHZW+ZDEL*(ZLHZI-ZLHZW)
    ZLH= FOLH (ZTITER,ZDEL)
    ZLSTC=ZLH/ZTITER+ZCWS*RV/(ZLH/ZTITER)
    ZRS=ZRDSV*ZEW/(PP-ZEW)
    ZF=ZS2-RCPD*LOG(ZTITER)+RD*LOG(PP-ZEW)-ZRS*ZLSTC
    ZDF=-RCPD/ZTITER-ZRS*((PP/(PP-ZEW))*ZLSTC*ZLH/(RV &
&     *ZTITER**2)+ZCWS*RV*ZLHZ/ZLH**2+(RCPV-ZCWS) &
&     /ZTITER)
    ZTITER=ZTITER-ZF/ZDF
ENDDO
!
! STOCKAGE DU RESULTAT.
! RESULT'S STORAGE.
!
WETPOINT_ARP=ZTITER
END FUNCTION WETPOINT_ARP
FUNCTION TRHO(PT,PQV,PQC)
! --------------------------------------------------------------
! **** *trho* Fonction température de densité: température virtuelle plus effets de flottabilité des condensats.
! --------------------------------------------------------------
! Sujet:
! Arguments explicites:
! Arguments implicites:
! Methode:
! Externes:
! Auteur:   2001-07, J.M. Piriou.
! Modifications:
! --------------------------------------------------------------
! En entree:
!   pt température en K.
!   pqv humidité spécifique de la vapeur d'eau (sans dimension).
!   pqc humidité spécifique de la somme de tous les condensats ayant atteint leur vitesse limite (qliq+qice+qrain+...) (sans dimension).
! En sortie:
!   trho en K.
! --------------------------------------------------------------
USE PARKIND1  ,ONLY : JPIM ,JPRB
USE CONST_THER, ONLY : RD
USE CONST_THER, ONLY : RV
IMPLICIT NONE
REAL(KIND=JPRB) :: PQV,PQC
REAL(KIND=JPRB) :: ZRV,ZRC
REAL(KIND=JPRB) :: PT,TRHO
!
!-------------------------------------------------
! Calcul des rapports de mélange.
!-------------------------------------------------
!
ZRV=PQV/(1._JPRB-PQV-PQC)
ZRC=PQC/(1._JPRB-PQV-PQC)
!
!-------------------------------------------------
! Expression de trho.
!-------------------------------------------------
!
TRHO=PT*(1._JPRB+ZRV*RV/RD)/(1._JPRB+ZRV+ZRC)
END FUNCTION TRHO
