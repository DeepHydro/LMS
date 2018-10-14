!==============================================================================!

   SUBROUTINE ADV_UV_EDGE_GCN

!==============================================================================!
! this subroutine calculate advective, coriolis, pressure gradient, etc in     !
! x and y momentum equations except vertical diffusion terms for internal mode ! 
!==============================================================================!

   USE ALL_VARS
   USE BCS

#  if defined (WET_DRY)
   USE MOD_WD
#  endif

   USE MOD_HEATFLUX

#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif

   IMPLICIT NONE
   REAL(SP) :: XFLUX(0:NT,KB),YFLUX(0:NT,KB)
   REAL(SP) :: PSTX_TM(0:NT,KB),PSTY_TM(0:NT,KB)
   REAL(SP) :: COFA1,COFA2,COFA3,COFA4,COFA5,COFA6,COFA7,COFA8
   REAL(SP) :: XADV,YADV,TXXIJ,TYYIJ,TXYIJ
   REAL(SP) :: VISCOF,VISCOF1,VISCOF2,TEMP,TPA,TPB
   REAL(SP) :: XIJA,YIJA,XIJB,YIJB,UIJ,VIJ
   REAL(SP) :: DIJ,ELIJ,TMPA,TMPB,TMP,XFLUXV,YFLUXV
   REAL(SP) :: FACT,FM1,EXFLUX,ISWETTMP
   INTEGER  :: I,IA,IB,J1,J2,K1,K2,K3,K4,K5,K6,K,II,J,I1,I2
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UIJ1,VIJ1,UIJ2,VIJ2,FXX,FYY
   REAL(SP),ALLOCATABLE,DIMENSION(:) :: UALFA,VALFA
   REAL(SP) :: UALFA_TMP,VALFA_TMP
   INTEGER :: ERROR
   REAL(SP) :: EPS,EPS_TMP

!------------------------------------------------------------------------------!

   FACT = 0.0_SP
   FM1  = 1.0_SP
   IF(HORZMIX == 'closure') THEN
     FACT = 1.0_SP
     FM1  = 0.0_SP
   END IF

!
!-----Initialize Flux Variables------------------------------------------------!
!
   XFLUX  = 0.0_SP
   YFLUX  = 0.0_SP
   PSTX_TM = 0.0_SP
   PSTY_TM = 0.0_SP
# if defined (SEMI_IMPLICIT)
   XFLUX3= 0.0_SP
   YFLUX3= 0.0_SP
# endif

!
!-----Loop Over Edges and Accumulate Flux--------------------------------------!
!
   ALLOCATE(UIJ1(NE),VIJ1(NE),UIJ2(NE),VIJ2(NE))
   ALLOCATE(UALFA(0:NT),VALFA(0:NT))
   ALLOCATE(FXX(NE),FYY(NE))
   
   DO K=1,KBM1
     UIJ1=0.0_SP;VIJ1=0.0_SP;UIJ2=0.0_SP;VIJ2=0.0_SP
     UALFA=1.0_SP;VALFA=1.0_SP
     FXX=0.0_SP;FYY=0.0_SP
     
     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)
       J1=IENODE(I,1)
       J2=IENODE(I,2)
!       DIJ=0.5_SP*(DT(J1)+DT(J2))

       K1=NBE(IA,1)
       K2=NBE(IA,2)
       K3=NBE(IA,3)
       K4=NBE(IB,1)
       K5=NBE(IB,2)
       K6=NBE(IB,3)

       XIJA=XIJC(I)-XC(IA)
       YIJA=YIJC(I)-YC(IA)
       XIJB=XIJC(I)-XC(IB)
       YIJB=YIJC(I)-YC(IB)


       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
#      if !defined (SEMI_IMPLICIT)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      else
       IF(ISWETCT(IA) == 1 .OR. ISWETCT(IB) == 1)THEN
#      endif
#      endif
       COFA1=A1U(IA,1)*U(IA,K)+A1U(IA,2)*U(K1,K)+A1U(IA,3)*U(K2,K)+A1U(IA,4)*U(K3,K)
       COFA2=A2U(IA,1)*U(IA,K)+A2U(IA,2)*U(K1,K)+A2U(IA,3)*U(K2,K)+A2U(IA,4)*U(K3,K)
       COFA5=A1U(IA,1)*V(IA,K)+A1U(IA,2)*V(K1,K)+A1U(IA,3)*V(K2,K)+A1U(IA,4)*V(K3,K)
       COFA6=A2U(IA,1)*V(IA,K)+A2U(IA,2)*V(K1,K)+A2U(IA,3)*V(K2,K)+A2U(IA,4)*V(K3,K)

       UIJ1(I)=COFA1*XIJA+COFA2*YIJA
       VIJ1(I)=COFA5*XIJA+COFA6*YIJA
       
        EPS_TMP=ABS(UIJ1(I)+EPSILON(EPS))
        IF(EPS_TMP ==0)THEN
            UALFA_TMP=1.0_SP
        ELSE
              UALFA_TMP=ABS(U(IA,K)-U(IB,K))/EPS_TMP
        ENDIF
        
        EPS_TMP=ABS(VIJ1(I)+EPSILON(EPS))
        IF(EPS_TMP ==0)THEN
            VALFA_TMP=1.0_SP
        ELSE
             VALFA_TMP=ABS(V(IA,K)-V(IB,K))/EPS_TMP
        ENDIF
        
       !UALFA_TMP=ABS(U(IA,K)-U(IB,K))/ABS(UIJ1(I)+EPSILON(EPS))
       !VALFA_TMP=ABS(V(IA,K)-V(IB,K))/ABS(VIJ1(I)+EPSILON(EPS))
       
       IF(UALFA_TMP > 1.0_SP)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1.0_SP)VALFA_TMP = 1.0_SP
       UALFA(IA)=MIN(UALFA(IA),UALFA_TMP)
       VALFA(IA)=MIN(VALFA(IA),VALFA_TMP)
       
       COFA3=A1U(IB,1)*U(IB,K)+A1U(IB,2)*U(K4,K)+A1U(IB,3)*U(K5,K)+A1U(IB,4)*U(K6,K)
       COFA4=A2U(IB,1)*U(IB,K)+A2U(IB,2)*U(K4,K)+A2U(IB,3)*U(K5,K)+A2U(IB,4)*U(K6,K)
       COFA7=A1U(IB,1)*V(IB,K)+A1U(IB,2)*V(K4,K)+A1U(IB,3)*V(K5,K)+A1U(IB,4)*V(K6,K)
       COFA8=A2U(IB,1)*V(IB,K)+A2U(IB,2)*V(K4,K)+A2U(IB,3)*V(K5,K)+A2U(IB,4)*V(K6,K)

       UIJ2(I)=COFA3*XIJB+COFA4*YIJB
       VIJ2(I)=COFA7*XIJB+COFA8*YIJB
       
       EPS_TMP=ABS(UIJ2(I)+EPSILON(EPS))
        IF(EPS_TMP ==0)THEN
            UALFA_TMP=1.0_SP
        ELSE
              UALFA_TMP=ABS(U(IA,K)-U(IB,K))/EPS_TMP
        ENDIF
        
          EPS_TMP=ABS(VIJ2(I)+EPSILON(EPS))
        IF(EPS_TMP ==0)THEN
            VALFA_TMP=1.0_SP
        ELSE
             VALFA_TMP=ABS(V(IA,K)-V(IB,K))/EPS_TMP
        ENDIF
        
        
       !UALFA_TMP=ABS(U(IA,K)-U(IB,K))/ABS(UIJ2(I)+EPSILON(EPS))
       !VALFA_TMP=ABS(V(IA,K)-V(IB,K))/ABS(VIJ2(I)+EPSILON(EPS))
       IF(UALFA_TMP > 1.0_SP)UALFA_TMP = 1.0_SP
       IF(VALFA_TMP > 1.0_SP)VALFA_TMP = 1.0_SP
       UALFA(IB)=MIN(UALFA(IB),UALFA_TMP)
       VALFA(IB)=MIN(VALFA(IB),VALFA_TMP)
       
!
!-------ADD THE VISCOUS TERM & ADVECTION TERM---------------------------------!
!

       VISCOF1=ART(IA)*SQRT(COFA1**2+COFA6**2+0.5_SP*(COFA2+COFA5)**2)
       VISCOF2=ART(IB)*SQRT(COFA3**2+COFA8**2+0.5_SP*(COFA4+COFA7)**2)

!       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON
       VISCOF=FACT*0.5_SP*HORCON*(VISCOF1+VISCOF2)/HPRNU + FM1*HORCON/HPRNU

       TXXIJ=(COFA1+COFA3)*VISCOF
       TYYIJ=(COFA6+COFA8)*VISCOF
       TXYIJ=0.5_SP*(COFA2+COFA4+COFA5+COFA7)*VISCOF
       FXX(I)=DIJ*(TXXIJ*DLTYC(I)-TXYIJ*DLTXC(I))
       FYY(I)=DIJ*(TXYIJ*DLTYC(I)-TYYIJ*DLTXC(I))
#      if defined (WET_DRY)
       END IF
#      endif
     END DO

     DO I=1,NE
       IA=IEC(I,1)
       IB=IEC(I,2)
       J1=IENODE(I,1)
       J2=IENODE(I,2)
!       DIJ=0.5_SP*(DT(J1)+DT(J2))


       ELIJ=0.5_SP*(EGF(J1)+EGF(J2))
#      if defined (HEAT_FLUX)
       ELIJ=ELIJ-0.5_SP*(EGF_AIR(J1)+EGF_AIR(J2))
#      endif
       DIJ=0.5_SP*(DT(J1)*DZ(J1,K)+DT(J2)*DZ(J2,K))
#      if defined (WET_DRY)
       IF(ISWETCT(IA)*ISWETC(IA) == 1 .OR. ISWETCT(IB)*ISWETC(IB) == 1)THEN
#      endif       
       UIJ1(I)=U(IA,K)+UALFA(IA)*UIJ1(I)
       VIJ1(I)=V(IA,K)+VALFA(IA)*VIJ1(I)
       UIJ2(I)=U(IB,K)+UALFA(IB)*UIJ2(I)
       VIJ2(I)=V(IB,K)+VALFA(IB)*VIJ2(I)

!      NORMAL VELOCITY              
       UIJ=0.5_SP*(UIJ1(I)+UIJ2(I))
       VIJ=0.5_SP*(VIJ1(I)+VIJ2(I))
       EXFLUX = DIJ*(-UIJ*DLTYC(I) + VIJ*DLTXC(I))

       XADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*UIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*UIJ1(I))*0.5_SP
       YADV=EXFLUX*((1.0_SP-SIGN(1.0_SP,EXFLUX))*VIJ2(I)+(1.0_SP+SIGN(1.0_SP,EXFLUX))*VIJ1(I))*0.5_SP

       !!CALCULATE BOUNDARY FLUX AUGMENTERS
       TPA = FLOAT(1-ISBC(I))*EPOR(IA)
       TPB = FLOAT(1-ISBC(I))*EPOR(IB)

       !!ACCUMULATE ADVECTIVE + DIFFUSIVE + BAROTROPIC PRESSURE GRADIENT TERMS
!       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+FXX*TPA
!       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+FYY*TPA
!       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-FXX*TPB
!       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-FYY*TPB
       XFLUX(IA,K)=XFLUX(IA,K)+XADV*TPA+(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I)))*EPOR(IA)
       YFLUX(IA,K)=YFLUX(IA,K)+YADV*TPA+(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I)))*EPOR(IA)
       XFLUX(IB,K)=XFLUX(IB,K)-XADV*TPB-(FXX(I)+3.0_SP*FXX(I)*FLOAT(ISBC(I)))*EPOR(IB)
       YFLUX(IB,K)=YFLUX(IB,K)-YADV*TPB-(FYY(I)+3.0_SP*FYY(I)*FLOAT(ISBC(I)))*EPOR(IB)

#  if defined (WET_DRY)
    END IF
#  endif


!        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV*DT1(IA)*ELIJ*DLTYC(I)
!        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV*DT1(IA)*ELIJ*DLTXC(I)
!        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV*DT1(IB)*ELIJ*DLTYC(I)
!        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV*DT1(IB)*ELIJ*DLTXC(I)
        PSTX_TM(IA,K)=PSTX_TM(IA,K)-GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTYC(I)
        PSTY_TM(IA,K)=PSTY_TM(IA,K)+GRAV_E(IA)*DT1(IA)*DZ1(IA,K)*ELIJ*DLTXC(I)
        PSTX_TM(IB,K)=PSTX_TM(IB,K)+GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTYC(I)
        PSTY_TM(IB,K)=PSTY_TM(IB,K)-GRAV_E(IB)*DT1(IB)*DZ1(IB,K)*ELIJ*DLTXC(I)

     END DO
   END DO

   DEALLOCATE(UIJ1,VIJ1,UIJ2,VIJ2)
   DEALLOCATE(UALFA,VALFA)
   DEALLOCATE(FXX,FYY)
  

    DO I=1,N
#     if defined (WET_DRY)
#       if !defined (SEMI_IMPLICIT)
        ISWETTMP = ISWETCT(I)*ISWETC(I)
#       else
        ISWETTMP = ISWETCT(I)
#       endif
        DO K=1,KBM1
	 XFLUX(I,K)  = XFLUX(I,K)*ISWETTMP
	 YFLUX(I,K)  = YFLUX(I,K)*ISWETTMP
         PSTX_TM(I,K)= PSTX_TM(I,K)*ISWETTMP
         PSTY_TM(I,K)= PSTY_TM(I,K)*ISWETTMP
        END DO
#     endif
       DO K=1,KBM1
        XFLUX(I,K)=XFLUX(I,K)+PSTX_TM(I,K)
        YFLUX(I,K)=YFLUX(I,K)+PSTY_TM(I,K)
       END DO
   END DO

!
!-------ADD VERTICAL CONVECTIVE FLUX, CORIOLIS TERM AND BAROCLINIC PG TERM----!
!
   DO I=1,N
#    if defined (WET_DRY)
#    if !defined (SEMI_IMPLICIT)
     IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#    else
     IF(ISWETCT(I) == 1)THEN
#    endif
#    endif
     DO K=1,KBM1
       IF(K == 1) THEN
         XFLUXV=-W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV=-W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       ELSE IF(K == KBM1) THEN
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))
       ELSE
         XFLUXV= W(I,K)*(U(I,K)*DZ1(I,K-1)+U(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(U(I,K)*DZ1(I,K+1)+U(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
         YFLUXV= W(I,K)*(V(I,K)*DZ1(I,K-1)+V(I,K-1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K-1))-&
                 W(I,K+1)*(V(I,K)*DZ1(I,K+1)+V(I,K+1)*DZ1(I,K))/&
                 (DZ1(I,K)+DZ1(I,K+1))
       END IF

!       XFLUX(I,K)=XFLUX(I,K)+XFLUXV/DZ(K)*ART(I)&
!                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*ART(I)
!       YFLUX(I,K)=YFLUX(I,K)+YFLUXV/DZ(K)*ART(I)&
!                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*ART(I)
       XFLUX(I,K)=XFLUX(I,K)+XFLUXV*ART(I)&
                 +DRHOX(I,K)-COR(I)*V(I,K)*DT1(I)*DZ1(I,K)*ART(I)
       YFLUX(I,K)=YFLUX(I,K)+YFLUXV*ART(I)&
                 +DRHOY(I,K)+COR(I)*U(I,K)*DT1(I)*DZ1(I,K)*ART(I)

     END DO
#  if defined (WET_DRY)
    END IF
#  endif
   END DO

      DO I=1,N
         IF(ISBCE(I) == 2) THEN
            DO K=1,KBM1
#              if !defined (SEMI_IMPLICIT)
               XFLUX(I,K)=0.0_SP
               YFLUX(I,K)=0.0_SP
#              else
               XFLUX(I,K)=PSTX_TM(I,K)
               YFLUX(I,K)=PSTY_TM(I,K)
#              endif
            END DO
         END IF
      END DO

   !ADJUST FLUX AT RIVER INFLOWS
   IF(NUMQBC >= 1) THEN
     IF(INFLOW_TYPE == 'node') THEN
       DO II=1,NUMQBC
         J=INODEQ(II)
         I1=NBVE(J,1)
         I2=NBVE(J,NTVE(J))
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=0.5_SP*QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ(J,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           XFLUX(I2,K)=XFLUX(I2,K)-TEMP/DZ(J,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
!           YFLUX(I2,K)=YFLUX(I2,K)-TEMP/DZ(J,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           XFLUX(I2,K)=XFLUX(I2,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
           YFLUX(I2,K)=YFLUX(I2,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE IF(INFLOW_TYPE == 'edge') THEN
       DO II=1,NUMQBC
         I1=ICELLQ(II)
         DO K=1,KBM1
           VLCTYQ(II)=QDIS(II)/QAREA(II)
!           TEMP=QDIS(II)*VQDIST(II,K)*VLCTYQ(II)
           TEMP=QDIS(II)*VQDIST(II,K)*VQDIST(II,K)*VLCTYQ(II)/DZ1(I1,K)
!           XFLUX(I1,K)=XFLUX(I1,K)-TEMP/DZ1(I1,K)*COS(ANGLEQ(II))
!           YFLUX(I1,K)=YFLUX(I1,K)-TEMP/DZ1(I1,K)*SIN(ANGLEQ(II))
           XFLUX(I1,K)=XFLUX(I1,K)-TEMP*COS(ANGLEQ(II))
           YFLUX(I1,K)=YFLUX(I1,K)-TEMP*SIN(ANGLEQ(II))
         END DO
       END DO
     ELSE
       PRINT*,'INFLOW_TYPE NOT CORRECT'
       CALL PSTOP
     END IF
   END IF

   DO I=1,N
#  if defined (WET_DRY)
#   if !defined (SEMI_IMPLICIT)
    IF(ISWETCT(I)*ISWETC(I) == 1)THEN
#   else
    IF(ISWETCT(I) == 1)THEN
#   endif
#  endif

     DO K=1,KBM1

!       UF(I,K)=U(I,K)*DT1(I,K)/D1(I,K)-DTI*XFLUX(I,K)/ART(I)/D1(I,K)
!       VF(I,K)=V(I,K)*DT1(I,K)/D1(I,K)-DTI*YFLUX(I,K)/ART(I)/D1(I,K)
       UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*XFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))
       VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*YFLUX(I,K)/ART(I)/(D1(I)*DZ1(I,K))

     END DO

#  if defined (WET_DRY)
    ELSE
     DO K=1,KBM1
       UF(I,K)=0.0_SP
       VF(I,K)=0.0_SP
     END DO
    END IF
#  endif
   END DO

   RETURN
   END SUBROUTINE ADV_UV_EDGE_GCN
!==============================================================================!
