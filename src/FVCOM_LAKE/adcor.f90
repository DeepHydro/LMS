   SUBROUTINE ADCOR

   USE ALL_VARS

#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif

   IMPLICIT NONE
   REAL(SP) :: UFC(0:NT,KB),VFC(0:NT,KB)
   REAL(SP),PARAMETER :: BETA0=0.5
   REAL(SP) ::CURCOR,PRECOR
   INTEGER :: I,K  
   REAL(SP) :: U_TMP,V_TMP,UF_TMP,VF_TMP  


   UFC=0.0
   VFC=0.0


#  if !defined (TWO_D_MODEL)
   DO I=1,N
#  if defined (SPHERICAL) && (NORTHPOLE)
     IF(CELL_NORTHAREA(I) == 1)THEN
       DO K=1,KBM1
#        if !defined (SEMI_IMPLICIT)
         U_TMP = -V(I,K)*COS(XC(I)*DEG2RAD)-U(I,K)*SIN(XC(I)*DEG2RAD)
         V_TMP = -V(I,K)*SIN(XC(I)*DEG2RAD)+U(I,K)*COS(XC(I)*DEG2RAD)

         UF_TMP=U_TMP*DT1(I)/D1(I)-DTI*UFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF_TMP=V_TMP*DT1(I)/D1(I)-DTI*VFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))

         UF(I,K)  = VF_TMP*COS(XC(I)*DEG2RAD)-UF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = UF_TMP*COS(XC(I)*DEG2RAD)+VF_TMP*SIN(XC(I)*DEG2RAD)
         VF(I,K)  = -VF(I,K)    
#        else
         XFLUX3(I,K) = UFC(I,K)
         YFLUX3(I,K) = VFC(I,K)
         XFLUX3_NP(I,K) = UFC_NP(I,K)
         YFLUX3_NP(I,K) = VFC_NP(I,K)
#        endif
       END DO
     ELSE
#  endif   
       DO K=1,KBM1
#        if !defined (SEMI_IMPLICIT)
         UF(I,K)=U(I,K)*DT1(I)/D1(I)-DTI*UFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
         VF(I,K)=V(I,K)*DT1(I)/D1(I)-DTI*VFC(I,K)/ART(I)/(D1(I)*DZ1(I,K))
#        else
         XFLUX3(I,K) = UFC(I,K)
         YFLUX3(I,K) = VFC(I,K) 
#        endif
       END DO
#  if defined (SPHERICAL) && (NORTHPOLE)
     END IF  
#  endif        
   END DO
#  endif

   RETURN
   END SUBROUTINE ADCOR
