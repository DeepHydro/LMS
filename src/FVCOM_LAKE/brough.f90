!==============================================================================|
!  Calculate Bottom Drag Coefficient based on Bottom Roughness                 !
!   note:                                                                      !
!   when the log function derived from the constant stress log-viscous         !
!   layer is applied to an estuary, if the value of z0 is close to             !
!   (zz(kbm1)-z(kb)*dt1, drag coefficient "cbc" could become a huge            !
!   number due to near-zero value of alog function. In our application         !
!   we simply cutoff at cbc=0.005. One could adjust this cutoff value          !
!   based on observations or his or her experiences.                           !   
!   CALCULATES:   WUBOT(N), WVBOT(N) : BOTTOM SHEAR STRESSES                   !
!==============================================================================|

   SUBROUTINE BOTTOM_ROUGHNESS

!==============================================================================!
   USE ALL_VARS
#  if defined (WET_DRY)
   USE MOD_WD
#  endif      
   IMPLICIT NONE
   INTEGER :: I,II
   REAL(SP), PARAMETER  :: KAPPA = .40_SP   !!VON KARMAN LENGTH SCALE
   REAL(SP), PARAMETER  :: VK2   = .160_SP  !!KAPPA SQUARED
   REAL(SP)             :: CBCMIN,Z0,ZTEMP,BTPS,RR,U_TAUB,Z0B_GOTM,Z0B_TEMP

#  if defined (TWO_D_MODEL)
!  FOR 2D MODEL ONLY
   REAL(SP), PARAMETER  :: CONST_CD=.0015_SP      !! CD SET CONSTANT TO THIS VALUE
   REAL(SP), PARAMETER  :: ALFA = .166667_SP,   & !! POWER OF WATER DEPTH
                           NN   = 0.02_SP         !! FACTOR TO DIVIDE
   REAL(SP), PARAMETER  :: CFMIN   = .0025_SP,  & !! DEEP WATER VALUE
                           H_BREAK = 1.0_SP,    & !! 
                           THETA   = 10._SP,    & !!
                           LAMB    = 0.3333333333_SP
#  endif
!==============================================================================!
!
!  SET CONSTANTS
!
   CBCMIN = BFRIC
   Z0     = Z0B

!==============================================================================|
   IF(BROUGH_TYPE == 'orig')THEN !USE ORIGINAL DOTFVM FORM FOR BOTTOM FRICTION  |
!==============================================================================|
     DO I=1,N
       IF(DT1(I) > 3.0)THEN
        ZTEMP=(ZZ1(I,KBM1)-Z1(I,KB))*DT1(I)/Z0
        CBC(I) = MAX(CBCMIN,VK2/(LOG(ZTEMP))**2)
       ELSE
        ZTEMP=(ZZ1(I,KBM1)-Z1(I,KB))*3.0/Z0
        CBC(I) = MAX(CBCMIN,VK2/(LOG(ZTEMP))**2)
       END IF
     END DO
!==============================================================================|
   ELSE IF(BROUGH_TYPE == 'gotm')THEN !GOTM FORMULATION FOR BOTTOM FRICTION    |
!==============================================================================|

!----Convert Input Z0B to GOTMS H0B
!     H0B = Z0B/.03  
! DAS fixed bug to match gotm's friction.f90
     DO I=1,N
     U_TAUB = 0.
     DO II=1,40       
        IF (UMOL <= 0.) THEN
!           Z0B_GOTM= 0.03*Z0B      
           Z0B_GOTM=Z0B   !0.03*H0B 
        ELSE
           Z0B_GOTM=0.1*UMOL/MAX(UMOL,U_TAUB)+Z0B !0.03*H0B
        END IF
          ztemp=(zz1(I,kbm1)-z1(I,kb))*dt1(i)
        RR=KAPPA/(LOG((Z0B_GOTM+ZTEMP)/Z0B_GOTM))
      U_TAUB = RR*SQRT( U(I,KBM1)*U(I,KBM1) + V(I,KBM1)*V(I,KBM1) )
     END DO
     CBC(I) =   RR*RR
     END DO

!==============================================================================|
   ELSE IF(BROUGH_TYPE == 'user_defined')THEN !Use User Defined broud_ud.F     | 
!==============================================================================|
   
     CALL BROUGH_UD

   END IF

!==============================================================================|
!  CALCULATE SHEAR STRESS ON BOTTOM  --> WUBOT/WVBOT                           |
!==============================================================================|
   DO  I = 1, N
     IF(D1(I) > 0.0_SP) THEN    
       BTPS= CBC(I)*SQRT(U(I,KBM1)**2+V(I,KBM1)**2)
       WUBOT(I) = -BTPS * U(I,KBM1)
       WVBOT(I) = -BTPS * V(I,KBM1)    
     ELSE
       WUBOT(I) = 0.0_SP
       WVBOT(I) = 0.0_SP
     END IF
   END DO
   

   RETURN
   END SUBROUTINE BOTTOM_ROUGHNESS
!==============================================================================|
