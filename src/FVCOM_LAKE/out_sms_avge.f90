!==============================================================================|
!   Write Data Averages over Interval INT_AVGE To Assess if Model has          |
!   Achieved Quasi-Periodic Behavior                                           |
!==============================================================================|

   SUBROUTINE out_sms_avge

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR 
#  endif
#  if defined (BALANCE_2D)
   USE MOD_BALANCE_2D
#  endif

   IMPLICIT NONE
   INTEGER :: I,K
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,WWTMP,KMTMP,KHTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1TMP,S1TMP,R1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP

#  if defined (BALANCE_2D)
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ADVUA2TMP,  ADVVA2TMP,  ADFX2TMP, ADFY2TMP  
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: DRX2D2TMP,  DRY2D2TMP,  CORX2TMP, CORY2TMP  
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: PSTX2TMP,   PSTY2TMP,   ADX2D2TMP,ADY2D2TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: WUSURBF2TMP,WVSURBF2TMP,DUDT2TMP, DVDT2TMP
#  endif

   CHARACTER(LEN=4)   :: FILENUMBER
   CHARACTER(LEN=120) :: DIR 
   CHARACTER(LEN=120) :: FNAME
   INTEGER  :: J1,J2,END_AVGE
   REAL(SP) :: TMP,FAC
    REAL(SP) ::uvwaver,uaver,vaver,waver
    REAL(SP) ::taver
!==============================================================================|
   
   END_AVGE = BEG_AVGE + NUM_AVGE*INT_AVGE - 1
 
!------------------------------------------------------------------------------!
!  ALLOCATE DATA FOR STORING AVERAGES                                          !
!------------------------------------------------------------------------------!

   IF(IINT == BEG_AVGE)THEN
     ALLOCATE(U_AVE(0:NT,KB))  ; U_AVE  = 0.0_SP
     ALLOCATE(V_AVE(0:NT,KB))  ; V_AVE  = 0.0_SP
     ALLOCATE(W_AVE(0:NT,KB))  ; W_AVE  = 0.0_SP
     ALLOCATE(KM_AVE(0:MT,KB)) ; KM_AVE = 0.0_SP
     ALLOCATE(KH_AVE(0:MT,KB)) ; KH_AVE = 0.0_SP
     ALLOCATE(S_AVE(0:MT,KB))  ; S_AVE  = 0.0_SP
     ALLOCATE(T_AVE(0:MT,KB))  ; T_AVE  = 0.0_SP
     ALLOCATE(R_AVE(0:MT,KB))  ; R_AVE  = 0.0_SP
     ALLOCATE(EL_AVE(0:MT))    ; EL_AVE = 0.0_SP

#  if defined (BALANCE_2D)
     ALLOCATE(ADVUA2_AVE(0:NT))   ; ADVUA2_AVE   = 0.0_SP
     ALLOCATE(ADVVA2_AVE(0:NT))   ; ADVVA2_AVE   = 0.0_SP
     ALLOCATE(ADFX2_AVE(0:NT))    ; ADFX2_AVE    = 0.0_SP  
     ALLOCATE(ADFY2_AVE(0:NT))    ; ADFY2_AVE    = 0.0_SP  
     ALLOCATE(DRX2D2_AVE(0:NT))   ; DRX2D2_AVE   = 0.0_SP
     ALLOCATE(DRY2D2_AVE(0:NT))   ; DRY2D2_AVE   = 0.0_SP
     ALLOCATE(CORX2_AVE(0:NT))    ; CORX2_AVE    = 0.0_SP  
     ALLOCATE(CORY2_AVE(0:NT))    ; CORY2_AVE    = 0.0_SP  
     ALLOCATE(PSTX2_AVE(0:NT))    ; PSTX2_AVE    = 0.0_SP
     ALLOCATE(PSTY2_AVE(0:NT))    ; PSTY2_AVE    = 0.0_SP
     ALLOCATE(ADX2D2_AVE(0:NT))   ; ADX2D2_AVE   = 0.0_SP
     ALLOCATE(ADY2D2_AVE(0:NT))   ; ADY2D2_AVE   = 0.0_SP
     ALLOCATE(WUSURBF2_AVE(0:NT)) ; WUSURBF2_AVE = 0.0_SP 
     ALLOCATE(WVSURBF2_AVE(0:NT)) ; WVSURBF2_AVE = 0.0_SP 
     ALLOCATE(DUDT2_AVE(0:NT))    ; DUDT2_AVE    = 0.0_SP    
     ALLOCATE(DVDT2_AVE(0:NT))    ; DVDT2_AVE    = 0.0_SP    
#  endif

   END IF

!------------------------------------------------------------------------------!
!  UPDATE AVERAGES                                                             !
!------------------------------------------------------------------------------!

   TMP = INT_AVGE
   FAC = 1.0_SP/TMP
   IF(IINT >= BEG_AVGE .AND. IINT <= END_AVGE)THEN 
     U_AVE  = U_AVE  + U*FAC
     V_AVE  = V_AVE  + V*FAC
     W_AVE  = W_AVE  + WW*FAC
     KM_AVE = KM_AVE + KM*FAC
     KH_AVE = KH_AVE + KH*FAC
     S_AVE  = S_AVE  + S1*FAC
     T_AVE  = T_AVE  + T1*FAC
     R_AVE  = R_AVE  + RHO1*FAC
     EL_AVE = EL_AVE + EL*FAC

#  if defined (BALANCE_2D)
     ADVUA2_AVE   = ADVUA2_AVE   + ADVUA2*FAC
     ADVVA2_AVE   = ADVVA2_AVE   + ADVVA2*FAC
     ADFX2_AVE    = ADFX2_AVE    + ADFX2*FAC
     ADFY2_AVE    = ADFY2_AVE    + ADFY2*FAC
     DRX2D2_AVE   = DRX2D2_AVE   + DRX2D2*FAC
     DRY2D2_AVE   = DRY2D2_AVE   + DRY2D2*FAC
     CORX2_AVE    = CORX2_AVE    + CORX2*FAC
     CORY2_AVE    = CORY2_AVE    + CORY2*FAC
     PSTX2_AVE    = PSTX2_AVE    + PSTX2*FAC
     PSTY2_AVE    = PSTY2_AVE    + PSTY2*FAC
     ADX2D2_AVE   = ADX2D2_AVE   + ADX2D2*FAC
     ADY2D2_AVE   = ADY2D2_AVE   + ADY2D2*FAC
     WUSURBF2_AVE = WUSURBF2_AVE + WUSURBF2*FAC
     WVSURBF2_AVE = WVSURBF2_AVE + WVSURBF2*FAC
     DUDT2_AVE    = DUDT2_AVE    + DUDT2*FAC
     DVDT2_AVE    = DVDT2_AVE    + DVDT2*FAC
#  endif

   END IF

!------------------------------------------------------------------------------!
!  OPEN FILE (Name Based on Iteration Number)                                  !
!------------------------------------------------------------------------------!

   J1 = MOD((IINT+1-BEG_AVGE),INT_AVGE)
   J2 = (IINT+1-BEG_AVGE)/INT_AVGE

   IF(IINT >= BEG_AVGE .AND.  J1 == 0 .AND. IINT <= END_AVGE)THEN

   
   IF(MSR)THEN
     WRITE(FILENUMBER,'(I4.4)') J2 
     OPEN(IOSMSV,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//'_avge'//FILENUMBER//'_uv.xy',STATUS='unknown')
     OPEN(IOSMST,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//'_avge'//FILENUMBER//'_elts.xy',STATUS='unknown')
     REWIND(IOSMSV)
     REWIND(IOSMST)

     REWIND(1)
     WRITE(IPT,*)'DUMPING SMS AVGES FILE: ',TRIM(FNAME)
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Single Processor Case)                                !
!------------------------------------------------------------------------------!

   IF(SERIAL)THEN

     WRITE(IOSMSV,10)
     WRITE(IOSMSV,20) N,10
     !! ELEMENT BASED VALUES

     DO I=1,N
       uaver = SUM(U_AVE(I,:))/KBM1
       vaver = SUM(V_AVE(I,:))/KBM1
       waver = SUM(W_AVE(I,:))/KBM1
       uvwaver=  SQRT (uaver*uaver+vaver*vaver+waver*waver)
       WRITE(IOSMSV,'(12E17.8)') VX(I)+VXMIN,VY(I)+VYMIN, U_AVE(I,1),V_AVE(I,1),W_AVE(I,1),&
        U_AVE(I,KBM1),V_AVE(I,KBM1),W_AVE(I,KBM1),uaver,vaver,waver,uvwaver
     END DO
     WRITE(IOSMSV,*) 'IINT==',IINT

   WRITE(IOSMSV,10)
   WRITE(IOSMSV,30) M,4
     DO I=1,M
       taver=  SUM(T_AVE(I,:))/KBM1
       WRITE(IOSMST,'(6E17.8)')VX(I)+VXMIN,VY(I)+VYMIN,BaseSurfaceElevation+EL_AVE(I),T_AVE(I,1),T_AVE(I,KBM1),taver
     END DO
    WRITE(IOSMST,*) 'IINT==',IINT

     
#  if defined (BALANCE_2D)
     
     DO I=1,N
       WRITE(1) ADVUA2_AVE(I),   ADVVA2_AVE(I),  ADFX2_AVE(I),  ADFY2_AVE(I),&
                DRX2D2_AVE(I),   DRY2D2_AVE(I),  CORX2_AVE(I),  CORY2_AVE(I),&
		PSTX2_AVE(I),    PSTY2_AVE(I),   ADX2D2_AVE(I), ADY2D2_AVE(I),&
		WUSURBF2_AVE(I), WVSURBF2_AVE(I),DUDT2_AVE(I),  DVDT2_AVE(I)
     ENDDO
#   endif

   END IF

   IF(MSR)CLOSE(IOSMSV)
   IF(MSR)CLOSE(IOSMST)
!------------------------------------------------------------------------------!
!  REINITIALIZE AVERAGING ARRAYS                                               !
!------------------------------------------------------------------------------!

   U_AVE  = 0.0_SP
   V_AVE  = 0.0_SP
   W_AVE  = 0.0_SP
   KM_AVE = 0.0_SP
   KH_AVE = 0.0_SP
   S_AVE  = 0.0_SP
   T_AVE  = 0.0_SP
   R_AVE  = 0.0_SP
   EL_AVE = 0.0_SP
   
#  if defined (BALANCE_2D) 
     ADVUA2_AVE   = 0.0_SP
     ADVVA2_AVE   = 0.0_SP
     ADFX2_AVE    = 0.0_SP  
     ADFY2_AVE    = 0.0_SP  
     DRX2D2_AVE   = 0.0_SP
     DRY2D2_AVE   = 0.0_SP
     CORX2_AVE    = 0.0_SP  
     CORY2_AVE    = 0.0_SP  
     PSTX2_AVE    = 0.0_SP
     PSTY2_AVE    = 0.0_SP
     ADX2D2_AVE   = 0.0_SP
     ADY2D2_AVE   = 0.0_SP
     WUSURBF2_AVE = 0.0_SP 
     WVSURBF2_AVE = 0.0_SP 
     DUDT2_AVE    = 0.0_SP    
     DVDT2_AVE    = 0.0_SP    
#   endif

   END IF

10 FORMAT('scat2d')
20 FORMAT('xyd ',I10,' suva ',I2,' su sv sw bu bv bw ua va wa uvwa')
30 FORMAT('xyd ',I10,' elts ',I2,' el st bt avt')

   RETURN
   END SUBROUTINE out_sms_avge   
!==============================================================================|
