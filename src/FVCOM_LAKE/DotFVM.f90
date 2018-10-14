!==============================================================================!
!  VERSION 1.0.1  
!==============================================================================!

!Preprocessor GCN;WET_DRY;DOS;WATER_QUALITY;HEAT_FLUX;ICE;MPDATA

   PROGRAM DotFVM

!==============================================================================!
!  INCLUDE MODULES                                                             !
!==============================================================================!

   USE ALL_VARS
   USE MOD_CLOCK
   USE MOD_WQM

#  if defined (WET_DRY)
   USE MOD_WD
#  endif

   USE BCS
!   USE MOD_LAG
#  if defined (GOTM)
   USE MOD_GOTM
#  endif

#  if defined (SEDIMENT)
   USE MOD_SED
#  endif

#  if defined (BioGen)
   USE MOD_BIO_3D
#  endif

   USE MOD_OBCS

#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif 

   USE MOD_HEATFLUX

#  if defined (SEMI_IMPLICIT)
   USE MOD_SEMI_IMPLICIT
#  endif

!  USE Mod_DotFVMProject
!------------------------------------------------------------------------------|
   IMPLICIT NONE

   REAL(SP) :: TMP1,TMP,UTMP,VTMP,TTIME
   INTEGER  :: I,K,J,IERR,N1,J1,J2,JN
   REAL(SP), ALLOCATABLE :: FTEMP(:),FTEMP2(:)
   CHARACTER(LEN=13) :: TSTRING

   DotFVM_VERSION     = 'DOTFVM 1.0'
   DotFVM_WEBSITE     = 'http://ww.hust.edu.cn'

!==============================================================================!
!   SETUP PARALLEL ENVIRONMENT                                                 !
!==============================================================================!

   SERIAL = .TRUE. 
   PAR    = .FALSE. 
   MSR    = .TRUE.
   MYID   = 1
   NPROCS = 1

     CALL GET_CASENAME(CASENAME)
     CALL DATA_RUN
  ! CALL DATA_RUN_HFX
!==============================================================================!
!   SETUP MODEL RUN                                                            !
!==============================================================================!

!  READ PARAMETERS CONTROLLING MODEL RUN
!    call OpenProject
!==============================================================================!
!  Delete existing log file                                                            !
!==============================================================================!
    IOLOG=97
   OPEN(IOLOG,FILE=TRIM(OUTDIR)//'_report.txt',STATUS='UNKNOWN')
   CLOSE(IOLOG,STATUS='DELETE')
   OPEN(IOLOG,FILE=TRIM(OUTDIR)//'_report.txt',STATUS='UNKNOWN')

!  READ PARAMETERS CONTROLLING DYE RELEASE
#  if defined (DYE_RELEASE)
   CALL SET_DYE_PARAM
#  endif   

!  OPEN INPUT/OUTPUT FILES
    CALL IOFILES
    CALL IOFILES_HFX

!  DETERMINE NUMBER OF ELEMENTS AND NODES IN THE MODEL
   CALL GETDIM
 
    
!  INPUT AND SETUP BOUNDARY FORCING (HEAT/RIVERS/WIND/etc)
   CALL BCS_FORCE
   CALL BCS_FORCE_HFX

!  INPUT WATER QUALITY MODEL VARIABLES
   IF(WQM_ON) then
       call WQMPARA
        CALL BCS_FORCE_WQM
    endif

!  ALLOCATE FLOWFIELD VARIABLES
   CALL ALLOC_VARS
   
!  ALLOCATE ELM1 AND ELM2 FOR ORLANSKI RADIATION OPEN BOUNDARY CONDITION
   CALL ALLOC_OBC_DATA

!  ALLOCATE AND WET/DRY CONTROL ARRAYS
#  if defined (WET_DRY)
   CALL ALLOC_WD_DATA
#  endif

!  ALLOCATE WATER QUALITY MODEL VARIABLES
   IF(WQM_ON)CALL ALLOC_WQM_VARS


!  ALLOCATE DYE VARS
#  if defined (DYE_RELEASE)
   CALL ALLOC_VARS_DYE
#  endif

#  if defined (BioGen)
   KBV=KB
   CALL GET_PARAMETER
   POINT_ST_TYPE='NONE' !No bio river yet
#  endif

!  SHIFT GRID/CORIOLIS/BATHYMETRY TO LOCAL DOMAIN
   CALL PDOMDEC

   
!  SET UP GRID METRICS (FLUX EDGES/CONTROL VOLUMES/ETC)
   CALL TRIANGLE_GRID_EDGE      !Set up fluxes and control Volumes
   CALL SET_SIGMA(INDEX_VERCOR) !Build Vertical Coordinate
   CALL CELL_AREA               !Calculate Element and Control Volume Areas
   CALL SHAPE_COEF_GCN          !Calc Shape Coefficients for Flux Construction


   CALL SET_BNDRY               !Boundary Condition Metrics

#  if defined (SEMI_IMPLICIT)
   CALL SET_IMPLICIT_PARAM
   CALL ALLOC_VARS_SEMI
#  endif
!
!  INITIALIZE FLOWFIELD  --> [T,S,U,V,EL,D,Q2,Q2L]
   CALL STARTUP

!  INITIALIZE GOTM
#  if defined (GOTM)
   CALL INIT_GOTM
#  endif
!
!  INITIALIZE SEDIMENT MODEL
#  if defined (SEDIMENT)
   IF(SEDIMENT_ON)CALL SETUP_SED(RESTART_SED)
#  endif


!  CALCULATE DEPTH HORIZONTAL DERIVATIVES
   CALL DEPTH_GRADIENT

# if defined (BioGen)
  IF(RESTART == 'hot_start')  CALL BIO_HOT_START
  IF(RESTART == 'cold_start') CALL BIO_INITIAL
!        CALL BIO_INITIAL
# endif

!
!  CALCULATE INTERVALS FOR SST DATA ASSIMILATION 
! IF(C_HFX)CALL BCOND_HFX

!  REPORT STATISTICS ON INITIAL VALUES
   !CALL REPORT('INITIAL VALUE INFORMATION')
   CALL REPORT('初始化信息')

!
!  CALCULATE INTERNAL TIME STEP AND SET INTEGRATION LIMITS
   ISTART=IINT+1

   CALL START_CLOCK

!  REPORT INTEGRATION INITIAL TIMES 
!
   IF(MSR)CALL REPORT_SIMTIME

  call out_binary_cfd(1,0)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  BEGIN MAIN LOOP OVER PHYSICAL TIME STEPS                                    |
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   DO IINT=ISTART,IEND
     TIME =DTI*FLOAT(IINT)/86400.0_SP
     THOUR=DTI*FLOAT(IINT)/3600.0_SP

!----SET RAMP FACTOR TO EASE SPINUP--------------------------------------------!
     RAMP=1.0_SP
     IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT)/FLOAT(IRAMP))
!     IF(IRAMP /= 0) RAMP=TANH(FLOAT(IINT-ISTART+1)/FLOAT(IRAMP))

!----SET UP WATER QUALITY MODEL COEFFICIENTS-----------------------------------!
     IF(WQM_ON)THEN
        TIME_R=MOD(IINT*DTI/3600.0_SP-14.743_SP, 24.0_SP)+6.0_SP
       CALL WQMCONST
     END IF

!----ADJUST CONSISTENCY BETWEEN 3-D VELOCITY AND VERT-AVERAGED VELOCITIES------!
     CALL ADJUST2D3D(1)

!----SPECIFY THE SOLID BOUNDARY CONDITIONS OF U&V INTERNAL MODES---------------!
     CALL BCOND_GCN(5,0)


!----SPECIFY THE SURFACE FORCING OF INTERNAL MODES-----------------------------!
     CALL BCOND_GCN(8,0)
     IF(C_HFX)CALL BCOND_NetHFX

!----SPECIFY BOTTOM FRESH WATER INPUT BOUNDARY CONDITION-----------------------!
     CALL BCOND_BFW(1)

!----SPECIFY THE BOTTOM ROUGHNESS AND CALCULATE THE BOTTOM STRESSES------------!
     CALL BOTTOM_ROUGHNESS

     IF(C_HFX)CALL BCOND_NetHFX

!==============================================================================!
!  CALCULATE DISPERSION (GX/GY) AND BAROCLINIC PRESSURE GRADIENT TERMS         !
!==============================================================================!
     CALL ADVECTION_EDGE_GCN(ADVX,ADVY)          !Calculate 3-D Adv/Diff       !


     IF(.NOT. BAROTROPIC)THEN                    !Barotropic Flow ?            !  
         !Calculate the rho mean
         IF(IRHO_MEAN > 0) THEN
           IF(MOD(IINT,IRHO_MEAN) == 0) CALL RHO_MEAN    
         END IF
         IF(C_BAROPG == 'sigma')    CALL BAROPG      !Sigma Level Pressure Gradient!
         IF(C_BAROPG == 's_levels') CALL PHY_BAROPG  !Z Level Pressure Gradient    !
     END IF                                      !                             !
                                                 !                             !
     ADX2D = 0.0_SP ; ADY2D = 0.0_SP             !Initialize GX/GY Terms       !
     DRX2D = 0.0_SP ; DRY2D = 0.0_SP             !Initialize BCPG for Ext Mode !

     DO K=1,KBM1
       DO I=1, N
          ADX2D(I)=ADX2D(I)+ADVX(I,K)   !*DZ1(I,K)
          ADY2D(I)=ADY2D(I)+ADVY(I,K)   !*DZ1(I,K)
          DRX2D(I)=DRX2D(I)+DRHOX(I,K)  !*DZ1(I,K)
          DRY2D(I)=DRY2D(I)+DRHOY(I,K)  !*DZ1(I,K)
        END DO
     END DO


     CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff
     ADX2D = ADX2D - ADVUA                       !Subtract to Form GX
     ADY2D = ADY2D - ADVVA                       !Subtract to Form GY

!----INITIALIZE ARRAYS USED TO CALCULATE AVERAGE UA/E  OVER EXTERNAL STEPS-----!
     UARD = 0.0_SP
     VARD = 0.0_SP
     EGF  = 0.0_SP

     EGF_AIR = 0.0_SP

!!#    if defined (WET_DRY)
!!     UARDS = 0.0_SP
!!     VARDS = 0.0_SP
!!#    endif


     IF(IOBCN > 0) THEN
         UARD_OBCN(1:IOBCN)=0.0_SP
     END IF


!==============================================================================!
!  LOOP OVER EXTERNAL TIME STEPS                                               !
!==============================================================================!
     DO IEXT=1,ISPLIT

       !!=======================================================

       TIME  =(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/86400.0_SP
       THOUR1=(DTI*FLOAT(IINT-1)+DTE*FLOAT(IEXT))/3600.0_SP

!
!----- USE RAMP VARIABLE TO EASE MODEL STARTUP---------------------------------!
!
!       TMP1 = FLOAT(IINT-ISTART)+FLOAT(IEXT)/FLOAT(ISPLIT)
       TMP1 = FLOAT(IINT-1)+FLOAT(IEXT)/FLOAT(ISPLIT)
       RAMP = 1.0_SP
       IF(IRAMP /= 0) RAMP = TANH(TMP1/FLOAT(IRAMP))
!
!------SURFACE BOUNDARY CONDITIONS FOR EXTERNAL MODEL--------------------------!
!
       CALL BCOND_GCN(9,0)

       CALL BCOND_BFW(2)
!
!------SAVE VALUES FROM CURRENT TIME STEP--------------------------------------!
!
       ELRK1 = EL1
       ELRK  = EL
       UARK  = UA
       VARK  = VA

       ELRK_AIR = EL_AIR

!------BEGIN MAIN LOOP OVER EXTERNAL MODEL 4 STAGE RUNGE-KUTTA INTEGRATION-----!
       DO K=1,4
         TIMERK = TIME + (ALPHA_RK(K)-1.)*DTE/86400.0_SP

         !FREE SURFACE AMPLITUDE UPDATE  --> ELF
         CALL EXTEL_EDGE(K)

        CALL BCOND_PA_AIR 
        ELF_AIR = ELRK_AIR +ALPHA_RK(K)*(ELF_AIR-ELRK_AIR) 


         CALL BCOND_GCN(1,K)

         DO I=1,IBCN(1)
           JN = OBC_LST(1,I)
           J=I_OBC_N(JN)
           ELF(J)=ELRK(J)+ALPHA_RK(K)*(ELF(J)-ELRK(J))
         END DO

         CALL N2E2D(ELF,ELF1)


#        if defined (WET_DRY)
         IF(WET_DRY_ON)CALL WET_JUDGE
#        endif

         CALL FLUX_OBN(K)

         !CALCULATE ADVECTIVE, DIFFUSIVE, AND BAROCLINIC MODES --> UAF ,VAF
         CALL ADVAVE_EDGE_GCN(ADVUA,ADVVA)           !Compute Ext Mode Adv/Diff

         CALL EXTUV_EDGE(K)

         CALL BCOND_GCN(2,K)

         !UPDATE WATER SURFACE ELEVATION
         CALL ASSIGN_ELM1_TO_ELM2

         EL  = ELF
         EL1 = ELF1

         EL_AIR = ELF_AIR    

         !!INTERPOLATE DEPTH FROM NODE-BASED TO ELEMENT-BASED VALUES
         CALL N2E2D(EL,EL1)

         !UPDATE DEPTH AND VERTICALLY AVERAGED VELOCITY FIELD
         D   = H + EL
         D1  = H1 + EL1
         UA  = UAF
         VA  = VAF
         DTFA = D

         !!ENSURE ALL CELLS ARE WET IN NO FLOOD/DRY CASE  
         CALL DEPTH_CHECK

         !SAVE VALUES FOR 3D MOMENTUM CORRECTION AND UPDATE
         IF(K == 3)THEN
           UARD = UARD + UA*D1
           VARD = VARD + VA*D1
           EGF  = EGF  + EL/ISPLIT
           EGF_AIR = EGF_AIR + EL_AIR/ISPLIT
         END IF

         !CALCULATE VALUES USED FOR SALINITY/TEMP BOUNDARY CONDITIONS
         IF(K == 4.AND.IOBCN > 0) THEN
           DO I=1,IOBCN
             J=I_OBC_N(I)
             TMP=-(ELF(J)-ELRK(J))*ART1(J)/DTE-XFLUX_OBCN(I)
             UARD_OBCN(I)=UARD_OBCN(I)+TMP/FLOAT(ISPLIT)
            END DO
         END IF

          !UPDATE WET/DRY FACTORS
#        if defined (WET_DRY)
         IF(WET_DRY_ON)CALL WD_UPDATE(1)
#        endif

       END DO     !! END RUNGE-KUTTA LOOP

     END DO     !! EXTERNAL MODE LOOP
!==============================================================================!
!  END LOOP OVER EXTERN STEPS                                                  !
!==============================================================================!


# if defined (SEMI_IMPLICIT)
        CALL BCOND_PA_AIR 
        CALL BCOND_GCN(1,0)
# endif

!==============================================================================!
!    ADJUST INTERNAL VELOCITY FIELD TO CORRESPOND TO EXTERNAL                  !
!==============================================================================!
     CALL ADJUST2D3D(2)

!==============================================================================!
!     CALCULATE INTERNAL VELOCITY FLUXES                                       |
!==============================================================================!

     CALL VERTVL_EDGE     ! Calculate/Update Sigma Vertical Velocity (Omega)   !

#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL WD_UPDATE(2)
#    endif

    CALL ADV_UV_EDGE_GCN ! Horizontal Advect/Diff + Vertical Advection        !
     CALL VDIF_UV         ! Implicit Integration of Vertical Diffusion of U/V  !

#    if defined (WET_DRY)
     DO I=1,N
       IF(H1(I) <= DJUST ) THEN
         DO K=1,KBM1
           UF(I,K)=UA(I)
           VF(I,K)=VA(I)
         END DO
       END IF
     END DO
#    endif

    CALL BCOND_GCN(3,0)    ! Boundary Condition on U/V At River Input           !          
     CALL WREAL           ! Calculate True Vertical Velocity (W)               !                                          !
     CALL VISCOF_H        ! Calculate horizontal diffusion coefficient for the scalar                                         !
!==============================================================================!
!    TURBULENCE MODEL SECTION                                                  |
!==============================================================================!
     IF(VERTMIX == 'closure')THEN
     !=================General Ocean Turbulence Model==========================!
#    if defined (GOTM)
        ! There is a bug in the advection of turbulent kinetic energy
        ! when using the GOTM MODULE
     CALL ADVANCE_GOTM            !!Solve TKE/TEPS eqns for KH/KM/KQ
#    else     
     !===================Original DOTFVM MY-2.5/Galperin 1988 Model=============!

     CALL ADV_Q(Q2,Q2F)       !!Advection of Q2 
     CALL ADV_Q(Q2L,Q2LF) 

     IF(TS_FCT) CALL FCT_Q2                              !Conservation Correction   !
     IF(TS_FCT) CALL FCT_Q2L                             !Conservation Correction   !

!end !defined (ONE_D_MODEL)
     CALL VDIF_Q                  !! Solve Q2,Q2*L eqns for KH/KM/KQ 
     Q2  = Q2F
     Q2L = Q2LF

#    endif
     ELSE
       KM = UMOL
       KH = UMOL*VPRNU
     END IF  

    CALL N2E3D(KM,KM1)
!==============================================================================!
!    SEDIMENT MODEL SECTION                                                    |
!==============================================================================!
#    if defined (SEDIMENT)
     ALLOCATE(FTEMP(0:MT),FTEMP2(0:NT))
     FTEMP2 = SQRT(WUBOT**2 + WVBOT**2)
     CALL E2N2D(FTEMP2,FTEMP)
     IF(SEDIMENT_ON)CALL ADVANCE_SED(DTI,THOUR*3600,FTEMP)
     DEALLOCATE(FTEMP,FTEMP2)
#    endif

#    if defined (BioGen)
     CALL BIO_3D1D                          
#    endif
!==============================================================================!
!    UPDATE TEMPERATURE IN NON-BAROTROPIC CASE                                 !
!==============================================================================!
     IF(TEMP_ON)THEN 
     
     CALL ADV_T                                     !Advection                 !

     IF(TS_FCT) CALL FCT_T                             !Conservation Correction   !

     IF(CASENAME == 'gom')THEN
       CALL VDIF_TS_GOM(1,TF1)
     ELSE  
       CALL VDIF_TS(1,TF1)                            !Vertical Diffusion        !
     END IF

     CALL BCOND_TS(1)                               !Boundary Conditions       !

     T1 = TF1                                       !Update to new time level  !
     CALL N2E3D(T1,T)                               !Shift to Elements         !

     END IF                                         !                          !
!==============================================================================!
!    UPDATE SALINITY IN NON-BAROTROPIC CASE                                    !
!==============================================================================!
     IF(SALINITY_ON)THEN                            !                          !   
     CALL ADV_S                                     !Advection                 !

     IF(TS_FCT) CALL FCT_S                             !Conservation Correction   !


     CALL VDIF_TS(2,SF1)                            !Vertical Diffusion        !
     CALL BCOND_TS(2)                               !Boundary Conditions       !


     S1 = SF1                                       !Update to new time level  !
     CALL N2E3D(S1,S)                               !Shift to Elements         !

     END IF                                         !                          !
!==============================================================================!
#  if defined (DYE_RELEASE)

!==============================================================================!
!    UPDATE DYE IN NON-BAROTROPIC CASE                                         !
!==============================================================================!
!     IF(DYE_ON.AND.IINT.GE.IINT_SPE_DYE_B) THEN     !                          !                          !
     IF(DYE_ON) THEN     !                          !                          !
     CALL ADV_DYE                                   !Advection                 !

!     IF(TS_FCT) CALL AVER_DYE                      !Conservation Correction   !

     CALL VDIF_DYE(DYEF)                            !Vertical Diffusion        !

!check
!!     DYE = DYEF                                     !Update to new time level  !
!!     IF(IINT.GE.IINT_SPE_DYE_B) CALL ARCHIVE
!check    
     CALL BCOND_DYE                                 !Boundary Conditions       !
     DYE = DYEF                                     !Update to new time level  !

!!     IF(MSR) WRITE(IPT,*) 'CALL Dye_on--iint=',iint,IINT_SPE_DYE_B
     END IF                                         !                          !
!==================================================================================!
#  endif

#    if !defined (MPDATA)
     IF(POINT_ST_TYPE == 'calculated')THEN
!    ADJUST TEMPERATURE AND SALINITY AT RIVER MOUTHS
       CALL ADJUST_TS
     END IF  
#    endif   

!==============================================================================!
!    CALCULATE WATER QUALITY VARIABLES CONCENTRATIONS                          |
!==============================================================================!
                                                    !                          !
     IF(WQM_ON)THEN                                 !Water Quality Active?     !
         !CALL ADV_WQM                                   !Advection                 !
         CALL   ADV_WASP
         CALL FCT_WASP
         !CALL VDIF_WQM(WQM_F)                           !Vertical Diffusion        !
         CALL VDIF_WASP(WQM_F)       
         CALL EXCHANGE_WQM                              !Interprocessor Exchange   !
         !CALL BCOND_WQM                                 !Boundary Conditions       !
         CALL BCOND_WASP
         !WQM(1:M,1:KBM2,1:NB) = WQM_F(1:M,1:KBM2,1:NB)  !Update                    !
         WQM=WQM_F
     END IF                                         !                          !
!==============================================================================!

!==============================================================================!
!     UPDATE THE DENSITY IN NON-BAROTROPIC CASE                                |
!==============================================================================!
     IF(.NOT.BAROTROPIC)THEN
       IF(CTRL_DEN == 'pdensity'   ) CALL DENS
       IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
       IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
     END IF
!==============================================================================!
!     MIMIC CONVECTIVE OVERTURNING TO STABILIZE VERTICAL DENSITY PROFILE       |
!==============================================================================!

     IF(VERT_STAB)THEN
       CALL CONV_OVER
       IF(.NOT.BAROTROPIC)THEN
         IF(CTRL_DEN == 'pdensity'   ) CALL DENS
         IF(CTRL_DEN == 'sigma-t'    ) CALL DENS2
         IF(CTRL_DEN == 'sigma-t_stp') CALL DENS3
       END IF
     END IF  

!==============================================================================!
!    LAGRANGIAN PARTICLE TRACKING                                              |
!==============================================================================!

!==============================================================================!
!     UPDATE VELOCITY FIELD (NEEDED TO WAIT FOR SALINITY/TEMP/TURB/TRACER)     |
!==============================================================================!
     U = UF
     V = VF

!----SHIFT SEA SURFACE ELEVATION AND DEPTH TO CURRENT TIME LEVEL---------------!
!
     ET  = EL  
     DT  = D 
     ET1 = EL1
     DT1 = D1
     
#    if defined (WET_DRY)
     IF(WET_DRY_ON) CALL WD_UPDATE(3)
#    endif

!==============================================================================!
!    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
!==============================================================================!

     IF(MSR)CALL REPORT_TIME(IINT,ISTART,IEND,TIME*86400,IPT) 
     IF(MOD(IINT,IREPORT)==0) THEN
            IF(Language=='zh-CN') THEN
                CALL REPORT("水流流场信息统计")
            ELSE
                CALL REPORT("FLOW FIELD STATS")
           ENDIF
      ENDIF

     CALL FLUSH(6)
!-------------UPDATE THE VISUALIZATION SERVER----------------------------------!
     CALL ARCHIVE
   
     CALL SHUTDOWN_CHECK

   END DO !!MAIN LOOP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    END MAIN LOOP OVER PHYSICAL TIME STEPS                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!==============================================================================!
!  END MAIN LOOP OVER DATA ASSIMILATION INTERVALS AND SWEEP #                  !
!==============================================================================!

!===============================================================================
!  MAIN LOOP OVER ENSEMBLE KALMAN FILTERS DATA ASSMILATION 
!===============================================================================

!==============================================================================!
!  CLOSE UP COMPUTATION                                                        !
!==============================================================================!
   IINT = IEND
   CALL CLOSEFILES
#  if defined(ENKF_ASSIM) 

#  else
   CALL ARCRST
   IF(MSR)THEN
     WRITE(IPT,*)'DUMPING RESTART'
   END IF
#  if defined (WET_DRY)
   IF(WET_DRY_ON) THEN
      write(CH8,'(I8.8)') IINT
      FNAME = 're_'//trim(CH8)//'_wd'
      CALL WD_DUMP(FNAME)
   ENDIF
#  endif
#  endif

    if(language=='zh-CN')THEN
        CALL REPORT('最终值信息')
   else
        CALL REPORT('FINAL VALUES INFORMATION')
   endif

   IF(MSR)THEN
     WRITE(IPT,*) ; WRITE(IPT,*)'Computation completed, congratulations!'
     CALL GET_CLOCK 
     CALL GETTIME(TSTRING,INT(TCURRENT-TINIT))
   END IF

   CALL PSTOP
!
!----------------------FORMAT STATEMENTS---------------------------------------|
!
7000 FORMAT(1X,A28,A13)  
7001 FORMAT(1X,A28,I8)  

   END PROGRAM DotFVM
!==============================================================================!
