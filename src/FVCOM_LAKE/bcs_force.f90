!==============================================================================|
!   Set Up the Following Boundary Conditions:                                  |
!     Bottom Freshwater (Groundwater) Info			               |
!     Freshwater River Discharge					       |
!     Meteorological Forcing						       |
!==============================================================================|

   SUBROUTINE BCS_FORCE           
!------------------------------------------------------------------------------|
   USE ALL_VARS
   USE BCS
   USE MOD_CLOCK
   USE MOD_UTILS
   USE MOD_OBCS

   IMPLICIT NONE
   CHARACTER(LEN=80) :: COMT
   REAL(SP) :: QPREC,QEVAP,WDS,WDD,HFLUX,HSHORT,RBUF
   REAL(SP) :: TX,TY,CD,TTIME,BFWTIME
   REAL(SP) :: FTEMP1,FTEMP2,FTEMP3,RBUF1,RBUF2,RBUF3
   REAL(SP), ALLOCATABLE :: RTEMP(:),RTEMP1(:,:),RTEMP2(:,:),RTEMP3(:,:)
   REAL(SP), ALLOCATABLE :: RTEMP11(:),RTEMP22(:)
   INTEGER,  ALLOCATABLE :: TEMP(:),TEMP2(:),TEMP3(:),TEMP4(:),NODE_SBC(:)
   INTEGER,  ALLOCATABLE :: TEMPD(:,:),TEMP2D(:,:),TEMP3D(:,:)
   INTEGER   I,J,K,NQTIME,NBFWTIME,ISBCN1,INMAX,IOS,NCNT,IGL,IERR,JN
   CHARACTER(LEN=13) :: TSTRING


!----------------------------REPORT--------------------------------------------!
   IF(MSR)WRITE(IPT,*  )'!'
   IF(MSR)WRITE(IPT,*)'!           SETTING UP PRESCRIBED BOUNDARY CONDITIONS   '
   IF(MSR)WRITE(IPT,*  )'!'

!==============================================================================|
!   Ground Water Information  BFWQDIS: m^3/s                                                 |
!==============================================================================|
!
!--------------determine global number of groundwater points and bcast---------!
!
   IF(MSR)THEN
     WRITE(IOPRT,*)'GROUNDWATER INFORMATION'
     READ(INBFW ,1000) COMT
     WRITE(IOPRT,1000) COMT
     READ(INBFW ,*) IBFW_GL
     WRITE(IOPRT,*) IBFW_GL
   END IF

   IBFW = 0
   IF(IBFW_GL > 0) THEN
     NCNT = 0

!
!--------------input node numbers for ground water inflow----------------------!
!
     ALLOCATE( NODE_BFW(IBFW_GL) )

     IF(MSR)THEN
     READ(INBFW ,*) (NODE_BFW(I),I=1,IBFW_GL)
     WRITE(IOPRT,*) (NODE_BFW(I),I=1,IBFW_GL)
!
!-----------------ensure all nodes exist in global domain----------------------!
!
     DO I=1,IBFW_GL
       IF(NODE_BFW(I) > MGL)THEN
         WRITE(IPT,*)'==================ERROR=================================='
         WRITE(IPT,*)'GROUND WATER NODE NUMBER',I,'IS NOT IN THE GLOBAL DOMAIN'
         WRITE(IPT,*)'ENSURE GROUNDWATER NODES <= ',MGL
         WRITE(IPT,*)'========================================================='
         CALL PSTOP
       END IF
     END DO
     END IF

!----------------------Shift To Local Domain If Parallel-----------------------!
!
     IF(SERIAL) IBFW = IBFW_GL
!
!----INPUT NUMBER OF DATA TIMES FOR GROUNDWATER DATA---------------------------!
!
     BFW_TM%LABEL = "Groundwater"
     IF(MSR)THEN
       READ(INBFW ,*) NBFWTIME
       WRITE(IOPRT,*) NBFWTIME
     END IF

     BFW_TM%NTIMES = NBFWTIME
     ALLOCATE(BFW_TM%TIMES(NBFWTIME))

!
!----READ IN FRESH WATER FLUX AT EACH TIME=BFWTIME-----------------------------!
!
     ALLOCATE(RTEMP(IBFW_GL))
     ALLOCATE(RTEMP11(IBFW_GL))  
     ALLOCATE(RTEMP22(IBFW_GL)) 

     ALLOCATE(BFWQDIS(IBFW,NBFWTIME))  ; BFWQDIS = 0.0_SP
     ALLOCATE(BFWQTDIS(IBFW,NBFWTIME)) ; BFWQTDIS = 0.0_SP 
     ALLOCATE(BFWQSDIS(IBFW,NBFWTIME)) ; BFWQSDIS = 0.0_SP 

     DO I=1,NBFWTIME
       IF(MSR)THEN
         READ(INBFW,*) BFWTIME
         WRITE(IOPRT,5000) BFWTIME
         BFW_TM%TIMES(I) = BFWTIME
         READ(INBFW,*) (RTEMP(J),J = 1,IBFW_GL)
         READ(INBFW,*) (RTEMP11(J),J = 1,IBFW_GL)
         READ(INBFW,*) (RTEMP22(J),J = 1,IBFW_GL)
       END IF


       IF(SERIAL)BFWQDIS(1:IBFW_GL,I)  = RTEMP(1:IBFW_GL)
       IF(SERIAL)BFWQTDIS(1:IBFW_GL,I) = RTEMP11(1:IBFW_GL)
       IF(SERIAL)BFWQSDIS(1:IBFW_GL,I) = RTEMP22(1:IBFW_GL)


       IF(MSR)WRITE(IOPRT,5000) (RTEMP(J),J = 1,IBFW_GL) 
       IF(MSR)WRITE(IOPRT,5000) (RTEMP11(J),J = 1,IBFW_GL) 
       IF(MSR)WRITE(IOPRT,5000) (RTEMP22(J),J = 1,IBFW_GL) 
     END DO
     DEALLOCATE(RTEMP,RTEMP11,RTEMP22)


   END IF !!IBFW_GL > 0

!
!--REPORT RESULTS--------------------------------------------------------------!
!
   ALLOCATE(TEMP(NPROCS))
   TEMP(1)  = IBFW

   IF(IBFW_GL == 0)THEN
     IF(MSR)WRITE(IPT,*)'!  GROUNDWATER FLUX      :    NONE'
   ELSE
   IF(MSR)WRITE(IPT,*)'!'
   IF(MSR)WRITE(IPT,100)'!  GROUNDWATER POINTS    :',IBFW_GL, (TEMP(I),I=1,NPROCS)
   IF(NBFWTIME > 0)THEN
     IF(MSR)WRITE(IPT,101)'!  GWATER DATA BEGIN     :',BFW_TM%TIMES(1)
     IF(MSR)WRITE(IPT,101)'!  GWATER DATA END       :',BFW_TM%TIMES(NBFWTIME)
   END IF
   END IF
   DEALLOCATE(TEMP)

!==============================================================================|
!   Input River/Dam/Intake/Outfall Boundary Values                             |
!==============================================================================|

!-------Check Selected Combination for Validity--------------------------------!
!
   REWIND(INRIV)
   READ(INRIV,'(A4,2X,A10)') INFLOW_TYPE,POINT_ST_TYPE
   IF(MSR)WRITE(IOPRT,*) 'River Inflow Information'
   IF(MSR)WRITE(IOPRT,*) 'INFLOW_TYPE==',INFLOW_TYPE
   IF(MSR)WRITE(IOPRT,*) 'POINT_ST_TYPE==',POINT_ST_TYPE

   IF(INFLOW_TYPE /= 'edge' .AND. INFLOW_TYPE /= 'node') THEN
     CALL PERROR(6,"INFLOW TYPE NOT CORRECT","SHOULD BE edge or node")
   END IF

   IF(POINT_ST_TYPE /= 'calculated' .AND. POINT_ST_TYPE /= 'specified') THEN
     CALL PERROR(6,"POINT_ST TYPE NOT CORRECT","SHOULD BE calculated or specified")
   END IF

!
!--Read in Number of Discharge Nodes/Edges-------------------------------------!
!
   IF(MSR)THEN
     READ(INRIV,*) NUMQBC_GL
   END IF

   NUMQBC = 0
   IF(NUMQBC_GL > 0)THEN
!
!--Shut off Temp/Salinity Averaging if River Flux is of type "specified"
!
!   IF(POINT_ST_TYPE == 'specified' .AND. TS_FCT)THEN
!     IF(MSR)THEN
!       WRITE(IPT,*)'=========WARNING================'
!       WRITE(IPT,*)'RIVER QUANTITIES ARE "specified"'
!       WRITE(IPT,*)'DEACTIVATING TS_FCT'
!       WRITE(IPT,*)'================================'
!     END IF
!     TS_FCT = .FALSE.
!   END IF
!
!--Read in Freshwater Discharge Nodes------------------------------------------!
!
     ALLOCATE(TEMP(NUMQBC_GL),TEMP2(NUMQBC_GL),TEMP3(NUMQBC_GL))
     IF(MSR)THEN
       DO I=1,NUMQBC_GL
         READ(INRIV,*) TEMP(I)
       END DO
     END IF
!
!--Determine Global--> Local Mapping of Freshwater Discharge Nodes
!
     IF(SERIAL)THEN
       NUMQBC = NUMQBC_GL
       IF(INFLOW_TYPE == 'node') THEN
         ALLOCATE(INODEQ(NUMQBC))
         INODEQ = TEMP
       ELSE IF(INFLOW_TYPE == 'edge') THEN
         ALLOCATE(ICELLQ(NUMQBC))
!         ICELLQ = TEMP(1:NCNT)
         ICELLQ = TEMP(1:NUMQBC)
       END IF
     END IF

     DEALLOCATE(TEMP,TEMP2,TEMP3)
!
!----Read in Freshwater Flux Vertical Distribution-----------------------------!
!

     ALLOCATE(RTEMP1(NUMQBC_GL,KBM1))
     IF(MSR)THEN
       DO I = 1, NUMQBC_GL
         READ(INRIV ,*) J,(RTEMP1(I,K),K = 1,KBM1)
         WRITE(IOPRT,*) J,(RTEMP1(I,K),K = 1,KBM1)
       END DO
     END IF

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     IF(NUMQBC > 0)THEN
         ALLOCATE(VQDIST(NUMQBC,KBM1))
         IF(SERIAL) VQDIST = RTEMP1
     END IF

     DEALLOCATE(RTEMP1)

!
!----Read in Time Dependent DataSets (DQDIS,DSDIS,DTDIS)------------------------!
!
     IF(MSR)READ(INRIV,*) NQTIME

     QBC_TM%NTIMES = NQTIME
     QBC_TM%LABEL  = "Freshwater Discharge"
     ALLOCATE(QBC_TM%TIMES(NQTIME))
     ALLOCATE(RTEMP1(NUMQBC_GL,NQTIME))
     ALLOCATE(RTEMP2(NUMQBC_GL,NQTIME))
     ALLOCATE(RTEMP3(NUMQBC_GL,NQTIME))

     IF(MSR)THEN
       DO I = 1, NQTIME
         READ(INRIV,*) TTIME
         QBC_TM%TIMES(I) = TTIME
         READ(INRIV,*) (RTEMP1(J,I),J = 1,NUMQBC_GL)
         READ(INRIV,*) (RTEMP2(J,I),J = 1,NUMQBC_GL)
         READ(INRIV,*) (RTEMP3(J,I),J = 1,NUMQBC_GL)
         WRITE(IOPRT,5000) TTIME
         WRITE(IOPRT,5000) (RTEMP1(J,I),J = 1,NUMQBC_GL)
         WRITE(IOPRT,5000) (RTEMP2(J,I),J = 1,NUMQBC_GL)
         WRITE(IOPRT,5000) (RTEMP3(J,I),J = 1,NUMQBC_GL)
       END DO
     END IF

!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     IF(NUMQBC > 0)THEN
       ALLOCATE(DQDIS(NUMQBC,NQTIME))
       ALLOCATE(DTDIS(NUMQBC,NQTIME))
       ALLOCATE(DSDIS(NUMQBC,NQTIME))

       IF(SERIAL)THEN
         DQDIS(1:NUMQBC_GL,:) = RTEMP1(1:NUMQBC_GL,:)
         DTDIS(1:NUMQBC_GL,:) = RTEMP2(1:NUMQBC_GL,:)
         DSDIS(1:NUMQBC_GL,:) = RTEMP3(1:NUMQBC_GL,:)
       END IF
     END IF

     DEALLOCATE(RTEMP1,RTEMP2,RTEMP3)

   CLOSE(INRIV)
!
!--REPORT RESULTS--------------------------------------------------------------!
!
   ALLOCATE(TEMP(NPROCS))
   TEMP(1)  = NUMQBC
   FTEMP1 = 0.0_SP; FTEMP2 = 0.0_SP; FTEMP3 = 0.0_SP;
   IF(NUMQBC > 0) FTEMP1 = MAXVAL(DQDIS)
   IF(NUMQBC > 0) FTEMP2 = MAXVAL(DTDIS)
   IF(NUMQBC > 0) FTEMP3 = MAXVAL(DSDIS)
   RBUF1 = FTEMP1 ; RBUF2 = FTEMP2 ; RBUF3 = FTEMP3

   END IF !! NUMQBC_GL > 0

   IF(MSR)WRITE(IPT,*)'!'
   IF(NUMQBC_GL == 0)THEN
     IF(MSR)WRITE(IPT,*)'!  FRESHWATER FLUX       :    NONE'
   ELSE
     IF(MSR)WRITE(IPT,100)'!  FRESHWATER POINTS     :',NUMQBC_GL, (TEMP(I),I=1,NPROCS)
     IF(MSR)CALL GETTIME(TSTRING,3600*INT(QBC_TM%TIMES(1)))
     IF(MSR)WRITE(IPT,102)'!  FWATER DATA BEGIN     :  ',TSTRING
     IF(MSR)CALL GETTIME(TSTRING,3600*INT(QBC_TM%TIMES(QBC_TM%NTIMES)))
     IF(MSR)WRITE(IPT,102)'!  FWATER DATA END       :  ',TSTRING
     IF(MSR)WRITE(IPT,101)'!  MAX DQDIS             :',RBUF1
     IF(MSR)WRITE(IPT,101)'!  MAX DTDIS             :',RBUF2
     IF(MSR)WRITE(IPT,101)'!  MAX DSDIS             :',RBUF3
     DEALLOCATE(TEMP)
   END IF


!==============================================================================|
!   Input Meteorological Boundary Conditions                                   |
!==============================================================================|
!    precipitation: mm/s       "qprec"                                         |
!    evaporation:   mm/s       "qevap"                                         |
!    wind:          wds (speed) wdd (direction)       			       |
!    heat flux:     w/m^2                              			       |
!==============================================================================|

   IF(M_TYPE == 'uniform')THEN

!==============================================================================|
!   UNIFORM METEOLOGICAL CONDITIONS                                            |
!==============================================================================|

     READ(INCWH,1000) COMT
     IF(MSR)WRITE(IOPRT,*)'Meteorological Forcing Info'
     IF(MSR)WRITE(IOPRT,1000) COMT

!
!----Determine Number of Data Times--------------------------------------------!
!
     NCNT = 0
     DO WHILE(.TRUE.)
       READ(INCWH,*,END=15,IOSTAT=IOS)
       READ(INCWH,*,END=15,IOSTAT=IOS)
       IF(IOS < 0)EXIT
       NCNT = NCNT + 1
     END DO
 15  CONTINUE
     IF(NCNT == 0)CALL PERROR(6,"NO UNIFORM METEO DATA PROVIDED")

     REWIND(INCWH) ; READ(INCWH,*)

!
!----Read in Precipitation/Evap/Wind/Heat Flux/Radiation Data at Each Time-----!
!

     UMF_TM%NTIMES = NCNT
     ALLOCATE(UMF_TM%TIMES(NCNT))
     ALLOCATE(UQPREC(NCNT),UQEVAP(NCNT))
     ALLOCATE(UWIND(NCNT),VWIND(NCNT))
     ALLOCATE(UHFLUX(NCNT),UHSHORT(NCNT))

     DO I=1,NCNT
       READ(INCWH ,*) TTIME
       IF(MSR)WRITE(IOPRT,*) TTIME
       UMF_TM%TIMES(I) = TTIME

       READ(INCWH ,*) QPREC, QEVAP, WDS, WDD, HFLUX,HSHORT

       IF(MSR)WRITE(IOPRT,5000) QPREC, QEVAP, WDS, WDD, HFLUX,HSHORT

!       UQPREC(I) = QPREC / (86400.0_SP*365.0_SP)
!       UQEVAP(I) = QEVAP / (86400.0_SP*365.0_SP)
!       UQPREC(I) = QPREC / 1000.0_SP
!       UQEVAP(I) = QEVAP / 1000.0_SP
       UQPREC(I) = QPREC
       UQEVAP(I) = QEVAP

       WDD = MOD(WDD,360.0_SP)
!       UWIND(I) = WDS * COS(6.28319_SP*WDD/360.0_SP)
!       VWIND(I) = WDS * SIN(6.28319_SP*WDD/360.0_SP)
       UWIND(I) = -WDS * SIN(WDD*DEG2RAD)
       VWIND(I) = -WDS * COS(WDD*DEG2RAD)

       UHFLUX(I)  = HFLUX
       UHSHORT(I) = HSHORT
     END DO

     IF(WINDTYPE /= 'speed' .AND. WINDTYPE /='stress') THEN
       WRITE(IPT,*)'==================ERROR=================================='
       WRITE(IPT,*)'NO UNIFORM METEO DATA PROVIDED'
       WRITE(IPT,*)'WINDTYPE IS NOT CORRECT, --->',WINDTYPE
       WRITE(IPT,*)'MUST BE "speed" or "stress"'
       WRITE(IPT,*)'========================================================='
       CALL PSTOP
     END IF

     CLOSE(INCWH)

!
!--REPORT RESULTS--------------------------------------------------------------!
!

   IF(MSR)THEN
     WRITE(IPT,*)'!'
     WRITE(IPT,*    )'!  UNIFORM METEO         :    SET'
      IF(UMF_TM%NTIMES > 0)THEN
        CALL GETTIME(TSTRING,INT(3600.*UMF_TM%TIMES(1)))
        WRITE(IPT,102)'!  METEO DATA BEGIN      :  ',TSTRING
        CALL GETTIME(TSTRING,INT(3600.*UMF_TM%TIMES(UMF_TM%NTIMES)))
        WRITE(IPT,102)'!  METEO DATA END        :  ',TSTRING
      END IF
    END IF

!==============================================================================|
!   NON-UNIFORM METEOLOGICAL CONDITIONS                                        |
!==============================================================================|

   ELSE IF (M_TYPE == 'non-uniform')THEN

!=====================HEAT FLUX/SHORT WAVE RADIATION===========================!
     REWIND(INHFX)
!
!----Input Number of Data Times for Heat Flux and Short Wave Radiation---------!
!
     IF(MSR)THEN
     NCNT = 0
     DO WHILE(.TRUE.)
       READ(INHFX,END=10)FTEMP1
       READ(INHFX)
       NCNT = NCNT + 1
     END DO
 10  CONTINUE
     REWIND(INHFX)

     IF(NCNT == 0)CALL PERROR(6,"NO DATA PROVIDED FOR HEAT FLUX AND SHORT WAVE RAD")
     END IF

     HFX_TM%NTIMES = NCNT 

!
!----Read in Data Times and Global Heat Flux/Short Wave Radiation Data---------!
!

     ALLOCATE(HFX_TM%TIMES(HFX_TM%NTIMES))
     ALLOCATE(RTEMP1(MGL,HFX_TM%NTIMES),RTEMP2(MGL,HFX_TM%NTIMES))

     IF(MSR)THEN
       DO J=1,HFX_TM%NTIMES
         READ(INHFX) HFX_TM%TIMES(J) 
         READ(INHFX) (RTEMP1(I,J),RTEMP2(I,J),I=1,MGL)
       END DO
     END IF
!
!----Broadcast Data------------------------------------------------------------!
!
!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     ALLOCATE(DHFLUX(M,HFX_TM%NTIMES),DHSHORT(M,HFX_TM%NTIMES))

     IF(SERIAL)THEN
       DHFLUX(1:MGL,:)  = RTEMP1(1:MGL,:)
       DHSHORT(1:MGL,:) = RTEMP2(1:MGL,:)
     END IF
     DEALLOCATE(RTEMP1,RTEMP2)
     IF(MSR)WRITE(IPT,101)'!  HFLUX/SWRAD READ      :    COMPLETE'
      

!=====================TIME DEPENDENT WIND FIELD================================!


     REWIND(INWND)
!
!----Input Number of Data Times for Wind Field---------------------------------!
!
     IF(MSR)THEN
     NCNT = 0
     DO WHILE(.TRUE.)
       READ(INWND,END=20)FTEMP1
       READ(INWND) 
       NCNT = NCNT + 1       
     END DO
 20  CONTINUE
     REWIND(INWND)

     IF(NCNT == 0)CALL PERROR(6,"NO DATA PROVIDED FOR SURFACE WIND FIELD")
     END IF

     WND_TM%NTIMES = NCNT 


!
!----Read in Data Times and Global Wind Data-----------------------------------!
!

     ALLOCATE(WND_TM%TIMES(WND_TM%NTIMES))
     ALLOCATE(RTEMP1(NGL,WND_TM%NTIMES),RTEMP2(NGL,WND_TM%NTIMES))

     IF(MSR)THEN
       DO J=1,WND_TM%NTIMES
         READ(INWND) WND_TM%TIMES(J) 
         READ(INWND) (RTEMP1(I,J),RTEMP2(I,J),I=1,NGL)
       END DO
     END IF

!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
     ALLOCATE(DTX(N,WND_TM%NTIMES),DTY(N,WND_TM%NTIMES))

     IF(SERIAL)THEN
       DTX(1:NGL,:)  = RTEMP1(1:NGL,:)
       DTY(1:NGL,:) = RTEMP2(1:NGL,:)
     END IF

     DEALLOCATE(RTEMP1,RTEMP2)
     IF(MSR)WRITE(IPT,101)'!  WIND FIELD READ       :    COMPLETE'
      


!=====================TIME DEPENDENT EVAPORATION AND PRECIPITATION=============!

     IF(EVP_FLAG)THEN
       REWIND(INEVP)
!
!----Input Number of Data Times for Evaporation and Precipitation--------------!
!
       IF(MSR)THEN
         NCNT = 0
         DO WHILE(.TRUE.)
           READ(INEVP,END=30)FTEMP1
           READ(INEVP) 
           NCNT = NCNT + 1       
         END DO
 30      CONTINUE
         REWIND(INEVP)

         IF(NCNT == 0)CALL PERROR(6,"NO DATA PROVIDED FOR EVAPORATION AND PRECIPITATION")
       END IF

       EVP_TM%NTIMES = NCNT 

!
!----Read in Data Times and Global Evaporation and Precipitation---------------!
!

       ALLOCATE(EVP_TM%TIMES(EVP_TM%NTIMES))
       ALLOCATE(RTEMP1(MGL,EVP_TM%NTIMES),RTEMP2(MGL,EVP_TM%NTIMES))

       IF(MSR)THEN
         DO J=1,EVP_TM%NTIMES
           READ(INEVP) EVP_TM%TIMES(J) 
           READ(INEVP) (RTEMP1(I,J),RTEMP2(I,J),I=1,MGL)
         END DO
       END IF
!
!----TRANSFORM TO LOCAL ARRAYS-------------------------------------------------|
!
       ALLOCATE(DQEVAP(M,EVP_TM%NTIMES),DQPREC(M,EVP_TM%NTIMES))

       IF(SERIAL)THEN
         DQEVAP(1:MGL,:) = RTEMP1(1:MGL,:)
         DQPREC(1:MGL,:) = RTEMP2(1:MGL,:)
       END IF

       DEALLOCATE(RTEMP1,RTEMP2)
       IF(MSR)WRITE(IPT,101)'!  EVAPORATION AND PRECIPITATION READ : COMPLETE'
     END IF 

!
!--REPORT RESULTS--------------------------------------------------------------!
!

     IF(MSR)WRITE(IPT,*)'!'
     IF(MSR)WRITE(IPT,*    )'!  NON-UNIFORM METEO     :    SET'
     IF(HFX_TM%NTIMES > 0)THEN
       CALL GETTIME(TSTRING,3600*INT(HFX_TM%TIMES(1)))
       IF(MSR)WRITE(IPT,102)'!  HEAT/RAD DATA BEGIN   :  ',TSTRING        
       CALL GETTIME(TSTRING,3600*INT(HFX_TM%TIMES(HFX_TM%NTIMES)))
       IF(MSR)WRITE(IPT,102)'!  HEAT/RAD DATA END     :  ',TSTRING
     END IF
     IF(WND_TM%NTIMES > 0)THEN
       CALL GETTIME(TSTRING,3600*INT(WND_TM%TIMES(1)))
       IF(MSR)WRITE(IPT,102)'!  WIND DATA BEGIN       :  ',TSTRING          
       CALL GETTIME(TSTRING,3600*INT(WND_TM%TIMES(WND_TM%NTIMES)))
       IF(MSR)WRITE(IPT,102)'!  WIND DATA END         :  ',TSTRING
     END IF
     IF(EVP_FLAG)THEN
       IF(EVP_TM%NTIMES > 0)THEN
         CALL GETTIME(TSTRING,3600*INT(EVP_TM%TIMES(1)))
         IF(MSR)WRITE(IPT,102)'!  EVAP/PREC DATA BEGIN       :  ',TSTRING          
         CALL GETTIME(TSTRING,3600*INT(EVP_TM%TIMES(EVP_TM%NTIMES)))
         IF(MSR)WRITE(IPT,102)'!  EVAP/PREC DATA END         :  ',TSTRING
       END IF
     END IF
     
     FTEMP1 = SUM(DHFLUX/FLOAT(M))/FLOAT(HFX_TM%NTIMES)
     FTEMP2 = SUM(DHSHORT)/FLOAT(M*HFX_TM%NTIMES)
     IF(SERIAL)THEN
       RBUF1 = FTEMP1
       RBUF2 = FTEMP2
     END IF

     IF(MSR)WRITE(IPT,101)'!  AVE HEAT FLUX         :',RBUF1/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,101)'!  AVE SHORT WAVE RAD    :',RBUF2/FLOAT(NPROCS)

!     FTEMP1 = SUM(DTX)/FLOAT(M*WND_TM%NTIMES)
!     FTEMP2 = SUM(DTY)/FLOAT(M*WND_TM%NTIMES)
     FTEMP1 = SUM(DTX)/FLOAT(N*WND_TM%NTIMES)
     FTEMP2 = SUM(DTY)/FLOAT(N*WND_TM%NTIMES)
     IF(SERIAL)THEN
       RBUF1 = FTEMP1
       RBUF2 = FTEMP2
     END IF

     IF(MSR)WRITE(IPT,101)'!  AVE WIND X-COMP       :',RBUF1/FLOAT(NPROCS)
     IF(MSR)WRITE(IPT,101)'!  AVE WIND Y-COMP       :',RBUF2/FLOAT(NPROCS)

     IF(EVP_FLAG)THEN
       FTEMP1 = SUM(DQEVAP)/FLOAT(M*EVP_TM%NTIMES)
       FTEMP2 = SUM(DQPREC)/FLOAT(M*EVP_TM%NTIMES)
       IF(SERIAL)THEN
         RBUF1 = FTEMP1
         RBUF2 = FTEMP2
       END IF

       IF(MSR)WRITE(IPT,101)'!  AVE EVAPORATION       :',RBUF1/FLOAT(NPROCS)
       IF(MSR)WRITE(IPT,101)'!  AVE PRECIPITATION     :',RBUF2/FLOAT(NPROCS)
     END IF
   ELSE
     WRITE(IPT,*)'==================ERROR=================================='
     WRITE(IPT,*)'M_TYPE NOT CORRECT, --->',M_TYPE
     WRITE(IPT,*)'MUST BE "uniform" or "non-uniform"'
     WRITE(IPT,*)'========================================================='
     CALL PSTOP
   END IF
!
!--Format Statements-----------------------------------------------------------!
!

   100  FORMAT(1X,A26,I6," =>",2X,4(I5,1H,))
   101  FORMAT(1X,A26,F10.4)  
   102  FORMAT(1X,A28,A13)  
   1000 FORMAT(A80)
   5000 FORMAT(8E14.5)

   RETURN
   END SUBROUTINE BCS_FORCE
!==============================================================================|
