
   SUBROUTINE out_binary_cfd(JSEXPLORER,IINTT)

!------------------------------------------------------------------------------|

   USE ALL_VARS
   USE MOD_WQM

   IMPLICIT NONE
   INTEGER:: JSEXPLORER
   INTEGER :: L1,K,DATATYPE
   INTEGER :: VAR_COUNT,AVR_COUNT,STEP
   INTEGER, INTENT(IN) :: IINTT
   CHARACTER(LEN=100) ::FILENAME_HD,FILENAME_TEMP,FILENAME_WQ,FILENAME_SALT
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: WQMEAN
   REAL(SP) ::uvaverage,taverage,angle, salt_avrg

!==============================================================================|

   FILENAME_HD = TRIM(OUTDIR)//"/out/"//trim(casename)//' HydroDynamics.csd'
   FILENAME_TEMP =  TRIM(OUTDIR)//"/out/"//trim(casename)//' Temperature.csd'
   FILENAME_WQ =  TRIM(OUTDIR)//"/out/"//trim(casename)//' WaterQuality.csd'
   FILENAME_SALT =  TRIM(OUTDIR)//"/out/"//trim(casename)//' Transport.csd'
!------------------------------------------------------------------------------!
!  OPEN FILE (Name Based on Iteration Number)                                  !
!------------------------------------------------------------------------------!
!C1-------------------------------------------------------------------------------------
      ! **  INITIAL CALL
      IF(JSEXPLORER == 1)THEN  
           STEP = IEND / IDMPSMS
       ! *** WATER QUALITY MODEL (HEM3D) RESULTS
          OPEN(95,FILE=FILENAME_HD,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='binary')
          CLOSE(95,STATUS='DELETE')

           OPEN(96,FILE=FILENAME_TEMP,STATUS='UNKNOWN', ACCESS='SEQUENTIAL',FORM='binary')
          CLOSE(96,STATUS='DELETE')
           
           DATATYPE=0
          OPEN(95,FILE=FILENAME_HD,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='binary')
          WRITE(95)DATATYPE,N,M,KBM1,STEP
          CLOSE(95,STATUS='KEEP')

          DATATYPE=1
          OPEN(96,FILE=FILENAME_TEMP,STATUS='UNKNOWN',  ACCESS='SEQUENTIAL',FORM='binary')
          WRITE(96)DATATYPE,N,M,KBM1,STEP
          CLOSE(96,STATUS='KEEP')

           IF(WQM_ON)THEN        
                    OPEN(97,FILE=FILENAME_WQ,STATUS='UNKNOWN', ACCESS='SEQUENTIAL',FORM='binary')
                    CLOSE(97,STATUS='DELETE')

                    DATATYPE= 2
                    OPEN(97,FILE=FILENAME_WQ,STATUS='UNKNOWN',  ACCESS='SEQUENTIAL',FORM='binary')
                    WRITE(97)DATATYPE,N,M,KBM1,STEP
                    CLOSE(97,STATUS='KEEP')
           ENDIF

          IF(SALINITY_ON)THEN
                  OPEN(98,FILE=FILENAME_SALT,STATUS='UNKNOWN', ACCESS='SEQUENTIAL',FORM='binary')
                  CLOSE(98,STATUS='DELETE')

                  DATATYPE=3
                  OPEN(98,FILE=FILENAME_SALT,STATUS='UNKNOWN',  ACCESS='SEQUENTIAL',FORM='binary')
                  WRITE(98)DATATYPE,N,M,KBM1,STEP
                  CLOSE(98,STATUS='KEEP')
          ENDIF
      ENDIF
!C1-------------------------------------------------------------------------------------
!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Single Processor Case)                                !
!------------------------------------------------------------------------------!

       IF(SERIAL)THEN

             OPEN(95,FILE=FILENAME_HD,STATUS='UNKNOWN', POSITION='APPEND',FORM='binary')     
             VAR_COUNT = 3
             AVR_COUNT =  2
             WRITE(95)IINTT,VAR_COUNT,AVR_COUNT

             DO L1=1,N
                call VectorAzimuth(UA(L1),VA(L1),uvaverage,angle)
                WRITE(95) uvaverage,angle
                DO K=1,KBM1
                       call VectorAzimuth(U(L1,K),V(L1,K),uvaverage,angle)
                       WRITE(95)uvaverage,angle,WW(L1,K)
                END DO
             END DO
             CALL FLUSH(95)
             CLOSE(95,STATUS='KEEP')

            OPEN(96,FILE=FILENAME_TEMP,STATUS='UNKNOWN', POSITION='APPEND',FORM='binary') 
             VAR_COUNT =  1
             AVR_COUNT =  2
            WRITE(96)IINTT,VAR_COUNT,AVR_COUNT

              DO L1=1,M
                 taverage = SUM(T1(L1,:))/KBM1
                 WRITE(96)BaseSurfaceElevation+EL(L1),taverage,(T1(L1,K),K=1,KBM1)
             END DO
            CALL FLUSH(96)
            CLOSE(96,STATUS='KEEP')

            IF(WQM_ON)THEN
                  ALLOCATE(WQMEAN(NB))
                  OPEN(97,FILE=FILENAME_WQ,STATUS='UNKNOWN', POSITION='APPEND',FORM='binary') 
                 VAR_COUNT =  0
                 AVR_COUNT =  8
                WRITE(97)IINTT,VAR_COUNT,AVR_COUNT

                  DO L1=1,M
                        do K = 1,NB
                            WQMEAN(K)= SUM(WQM(L1,1:KBM2,K)) / KBM2
                        end do
                        WRITE(97) (WQMEAN(K),K=1,NB)
                 END DO
                CALL FLUSH(97)
                CLOSE(97,STATUS='KEEP')
            ENDIF

             IF(SALINITY_ON)THEN
                       OPEN(98,FILE=FILENAME_SALT,STATUS='UNKNOWN', POSITION='APPEND',FORM='binary') 
                         VAR_COUNT =  1
                         AVR_COUNT =  1
                        WRITE(98)IINTT,VAR_COUNT,AVR_COUNT

                        DO L1=1,M
                                salt_avrg = SUM(S1(L1,:))/KBM1
                                WRITE(98) salt_avrg,(S1(L1,K),K=1,KBM1)
                        END DO
                        CALL FLUSH(98)
                        CLOSE(98,STATUS='KEEP')
             ENDIF

       END IF

   RETURN
   END SUBROUTINE out_binary_cfd   
!==============================================================================|
