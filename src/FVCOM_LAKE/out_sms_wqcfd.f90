!==============================================================================|
!   Write Output Files Used by SMS Post-Processing Software                    |
!==============================================================================|

   SUBROUTINE out_sms_wqcfd(JSEXPLORER)            

!------------------------------------------------------------------------------|
   USE MOD_WQM
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IINTT
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: WQMEAN
   INTEGER :: I1,I,IOSMSWQ,K,L1
   CHARACTER(LEN=4) :: FILENUMBER
   INTEGER:: JSEXPLORER
   CHARACTER(LEN=100) ::FILENAME_WQ
   INTEGER  :: J1,J2,END_AVGE,DATATYPE
!==============================================================================|
   
   FILENAME_WQ = TRIM(OUTDIR)//"/sms/"//trim(casename)//'WaterQuality.csd'
   END_AVGE = BEG_AVGE + NUM_AVGE*INT_AVGE - 1

    IF(JSEXPLORER == 1)THEN  
       ! *** WATER QUALITY MODEL (HEM3D) RESULTS
          OPEN(95,FILE=FILENAME_WQ,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='binary')
          CLOSE(95,STATUS='DELETE')
           
           DATATYPE=0
          OPEN(95,FILE=FILENAME_WQ,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='binary')
          WRITE(95)DATATYPE,N,M,KBM1,IEND
          CLOSE(95,STATUS='KEEP')
      ENDIF

   J1 = MOD((IINT+1-BEG_AVGE),INT_AVGE)
   J2 = (IINT+1-BEG_AVGE)/INT_AVGE

   IF(IINT >= BEG_AVGE  .AND. IINT <= END_AVGE)THEN
     OPEN(95,FILE=FILENAME_WQ,STATUS='UNKNOWN', POSITION='APPEND',FORM='binary')     
     WRITE(95)IINT
      do I=1,NB
           WQMEAN(I)= SUM(WQM(I1,:,I))/KBM1
        end do

     DO L1=1,N
         uvaverage= SQRT (UA(L1)*UA(L1)+VA(L1)*VA(L1))
         WRITE(95)uvaverage
        DO K=1,KBM1
               WRITE(95)U(L1,K),V(L1,K),WW(L1,K)
        END DO
     END DO
     CALL FLUSH(95)
     CLOSE(95,STATUS='KEEP')

   ENDIF

  ALLOCATE(WQMEAN(NB))
!------------------------------------------------------------------------------!
!  WRITE TO FILES (SERIAL EXECUTION)                                           !
!------------------------------------------------------------------------------!

!  WRITE NODAL SURFACE GRID COORDINATES, ELEVATION, SALINITY, AND TEMPERATURE-
   IF(SERIAL)THEN
     WRITE(IOSMSWQ,10)
     WRITE(IOSMSWQ,30) M,IINTT,8
     DO I1=1,M
        do I=1,NB
           WQMEAN(I)= SUM(WQM(I1,:,I))/KBM1
        end do
!       WRITE(IOSMSWQ,'(10E17.8)')VX(I1)+VXMIN,VY(I1)+VYMIN,WQMEAN(1),WQMEAN(2),WQMEAN(3),WQMEAN(4)&
!       ,WQMEAN(5),WQMEAN(6),WQMEAN(7),WQMEAN(8
    WRITE(IOSMSWQ,'(10E17.8)')VX(I1)+VXMIN,VY(I1)+VYMIN,(WQMEAN(K),K=1,NB)
     END DO
     WRITE(IOSMSWQ,*) 'IINT==',IINT
     WRITE(IOSMSWQ,*) 'T1EL==',TIME*86400.0_SP
     WRITE(IOSMSWQ,*) 'T2EL==',T2EL+DELTT
  ENDIF   
!------------------------------------------------------------------------------!
!  CLOSE UP FILES                                                              !
!------------------------------------------------------------------------------!

   IF(MSR)CLOSE(IOSMSWQ)
!
!--FORMATS---------------------------------------------------------------------!
!

10 FORMAT('scat2d')
!30 FORMAT('xyd ',I10,' wq',I2,'DO CBOD PHYT NH3 NO3 ON OPO4 OP')
!30 FORMAT('xyd ',I10,' wq',I2,'do cbod')
30 FORMAT('xyd ',I10,' wq',I4.4,' ', I2,' DO CBOD PHYT NH3 NO3 ON OPO4 OP')
   RETURN
   END SUBROUTINE out_sms_wqcfd
!==============================================================================|