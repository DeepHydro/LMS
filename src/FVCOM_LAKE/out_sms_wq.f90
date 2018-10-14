!==============================================================================|
!   Write Output Files Used by SMS Post-Processing Software                    |
!==============================================================================|

   SUBROUTINE OUT_SMS_WQ(IINTT)            

!------------------------------------------------------------------------------|
   USE MOD_WQM
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IINTT
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: WQMEAN
    REAL(SP), ALLOCATABLE, DIMENSION(:)::TN,TP
   INTEGER :: I1,I,IOSMSWQ,K, MINDEX
   CHARACTER(LEN=4) :: FILENUMBER

!==============================================================================|
   
!------------------------------------------------------------------------------!
!  OPEN AND REWIND FILES                                                       !
!------------------------------------------------------------------------------!
  IOSMSWQ = 100
   IF(MSR)THEN
      WRITE(FILENUMBER,'(I4.4)') IINTT
     OPEN(IOSMSWQ,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//FILENUMBER//'_wq.xy',STATUS='unknown')
     REWIND(IOSMSWQ)
   END IF

  ALLOCATE(WQMEAN(M,NB))
  ALLOCATE(TN(M))
  ALLOCATE(TP(M))
!------------------------------------------------------------------------------!
!  WRITE TO FILES (SERIAL EXECUTION)                                           !
!------------------------------------------------------------------------------!

!  WRITE NODAL SURFACE GRID COORDINATES, ELEVATION, SALINITY, AND TEMPERATURE-
   IF(SERIAL)THEN
     WRITE(IOSMSWQ,10)
     WRITE(IOSMSWQ,30) M,IINTT,10
     DO MINDEX=1,M
        do I=1,NB
!           WQMEAN(MINDEX,I)= SUM(WQM(MINDEX,:,I))/KBM1
            WQMEAN(MINDEX,I)= SUM(WQM(MINDEX,1:KBM2,I))/KBM2
        end do
          TN(MINDEX)= WQM(MINDEX,1,4)+WQM(MINDEX,1,5)+WQM(MINDEX,1,6)
          TP(MINDEX)= WQM(MINDEX,1,7)+WQM(MINDEX,1,8)
!       WRITE(IOSMSWQ,'(10E17.8)')VX(I1)+VXMIN,VY(I1)+VYMIN,WQMEAN(1),WQMEAN(2),WQMEAN(3),WQMEAN(4)&
!       ,WQMEAN(5),WQMEAN(6),WQMEAN(7),WQMEAN(8
        !WRITE(IOSMSWQ,'(10E17.8)')VX(MINDEX)+VXMIN,VY(MINDEX)+VYMIN,(WQMEAN(MINDEX,K),K=1,NB)
        WRITE(IOSMSWQ,'(12E17.8)')VX(MINDEX)+VXMIN,VY(MINDEX)+VYMIN,(WQM(MINDEX,1,K),K=1,NB),TN(MINDEX),TP(MINDEX)
     END DO
     WRITE(IOSMSWQ,*) 'IINT==',IINT
     WRITE(IOSMSWQ,*) 'T1EL==',TIME*86400.0_SP
     WRITE(IOSMSWQ,*) 'T2EL==',T2EL+DELTT
  ENDIF   

  DEALLOCATE(WQMEAN)
  DEALLOCATE(TN)
  DEALLOCATE(TP)
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
30 FORMAT('xyd ',I10,' wq',I4.4,' ', I2,' DO CBOD PHYT NH3 NO3 ON OPO4 OP TN TP')
   RETURN
   END SUBROUTINE OUT_SMS_WQ
!==============================================================================|