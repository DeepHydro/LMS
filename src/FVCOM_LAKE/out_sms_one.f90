!==============================================================================|
!   Write Output Files Used by SMS Post-Processing Software                    |
!==============================================================================|

   SUBROUTINE OUT_SMS_ONE(IINTT)            

!------------------------------------------------------------------------------|
   USE ALL_VARS

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IINTT
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,T1TMP,S1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: UATMP,VATMP,UUWINDTMP,VVWINDTMP,ELTMP 
   INTEGER :: I1,I,K
   CHARACTER(LEN=4) :: FILENUMBER
   CHARACTER(LEN=60) ::ELT_FORMAT
   CHARACTER(LEN=60) ::UV_FORMAT
     CHARACTER(LEN=4) ::UNAME
     CHARACTER(LEN=4) ::VNAME
     CHARACTER(LEN=4) ::TNAME
    REAL(SP) ::uv_avrg
    REAL(SP) ::temp_avrg,salt_avrg
!==============================================================================|
   
!------------------------------------------------------------------------------!
!  OPEN AND REWIND FILES                                                       !
!------------------------------------------------------------------------------!

   IF(MSR)THEN
      WRITE(FILENUMBER,'(I4.4)') IINTT
!     OPEN(IOSMSV,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//'_uvi_uva.xy',STATUS='unknown')
!     OPEN(IOSMST,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//'_elts.xy',STATUS='unknown')
     OPEN(IOSMSV,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//FILENUMBER//'_uv.xy',STATUS='unknown')
     OPEN(IOSMST,FILE=TRIM(OUTDIR)//"/sms/"//trim(casename)//FILENUMBER//'_elts.xy',STATUS='unknown')
     REWIND(IOSMSV)
     REWIND(IOSMST)
   END IF

!------------------------------------------------------------------------------!
!  WRITE TO FILES (SERIAL EXECUTION)                                           !
!------------------------------------------------------------------------------!


   IF(SMS_OUT_Mode == 'Simple') THEN
       IF(SERIAL)THEN
        !  WRITE NODAL SURFACE GRID COORDINATES, ELEVATION, SALINITY, AND TEMPERATURE-
         WRITE(IOSMST,10)
         WRITE(IOSMST,30) M,4
         DO I1=1,M
           temp_avrg = SUM(T1(I1,:))/KB
           WRITE(IOSMST,'(6E17.8)')VX(I1)+VXMIN,VY(I1)+VYMIN,BaseSurfaceElevation+EL(I1),T1(I1,1),S1(I1,1),temp_avrg
         END DO
         WRITE(IOSMST,*) 'IINT==',IINT
         WRITE(IOSMST,*) 'T1EL==',TIME*86400.0_SP
         WRITE(IOSMST,*) 'T2EL==',T2EL+DELTT

        !  WRITE ELEMENT SURFACE GRID, VELOCITY, VERT-AVGED VELOCITY, AND WIND DATA
           WRITE(IOSMSV,10)
           WRITE(IOSMSV,20) N,7
           DO I1=1,N
        !     WRITE(IOSMSV,'(8E17.8)') &
        !       XC(I1)+VXMIN,YC(I1)+VYMIN,U(I1,1),V(I1,1),UA(I1),VA(I1),UUWIND(I1),VVWIND(I1)
                uv_avrg= SQRT (UA(I1)*UA(I1)+VA(I1)*VA(I1))
                WRITE(IOSMSV,'(9E17.8)') &
                XC(I1)+VXMIN,YC(I1)+VYMIN,U(I1,1),V(I1,1),UA(I1),VA(I1),U(I1,KB),V(I1,KB),uv_avrg
           END DO
           WRITE(IOSMSV,*) 'IINT==',IINT
      END IF

    ELSE
        !  WRITE NODAL SURFACE GRID COORDINATES, ELEVATION, SALINITY, AND TEMPERATURE-
        WRITE(IOSMST,10)
        WRITE(ELT_FORMAT,50) M,IINTT,KBM1+3
          DO K=1,KBM1
             WRITE(TNAME,51) K
             WRITE(ELT_FORMAT,*)  TRIM(ELT_FORMAT)//TRIM(TNAME)
         ENDDO
        WRITE(IOSMST,*)ELT_FORMAT

         DO I1=1,M
           temp_avrg = SUM(T1(I1,:))/KBM1
           salt_avrg = SUM(S1(I1,:))/KBM1
           WRITE(IOSMST,'(11E17.8)')VX(I1)+VXMIN,VY(I1)+VYMIN,BaseSurfaceElevation+EL(I1),temp_avrg, salt_avrg, (T1(I1,K),K=1,KBM1)
         END DO
         WRITE(IOSMST,*) 'IINT==',IINT
         WRITE(IOSMST,*) 'T1EL==',TIME*86400.0_SP
         WRITE(IOSMST,*) 'T2EL==',T2EL+DELTT

          !  WRITE ELEMENT SURFACE GRID, VELOCITY, VERT-AVGED VELOCITY, AND WIND DATA
         WRITE(IOSMSV,10)
         WRITE(UV_FORMAT,40) N,IINTT,KBM1*2+3
         DO K=1,KBM1
             WRITE(UNAME,41) K
             WRITE(VNAME,42) K
             WRITE(UV_FORMAT,*)  TRIM(UV_FORMAT)//TRIM(UNAME)//TRIM(VNAME)
         ENDDO
        WRITE(IOSMSV,*)UV_FORMAT

           DO I1=1,N
                uv_avrg= SQRT (UA(I1)*UA(I1)+VA(I1)*VA(I1))
                WRITE(IOSMSV,'(15E17.8)') &
                XC(I1)+VXMIN,YC(I1)+VYMIN,UA(I1),VA(I1),uv_avrg,(U(I1,K),V(I1,K),K=1,KBM1)
           END DO
           WRITE(IOSMSV,*) 'IINT==',IINT
ENDIF
!------------------------------------------------------------------------------!
!  CLOSE UP FILES                                                              !
!------------------------------------------------------------------------------!

   IF(MSR)CLOSE(IOSMSV)
   IF(MSR)CLOSE(IOSMST)

!
!--FORMATS---------------------------------------------------------------------!
!

10 FORMAT('scat2d')
20 FORMAT('xyd ',I10,' suva ',I2,' su sv ua va bu bv uva')
30 FORMAT('xyd ',I10,' elts ',I2,' el t s avt')

40 FORMAT('xyd ',I10,' suva',I4.4, '  ', I2,' ua va uva ')
41 FORMAT(' u', I1)
42 FORMAT(' v', I1)

50 FORMAT('xyd ',I10,' elts', I4.4, '  ', I2,' el avt avs')
51 FORMAT(' t', I1)
   RETURN
   END SUBROUTINE OUT_SMS_ONE
!==============================================================================|




