!==============================================================================|
!   DUMP DATA FILE FOR RESTART                                                 |
!==============================================================================|

   SUBROUTINE ARCRST            
   USE ALL_VARS

#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (DYE_RELEASE)
   USE MOD_DYE
#  endif   
   IMPLICIT NONE
   REAL(SP), ALLOCATABLE :: RTP(:)
   INTEGER I,K,ME,NPC,N1
   CHARACTER(LEN=6) :: FILENUMBER
   CHARACTER(LEN=120) :: DIR,FNAME
   CHARACTER(LEN=8) :: CH8
!==============================================================================|
   
   ME = MYID ; NPC = NPROCS 
   ALLOCATE(RTP(N)) ; RTP = 0.0_SP

   IF(MSR)THEN
     if (IRESTART .NE. 0) then
        WRITE(FILENUMBER,'(I6.6)') IINT/IRESTART
        OPEN(1,FILE='re_'//FILENUMBER//'.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')      
     else
        WRITE(FILENUMBER,'(I6.6)') IINT
        OPEN(1,FILE='re_'//FILENUMBER//'.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')      
     end if
     REWIND(1)
     WRITE(1) IINT
   END IF

   IF(SERIAL)THEN
     WRITE(1) ((U(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((V(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((W(I,K),    K=1,KB),I=0,N)
#    if defined (GOTM)
     WRITE(1) ((TKE(I,K),   K=1,KB),I=0,M)
     WRITE(1) ((TEPS(I,K),  K=1,KB),I=0,M)
#    else
     WRITE(1) ((Q2(I,K),   K=1,KB),I=0,M)
     WRITE(1) ((Q2L(I,K),  K=1,KB),I=0,M)
     WRITE(1) ((L(I,K)  ,  K=1,KB),I=0,M)
#    endif
     WRITE(1) ((S(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((T(I,K),    K=1,KB),I=0,N)
     WRITE(1) ((RHO(I,K),  K=1,KB),I=0,N)
     WRITE(1) ((TMEAN(I,K),K=1,KB),I=0,N)
     WRITE(1) ((SMEAN(I,K),K=1,KB),I=0,N)
     WRITE(1) ((RMEAN(I,K),K=1,KB),I=0,N)

     WRITE(1) ((S1(I,K),    K=1,KB),I=1,M)
     WRITE(1) ((T1(I,K),    K=1,KB),I=1,M)
     WRITE(1) ((RHO1(I,K),  K=1,KB),I=1,M)
     WRITE(1) ((TMEAN1(I,K),K=1,KB),I=1,M)
     WRITE(1) ((SMEAN1(I,K),K=1,KB),I=1,M)
     WRITE(1) ((RMEAN1(I,K),K=1,KB),I=1,M)

     WRITE(1) ((KM(I,K),K=1,KB),I=1,M)
     WRITE(1) ((KH(I,K),K=1,KB),I=1,M)
     WRITE(1) ((KQ(I,K),K=1,KB),I=1,M)

     WRITE(1) (UA(I), I=0,N)
     WRITE(1) (VA(I), I=0,N)

     WRITE(1) (EL1(I), I=1,N)
     WRITE(1) (ET1(I), I=1,N)
     WRITE(1) (H1(I),  I=1,N)
     WRITE(1) (D1(I),  I=1,N)
     WRITE(1) (DT1(I), I=1,N)
     WRITE(1) (RTP(I), I=1,N)

     WRITE(1) (EL(I), I=1,M)
     WRITE(1) (ET(I), I=1,M)
     WRITE(1) (H(I),  I=1,M)
     WRITE(1) (D(I),  I=1,M)
     WRITE(1) (DT(I), I=1,M)

#    if defined (WATER_QUALITY)
     DO N1 = 1, NB
       WRITE(1) ((WQM(I,K,N1),K=1,KB),I=1,M)
     END DO
#    endif

#    if defined (DYE_RELEASE)
     IF(IINT.GT.IINT_SPE_DYE_B) THEN
     WRITE(1) ((DYE(I,K),K=1,KB),I=1,M)
     WRITE(1) ((DYEMEAN(I,K),K=1,KB),I=1,M)
     ENDIF
#    endif

   END IF
   IF(MSR) CLOSE(1)
   DEALLOCATE(RTP)

   RETURN
   END SUBROUTINE ARCRST
!==============================================================================|



