!==============================================================================|
!   Write Binary Data For Iteration Number "IINTT"                             |
!==============================================================================|

   SUBROUTINE OUT_TXT(IINTT)   

!------------------------------------------------------------------------------|

   USE ALL_VARS

#  if defined (WET_DRY)
   USE MOD_WD
#  endif
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IINTT
   INTEGER :: I,K
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: UTMP,VTMP,WWTMP,KMTMP
   REAL(SP), ALLOCATABLE, DIMENSION(:) :: UATMP,VATMP
   REAL(SP), ALLOCATABLE, DIMENSION(:,:) :: T1TMP,S1TMP,R1TMP
   REAL(SP), ALLOCATABLE, DIMENSION(:)   :: ELTMP

#  if defined (DYE_RELEASE)
   REAL(SP), ALLOCATABLE, DIMENSION(:,:)   :: DYETMP
#  endif

   CHARACTER(LEN=4) :: FILENUMBER
   CHARACTER(LEN=120) :: DIR 
   CHARACTER(LEN=120) :: FNAME,FNAME1,FNAME2

!==============================================================================|
!------------------------------------------------------------------------------!
!  OPEN FILE (Name Based on Iteration Number)                                  !
!------------------------------------------------------------------------------!

   IF(MSR)THEN
     WRITE(FILENUMBER,'(I4.4)') IINTT
     FNAME1 = TRIM(CASENAME)//'_assim'//FILENUMBER//'.dat'
     FNAME2 = TRIM(CASENAME)//'_sim'//FILENUMBER//'.dat'
     DIR = TRIM(OUTDIR)//"/medm"
     FNAME = FNAME2 
!     OPEN(1,FILE=TRIM(DIR)//"/"//TRIM(FNAME),STATUS='unknown',FORM='unformatted') 
     OPEN(1,FILE=TRIM(DIR)//"/"//TRIM(FNAME),STATUS='unknown') 
     REWIND(1)
     WRITE(IPT,*)'DUMPING MEDM FILE: ',TRIM(FNAME)
   END IF

!------------------------------------------------------------------------------!
!  WRITE VALUES TO FILE (Single Processor Case)                                !
!------------------------------------------------------------------------------!

   IF(SERIAL)THEN

     !! ELEMENT BASED VALUES
     WRITE(1,*) IINTT,NGL,MGL,THOUR
#    if !defined (TWO_D_MODEL)
     DO I=1,NGL
       WRITE(1,10) (U(I,K),V(I,K),WW(I,K),KM1(I,K),K=1,KBM1)
     END DO
     
     !! NODE BASED VALUES
#       if !defined (DYE_RELEASE)     
     DO I=1,M
       WRITE(1,40) EL(I),(T1(I,K),S1(I,K),RHO1(I,K),K=1,KBM1)
     END DO
#       else
           DO I=1,M
              WRITE(1,50) EL(I),(T1(I,K),S1(I,K),RHO1(I,K),DYE(I,K),K=1,KBM1)
           END DO
#       endif   
#    else
     DO I=1,NGL
       WRITE(1,20) UA(I),VA(I)
     END DO
     
     !! NODE BASED VALUES
     DO I=1,M
       WRITE(1,30) EL(I)
     END DO
#    endif
   END IF


!
!--FORMATS---------------------------------------------------------------------!
!

10 FORMAT(12E14.5)
20 FORMAT(2E14.5)
30 FORMAT(E14.5)
40 FORMAT(13E14.5)
50 FORMAT(17E14.5)
   IF(MSR) CLOSE(1)

   RETURN
   END SUBROUTINE OUT_TXT
!==============================================================================|