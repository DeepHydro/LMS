MODULE MOD_WD
!#if defined (WET_DRY)

   USE MOD_PREC
   IMPLICIT NONE
   SAVE
!
!--Parameters for Wet/Dry Treatment                 
!
   LOGICAL  :: WET_DRY_ON                   !!TRY IF WET/DRY ACTIVE 

!-----variables controlling porosities through wet/dry determination----------------!
                                                                                                                          
   INTEGER , ALLOCATABLE :: ISWETN(:)           !!NODE POROSITY AT NODES FOR TIME N
   INTEGER , ALLOCATABLE :: ISWETC(:)           !!CELL POROSITY AT CELLS FOR TIME N
   INTEGER , ALLOCATABLE :: ISWETNT(:)          !!NODE POROSITY AT NODES FOR TIME N-1 INTERNAL
   INTEGER , ALLOCATABLE :: ISWETCT(:)          !!CELL POROSITY AT CELLS FOR TIME N-1 INTERNAL
   INTEGER , ALLOCATABLE :: ISWETCE(:)          !!CELL POROSITY AT CELLS FOR TIME N-1 EXTERNAL
   REAL(SP), ALLOCATABLE :: FWET_N_N(:)         !!NODE POROSITY AT NODES FOR TIME N
   REAL(SP), ALLOCATABLE :: FWET_C_C(:)         !!CELL POROSITY AT CELLS FOR TIME N
!!   REAL(SP), ALLOCATABLE :: UAS(:)          !!VERT-AVGD X-VELOC USED FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: VAS(:)          !!VERT-AVGD Y-VELOC USED FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: UARDS(:)        !!UA AVGD OVER EXT INT FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: VARDS(:)        !!VA AVGD OVER EXT INT FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: US(:,:)         !!X-VELOCITY FOR MASS CONS IN FLOODING/DRYING PROCESS
!!   REAL(SP), ALLOCATABLE :: VS(:,:)         !!Y-VELOCITY FOR MASS CONS IN FLOODING/DRYING PROCESS

   CONTAINS !------------------------------------------------------------------!
            ! SET_WD_PARAM        :   READ WET/DRY PARAMETERS FROM INPUT       !
            ! ALLOC_WD_DATA       :   ALLOCATE AND INITIALIZE WET/DRY ARRAYS   !
            ! SET_WD_DATA         :   SET VALUES IN WET/DRY ARRAYS             !
            ! WET_JUDGE           :   DETERMINE IF NODES/ELEMENTS ARE WET/DRY  ! 
            ! WD_UPDATE           :   SWAP WET/DRY VARIABLES BETWEEN TIME LEVS ! 
            ! WD_DUMP             :   DUMP WET/DRY FLAGS FOR RESTART           !
            ! WD_READ             :   READ WET/DRY FLAGS FOR RESTART           !
            ! -----------------------------------------------------------------!

!==============================================================================|
!==============================================================================|


!==============================================================================|
!==============================================================================|

                                                                                                                             
   SUBROUTINE SET_WD_DATA 
                                                                                                                             
!------------------------------------------------------------------------------|
!  INITIALIZE ARRAYS USED FOR WET/DRY TREATMENT                                |
!------------------------------------------------------------------------------|
                                                                                                                             
   USE ALL_VARS

   IMPLICIT NONE
   INTEGER :: I
                                                                                                                             
   IF(RESTART == 'cold_start') THEN

!-------- SET WET/DRY FLAGS AND MODIFY WATER SURFACE ELEVATION-----------------!

     CALL WET_JUDGE

!-------- EXCHANGE MODIFIED FREE SURFACE ELEVATION ACROSS PROCESSOR BOUNDS-----!

!-------- TRANSFER ELEVATION FIELD TO DEPTH AND OLD TIME LEVELS----------------!

     EL1 = ELF1
     D1  = H1 + EL1
     EL = ELF
     ET = EL
     D  = EL + H
     DT = D
     DTFA = D
     ET1 = EL1
     DT1 = D1

   END IF 

   RETURN
   END SUBROUTINE SET_WD_DATA

!==============================================================================|
!==============================================================================|

   SUBROUTINE ALLOC_WD_DATA  

!------------------------------------------------------------------------------|
!  ALLOCATE AND INITIALIZE WET/DRY TREATMENT ARRAYS                            |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   IMPLICIT NONE

!-----variables controlling porosities through wet/dry determination----------------!
                                                                                                                          
   ALLOCATE(ISWETN(0:MT))        ; ISWETN     = 1
   ALLOCATE(ISWETC(0:NT))        ; ISWETC     = 1
   ALLOCATE(ISWETNT(0:MT))       ; ISWETNT    = 1
   ALLOCATE(ISWETCT(0:NT))       ; ISWETCT    = 1
   ALLOCATE(ISWETCE(0:NT))       ; ISWETCE    = 1
   ALLOCATE(FWET_N_N(0:MT))      ; FWET_N_N   = 1.0_SP
   ALLOCATE(FWET_C_C(0:NT))      ; FWET_C_C   = 1.0_SP

!!   ALLOCATE(US(0:NT,KB))         ;US    = ZERO   !!X-VELOCITY FOR MASS CONSERVATION
!!   ALLOCATE(VS(0:NT,KB))         ;VS    = ZERO   !!Y-VELOCITY FOR MASS CONSERVATION

!!   ALLOCATE(UAS(0:NT))           ;UAS       = ZERO  !!VERT AVGD X-VELOC FOR MASS CONSERVATION
!!   ALLOCATE(VAS(0:NT))           ;VAS       = ZERO  !!VERT AVGD Y-VELOC FOR MASS CONSERVATION
!!   ALLOCATE(UARDS(0:NT))         ;UARDS     = ZERO  !!UA AVGD OVER EXTERNAL INT FOR MASS CONSERVATION
!!   ALLOCATE(VARDS(0:NT))         ;VARDS     = ZERO  !!VA AVGD OVER EXTERNAL INT FOR MASS CONSERVATION

   RETURN
   END SUBROUTINE ALLOC_WD_DATA

!==============================================================================|
!==============================================================================|

   SUBROUTINE WET_JUDGE

!------------------------------------------------------------------------------|
!  DETERMINE IF NODES/ELEMENTS ARE WET OR DRY                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   REAL(SP) :: DTMP
   INTEGER  :: ITA_TEMP
   INTEGER  :: I,IL,IA,IB,K1,K2,K3,K4,K5,K6

!
!--Determine If Node Points Are Wet/Dry Based on Depth Threshold---------------!
!
   ISWETN = 1
   DO I = 1, M
     DTMP = H(I) + ELF(I)
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETN(I) = 0
   END DO

!
!--Determine if Cells are Wet/Dry Based on Depth Threshold---------------------!
!
   ISWETC = 1
   DO I = 1, N
     DTMP =  MAX(ELF(NV(I,1)),ELF(NV(I,2)),ELF(NV(I,3)))  + &
             MIN(  H(NV(I,1)),  H(NV(I,2)),  H(NV(I,3)))
     IF((DTMP - MIN_DEPTH) < 1.0E-5_SP) ISWETC(I) = 0
   END DO

!
!--A Secondary Condition for Nodal Dryness-(All Elements Around Node Are Dry)--!
!
   DO I = 1, M
     IF(SUM(ISWETC(NBVE(I,1:NTVE(I)))) == 0)  ISWETN(I) = 0
   END DO

!
!--Adjust Water Surface So It Does Not Go Below Minimum Depth------------------!
!
   ELF = MAX(ELF,-H + MIN_DEPTH)

!
!--Recompute Element Based Depths----------------------------------------------!
!
   DO I = 1, N
     ELF1(I) = ONE_THIRD*(ELF(NV(I,1))+ELF(NV(I,2))+ELF(NV(I,3)))
   END DO

!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     FWET_N_N = FLOAT(ISWETN)
     FWET_C_C = FLOAT(ISWETC)
     CALL EXCHANGE(EC,NT,1,MYID,NPROCS,FWET_C_C)
     CALL EXCHANGE(NC,MT,1,MYID,NPROCS,FWET_N_N)
     ISWETN = INT(FWET_N_N+.5)
     ISWETC = INT(FWET_C_C+.5)
   END IF
#  endif

   RETURN
   END SUBROUTINE WET_JUDGE

!==============================================================================|
!==============================================================================|

   SUBROUTINE WD_UPDATE(INCASE)

!------------------------------------------------------------------------------|
!  SHIFT WET/DRY VARIABLES TO NEW TIME LEVELS                                  |
!------------------------------------------------------------------------------|

   USE MOD_PREC
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: INCASE
   INTEGER :: I


   SELECT CASE(INCASE)

!------------------------------------------------------------------------------!
   CASE(1)    !! SHIFT AT END OF EXTERNAL MODE
!------------------------------------------------------------------------------!
   ISWETCE=ISWETC
!------------------------------------------------------------------------------!
   CASE(2)    !! UPDATE NODE WET/DRY AFTER DEPTH ADJUSTMENT
!------------------------------------------------------------------------------!
   DO I = 1,M
     IF(DTFA(I)-MIN_DEPTH <= 1.0E-5_SP) THEN
       ISWETN(I) = 0
     END IF
   END DO
!------------------------------------------------------------------------------!
   CASE(3)    !! SHIFT VARIABLES AT END OF INTERNAL MODE
!------------------------------------------------------------------------------!

   ISWETCT=ISWETC
   ISWETNT=ISWETN

   END SELECT

   RETURN
   END SUBROUTINE WD_UPDATE

!==============================================================================|
!==============================================================================|

!==============================================================================|
!   DUMP WET/DRY FLAG DATA FOR RESTART                                         |
!==============================================================================|

   SUBROUTINE WD_DUMP(FNAME)

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
   INTEGER I
   CHARACTER(LEN=120) :: FNAME
!==============================================================================|

   IF(MSR)THEN
     OPEN(1,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')

     REWIND(1)
     WRITE(1,*) IINT
     WRITE(1,*) NGL,MGL
   END IF

   IF(SERIAL)THEN
     WRITE(1,*) (ISWETC(I), I=1,N)
     WRITE(1,*) (ISWETN(I), I=1,M)
   ELSE
   ALLOCATE(NTEMP1(NGL),NTEMP2(MGL))
#  if defined (MULTIPROCESSOR)
   CALL IGATHER(LBOUND(ISWETC,1),UBOUND(ISWETC,1),N,NGL,1,MYID,NPROCS,EMAP,ISWETC,NTEMP1)
   CALL IGATHER(LBOUND(ISWETN,1),UBOUND(ISWETN,1),M,MGL,1,MYID,NPROCS,NMAP,ISWETN,NTEMP2)
   IF(MSR)THEN
     WRITE(1,*) (NTEMP1(I), I=1,NGL)
     WRITE(1,*) (NTEMP2(I), I=1,MGL)
   END IF
   DEALLOCATE(NTEMP1,NTEMP2)
#  endif
   END IF

   IF(MSR) CLOSE(1)

   RETURN
   END SUBROUTINE WD_DUMP
!==============================================================================|


!==============================================================================|
!   READ WET/DRY FLAG DATA FOR RESTART                                         |
!==============================================================================|

   SUBROUTINE WD_READ(FNAME)

!------------------------------------------------------------------------------|

   USE ALL_VARS
#  if defined (MULTIPROCESSOR)
   USE MOD_PAR
#  endif
   IMPLICIT NONE
   INTEGER, ALLOCATABLE,DIMENSION(:) :: NTEMP1,NTEMP2
   INTEGER I,IINT_TMP
   CHARACTER(LEN=120) :: FNAME
   LOGICAL :: FEXIST
!==============================================================================|


!   FNAME = "./"//TRIM(INPDIR)//"/"//trim(casename)//"_restart_wd.dat"
!
!--Make Sure Wet/Dry Flag Data Exists------------------------------------------!
!
   INQUIRE(FILE=TRIM(FNAME),EXIST=FEXIST)
   IF(MSR .AND. .NOT.FEXIST)THEN
     WRITE(IPT,*)'WET/DRY RESTART FILE: ',FNAME,' DOES NOT EXIST'
     WRITE(IPT,*)'HALTING.....'
     CALL PSTOP
   END IF

!
!--Ensure File Header is Consistent with Main Flow Variable Restart File-------!
!
   OPEN(1,FILE=FNAME,FORM='FORMATTED',STATUS='UNKNOWN')
   REWIND(1)
   READ(1,*) IINT_TMP
   READ(1,*)
   IF(IINT_TMP /= IINT .AND. MSR)THEN
     WRITE(IPT,*)'IINT IN ',FNAME,' NOT EQUAL TO IINT'
     WRITE(IPT,*)'FROM MAIN RESTART FILE',IINT,IINT_TMP
     CALL PSTOP
   END IF

!
!--Read Variables--------------------------------------------------------------!
!
   ALLOCATE(NTEMP1(NGL),NTEMP2(MGL))
   READ(1,*) (NTEMP1(I), I=1,NGL)
   READ(1,*) (NTEMP2(I), I=1,MGL)

!
!--Transfer Variables to Local Domains-----------------------------------------!
!
   IF(SERIAL)THEN
     ISWETC(1:N) = NTEMP1(1:N)
     ISWETN(1:M) = NTEMP2(1:M)
   END IF

#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     DO I=1,N
       ISWETC(I) = NTEMP1(EGID(I))
     END DO
     DO I=1,M
       ISWETN(I) = NTEMP2(NGID(I))
     END DO
   END IF
#  endif

   DEALLOCATE(NTEMP1,NTEMP2)
   CLOSE(1)
!
!--Extend Element/Node Based Wet/Dry Flags to Domain Halo----------------------!
!
#  if defined (MULTIPROCESSOR)
   IF(PAR)THEN
     FWET_C_C = ISWETC
     FWET_N_N = ISWETN
     CALL EXCHANGE(EC,NT,1,MYID,NPROCS,FWET_C_C)
     CALL EXCHANGE(NC,MT,1,MYID,NPROCS,FWET_N_N)
     ISWETN = INT(FWET_N_N+.5)
     ISWETC = INT(FWET_C_C+.5)
   END IF
#  endif

   ISWETNT = ISWETN
   ISWETCT = ISWETC
   ISWETCE = ISWETC
!
!--Extend Element/Node Based Wet/Dry Flags to Previous Time Levels-------------!
!

   RETURN
   END SUBROUTINE WD_READ
!==============================================================================|


!# endif
END MODULE MOD_WD
