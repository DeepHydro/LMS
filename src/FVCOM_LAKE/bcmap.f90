
   SUBROUTINE BCMAP

!==============================================================================!
   USE ALL_VARS
   USE MOD_OBCS

   IMPLICIT NONE

   INTEGER              :: I,I1,I2,NCNT,IERR,J
   INTEGER, ALLOCATABLE :: TEMP1(:),TEMP2(:),TEMP3(:),TEMP4(:),&
                           TEMP5(:),TEMP6(:),TEMP7(:),ITEMP(:)
!------------------------------------------------------------------------------!

!==============================================================================|
!   OPEN BOUNDARY CONDITION NODES                                              |
!==============================================================================|

!----------------------------REPORT--------------------------------------------!

   IF(MSR)WRITE(IPT,*  )'!'
   IF(MSR)WRITE(IPT,*)'!           SETTING UP OPEN BOUNDARY NODES  '
   IF(MSR)WRITE(IPT,*  )'!'

   IOBCN = 0
   IBCN  = 0

   IF(IOBCN_GL > 0)THEN

!------------Read in Open Boundary Nodes and Temperature/Salinity Conditions---!

     ALLOCATE(I_OBC_GL(IOBCN_GL))
     ALLOCATE(TYPE_OBC_GL(IOBCN_GL))
!GWC ADD OBC FORCING FOR TEMP/SALT
     ALLOCATE(TEMP_OBC_GL(IOBCN_GL)) ; TEMP_OBC_GL = 0.
     ALLOCATE(SALT_OBC_GL(IOBCN_GL)) ; SALT_OBC_GL = 0.
!GWC
     READ(INOBC,*)

     IF(TS_NUDGING_OBC)THEN
!      GWC READ IN SALT/TEMP OBC FORCING
       DO I=1,IOBCN_GL
         READ(INOBC,*) I1,I_OBC_GL(I),TYPE_OBC_GL(I),TEMP_OBC_GL(I),SALT_OBC_GL(I) 
       END DO
!      GWC
     ELSE
       DO I=1,IOBCN_GL
         READ(INOBC,*) I1,I_OBC_GL(I),TYPE_OBC_GL(I) 
       END DO
     END IF

     CLOSE(INOBC)

!----------------------Make Sure It Is In Global Domain------------------------!

     DO I=1,IOBCN_GL
       IF((I_OBC_GL(I) > MGL))THEN
         WRITE(IPT,*)'==================ERROR=================================='
         WRITE(IPT,*)'OPEN BOUNDARY NODE NUMBER',I,'IS NOT IN THE'
         WRITE(IPT,*)'GLOBAL DOMAIN'
         WRITE(IPT,*)'CHECK INPUT FILE AND ENSURE OPEN BOUNDARY NODES <= ',MGL
         WRITE(IPT,*)'========================================================='
         CALL PSTOP
       END IF
     END DO

!-------------------Ensure OBC Types are Valid-[1-->10]------------------------!

     DO I=1,IOBCN_GL
       IF(TYPE_OBC_GL(I) < 1 .OR. TYPE_OBC_GL(I) > 10)THEN
         IF(MSR)THEN
           WRITE(*,*)'ERROR: Outer Boundary Node Type in File:'
           WRITE(*,*)'casename_obc.dat Must be >0 and <= 11'
         END IF
         CALL PSTOP
       END IF
     END DO

!----------Shift Open Boundary Node List,Type,Salt,and Temp to Local-----------!

     IF(SERIAL)THEN

       IOBCN    = IOBCN_GL

       ALLOCATE(I_OBC_N(IOBCN))
       I_OBC_N = I_OBC_GL
       ALLOCATE(TYPE_OBC(IOBCN))
       TYPE_OBC = TYPE_OBC_GL
!      GWC ADD TEMP/SALT
       ALLOCATE(TEMP_OBC(IOBCN)) ; TEMP_OBC = 0.0
       TEMP_OBC = TEMP_OBC_GL
       ALLOCATE(SALT_OBC(IOBCN)) ; SALT_OBC = 0.0
       SALT_OBC = SALT_OBC_GL
!      GWC
     END IF



!----------------------Set 11 Types Open Boundary Nodes Arrays------------------!

     CALL SEPARATE_OBC


   END IF !!IOBCN_GL > 0

!==============================================================================|
!   NODES USED TO CORRECT INFLOW FOR FRICTIONALLY ADJUSTED GEOSTROPHIC FLOW    |
!==============================================================================|

!----------------------Read In Nodes-------------------------------------------!
   NOBCGEO_GL = 0 ; NOBCGEO = 0
   IF(JMPOBC)THEN
     READ(INJMP,*)NOBCGEO_GL
     ALLOCATE(IBCGEO(NOBCGEO_GL))
     DO I=1,NOBCGEO_GL
       READ(INJMP,*)IBCGEO(I)
     END DO
     CLOSE(INJMP)

!----------------------Make Sure It Is In Global Domain------------------------!
     DO I=1,NOBCGEO_GL
       IF(IBCGEO(I) > MGL)THEN
         WRITE(IPT,*)'==================ERROR=================================='
         WRITE(IPT,*)'JMP BOUNDARY NODE NUMBER',I,'IS NOT IN THE'
         WRITE(IPT,*)'GLOBAL DOMAIN'
         WRITE(IPT,*)'CHECK INPUT FILE AND ENSURE JMP NODES <= ',MGL
         WRITE(IPT,*)'========================================================='
         CALL PSTOP
       END IF
     END DO

!----------------------Shift To Local Domain If Parallel-----------------------!


     IF(SERIAL) NOBCGEO = NOBCGEO_GL


#   if defined (MULTIPROCESSOR)
     IF(PAR)THEN
       ALLOCATE(TEMP1(NOBCGEO_GL))
       NCNT = 0
       DO I=1,NOBCGEO_GL
         I1 = NLID_X(IBCGEO(I))
         IF(I1 /= 0)THEN
           NCNT = NCNT + 1
           TEMP1(NCNT) = I1
         END IF
       END DO

       NOBCGEO = NCNT
       DEALLOCATE(IBCGEO)
       IF(NOBCGEO > 0)THEN
         ALLOCATE(IBCGEO(NOBCGEO))
         IBCGEO = TEMP1(1:NCNT)
       END IF
       DEALLOCATE(TEMP1)
     END IF
#   endif

   END IF  !!JMPOBC = .TRUE.

!==============================================================================|
!   REPORT AND CHECK RESULTS                                                   |
!==============================================================================|
   ALLOCATE(TEMP1(NPROCS),TEMP2(NPROCS),TEMP3(NPROCS),TEMP4(NPROCS))
   ALLOCATE(TEMP5(NPROCS),TEMP6(NPROCS),TEMP7(NPROCS))
   TEMP1(1)  = IOBCN
   TEMP2(1) = IBCN(1)
   TEMP3(1) = IBCN(2)
   TEMP4(1) = IBCN(3)
   TEMP5(1) = IBCN(4)
   TEMP6(1) = IBCN(5)
   TEMP7(1) = NOBCGEO

# if defined (MULTIPROCESSOR)
   CALL MPI_GATHER(IOBCN,  1,MPI_INTEGER,TEMP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(IBCN(1),1,MPI_INTEGER,TEMP2,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(IBCN(2),1,MPI_INTEGER,TEMP3,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(IBCN(3),1,MPI_INTEGER,TEMP4,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(IBCN(4),1,MPI_INTEGER,TEMP5,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(IBCN(5),1,MPI_INTEGER,TEMP6,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(NOBCGEO,1,MPI_INTEGER,TEMP7,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
# endif

   
   IF(MSR)WRITE(IPT,100)'!  IOBCN                 :',IOBCN_GL,   (TEMP1(I),I=1,NPROCS)
   IF(MSR)WRITE(IPT,100)'!  IBCN(1)               :',IBCN_GL(1), (TEMP2(I),I=1,NPROCS)
   IF(MSR)WRITE(IPT,100)'!  IBCN(2)               :',IBCN_GL(2), (TEMP3(I),I=1,NPROCS)
   IF(MSR)WRITE(IPT,100)'!  IBCN(3)               :',IBCN_GL(3), (TEMP4(I),I=1,NPROCS)
   IF(MSR)WRITE(IPT,100)'!  IBCN(4)               :',IBCN_GL(4), (TEMP5(I),I=1,NPROCS)
   IF(MSR)WRITE(IPT,100)'!  IBCN(5)               :',IBCN_GL(5), (TEMP6(I),I=1,NPROCS)
   IF(MSR)WRITE(IPT,100)'!  NOBCGEO               :',NOBCGEO_GL, (TEMP7(I),I=1,NPROCS)
   DEALLOCATE(TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7)

   RETURN
   100 FORMAT(1X,A26,I6," =>",2X,4(I5,1H,))
   END SUBROUTINE BCMAP
!==============================================================================|
