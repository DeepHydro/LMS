!==============================================================================|
!     ENSURE DEPTH IS GREATER THAN MIN_DEPTH                                   |
!     IF THIS CONDITION IS VIOLATED, HALT THE PROGRAM AND WRITE A WARNING      |
!     MESSAGE TO THE SCREEN                                                    |
!==============================================================================|

   SUBROUTINE DEPTH_CHECK

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   INTEGER, DIMENSION(NPROCS) :: SBUF,RBUF
   INTEGER  :: I,II,MLOC,IERR
   REAL(SP) :: DMIN
!==============================================================================|

!--Calculate Minimum Depth and Set Global Node Number if Min Depth < MIN_DEPTH
   SBUF = 0
   MLOC = 0
   DMIN = MINVAL(D(1:M))
   IF(DMIN < MIN_DEPTH)THEN
    MLOC = MINLOC(D(1:M),DIM=1)
    ENDIF
!--Reduce in Master Processor Array and Dump To Screen
   SBUF(MYID) = MLOC
   RBUF = SBUF

!--If Depth Condition is Violated Write Warning and Halt User
   IF(MSR)THEN
     IF(SUM(RBUF) == 0)RETURN
     DO I=1,NPROCS
       II = RBUF(I)
       IF(II /= 0)THEN
         WRITE(*,*)'DEPTH IN NODE: ',II,' AT ',XG(II)+VXMIN,YG(II)+VYMIN, &
                        ' IS LESS THAN MIN_DEPTH. Minimum depth is:',DMIN
         WRITE(*,*)'ADJUST BATHYMETRY AT THIS (THESE) LOCATION(S) OR'
         WRITE(*,*)'RECOMPILE DOTFVM WITH FLOODING/DRYING FORMULATION'                 
         WRITE(*,*)'STOPPING....'
         CALL PSTOP
       END IF
     END DO
   END IF

   RETURN 
   END SUBROUTINE DEPTH_CHECK
!==============================================================================|
