
!==============================================================================|
!     CALCULATE THE GRADIENT OF THE WATER DEPTH WITH X ANDY                    |
!==============================================================================|

   SUBROUTINE DEPTH_GRADIENT      

!==============================================================================|
   USE ALL_VARS
   IMPLICIT NONE
   REAL(SP) :: PHPX,PHPY,F1
   INTEGER  :: I,J,I1,I2

!==============================================================================|

!----------CALCULATE DERIVATIVES OF DEPTH WITH X AND Y AT NODES----------------!

   DO I=1,M
     PHPX = 0.0_SP ; PHPY = 0.0_SP
     DO J=1,NTSN(I)-1
       I1 = NBSN(I,J)
       I2 = NBSN(I,J+1)
       F1 = 0.50_SP*(H(I1)+H(I2)) 

       PHPX = PHPX + F1*(VY(I1)-VY(I2))
       PHPY = PHPY + F1*(VX(I2)-VX(I1))

     END DO
     PHPX = PHPX/ART2(I)
     PHPY = PHPY/ART2(I)

     IF(PHPX==0.0_SP .AND.  PHPY ==0.0_SP)THEN
       SITA_GD(I) = 0.0_SP ; PHPN(I) = 0.0_SP
     ELSE
       SITA_GD(I) = ATAN2(PHPY,PHPX)
       PHPN(I)    = SQRT(PHPY*PHPY+PHPX*PHPX)
     END IF

   END DO

   RETURN
   END SUBROUTINE DEPTH_GRADIENT
!==============================================================================|
