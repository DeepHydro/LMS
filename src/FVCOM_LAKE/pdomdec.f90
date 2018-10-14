!==============================================================================|
!  SET UP LOCAL PHYSICAL DOMAIN (CONNECTIVITY/MESH)                            |
!==============================================================================|

   SUBROUTINE PDOMDEC

!==============================================================================!
   USE ALL_VARS

   IMPLICIT NONE
   INTEGER I,EGL,J,IERR,I1,I2,N_SPONGE
   REAL(SP), ALLOCATABLE :: CORRG(:),CORR(:)
   REAL(SP), ALLOCATABLE :: R_SPG(:),C_SPG(:) 
   INTEGER, ALLOCATABLE  :: N_SPG(:)
   REAL(SP)  TEMP,DTMP,C_SPONGE

!!==============================================================================|
!!  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
!!  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
!!  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
!!==============================================================================|
!
   IF(SERIAL) NV = NVG
!
!!==============================================================================|
!!   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
!!==============================================================================|
!

   VXMIN = MINVAL(XG(1:MGL)) ; VXMAX = MAXVAL(XG(1:MGL))
   VYMIN = MINVAL(YG(1:MGL)) ; VYMAX = MAXVAL(YG(1:MGL))

!--------------SHIFT GRID TO UPPER RIGHT CARTESIAN-----------------------------!

   XG = XG - VXMIN
   YG = YG - VYMIN
   XG(0) = 0.0_SP ; YG(0) = 0.0_SP

!--------------CALCULATE GLOBAL ELEMENT CENTER GRID COORDINATES----------------!

   ALLOCATE(XCG(0:NGL),YCG(0:NGL)) ; XCG = 0.0_SP ; YCG = 0.0_SP
   DO I=1,NGL   
     XCG(I)  = (XG(NVG(I,1)) + XG(NVG(I,2)) + XG(NVG(I,3)))/3.0_SP
     YCG(I)  = (YG(NVG(I,1)) + YG(NVG(I,2)) + YG(NVG(I,3)))/3.0_SP
   END DO

     XCG(0) = 0.0_SP ; YCG(0) = 0.0_SP

!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL)THEN
       VX = XG
       VY = YG
     END IF

!==============================================================================|
!   SET UP LOCAL MESH (BATHYMETRIC DEPTH)                                      |
!==============================================================================|


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!

     IF(SERIAL) H = HG

!--------------CALCULATE EXTREMUMS---------------------------------------------!

     HMAX = MAXVAL(ABS(HG(1:MGL)))
     HMIN = MINVAL(HG(1:MGL))

!==============================================================================|
!   SET UP LOCAL CORIOLIS FORCE                                                |
!==============================================================================|

!--------------READ IN CORIOLIS PARAMETER--------------------------------------!

     ALLOCATE(CORRG(0:MGL))  ; CORRG = 0.0_SP
       DO I=1,MGL
         READ(INCOR,*) TEMP,TEMP,CORRG(I)
       END DO
     CLOSE(INCOR)


!--------------TRANSFORM TO LOCAL DOMAINS IF PARALLEL--------------------------!
     ALLOCATE(CORR(0:MT)) ; CORR = 0.0_SP
     IF(SERIAL) CORR = CORRG

!==============================================================================|
!   COMPUTE FACE CENTER VALUES FOR GRID, DEPTH, AND CORIOLIS PARAMETER         |
!==============================================================================|

     DO I=1,NT
!       XC(I)  = SUM(VX(NV(I,1:3)))/3.0
       XC(I)  = (VX(NV(I,1)) + VX(NV(I,2)) + VX(NV(I,3)))/3.0_SP
       YC(I)  = (VY(NV(I,1)) + VY(NV(I,2)) + VY(NV(I,3)))/3.0_SP
!       YC(I)  = SUM(VY(NV(I,1:3)))/3.0
       H1(I)  = SUM( H(NV(I,1:3)))/3.0_SP
       COR(I) = CORR(NV(I,1)) + CORR(NV(I,2)) + CORR(NV(I,3))
       COR(I) = COR(I)/3.0_SP
!       COR(I) = SUM(CORR(NV(I,1:3)))/3.0
       COR(I) = 2.*7.292e-5_SP*SIN(COR(I)*2.0_SP*3.14159_SP/360.0_SP)
     END DO

!==============================================================================|
!   COMPUTE GRAVITY VARIED WITH LATITUDE                                       |
!==============================================================================|

     ALLOCATE(GRAV_N(0:MT),GRAV_E(0:NT))

     GRAV_N = GRAV
     GRAV_E = GRAV

!==============================================================================|
!   WRITE TO SMS GRID FILE WHILE GLOBAL VALUES EXIST                           |
!==============================================================================|

   IF(MSR)THEN
     WRITE(IOSMSD,*)'scat2d'
     WRITE(IOSMSD,*)'xyd ',MGL,' dep ',1,' dep '
     DO I=1,MGL
       WRITE(IOSMSD,*) XG(I),YG(I),HG(I)
     END DO
     CLOSE(IOSMSD)
   END IF
   DEALLOCATE(CORR,CORRG)

   RETURN
   END SUBROUTINE PDOMDEC
!==============================================================================|
