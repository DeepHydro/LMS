
!==============================================================================|
!     READ MESH AND DETERMINE NUMBER OF ELEMENTS AND NODES                     |
!==============================================================================|

   SUBROUTINE GETDIM

!==============================================================================|
   USE ALL_VARS
   USE MOD_OBCS
   IMPLICIT NONE
   INTEGER, ALLOCATABLE :: NVTMP(:,:)
   INTEGER              :: I,LM1,J,IOS,LMAX
!==============================================================================|

   IF(MSR)WRITE(IPT,*)'!                  GLOBAL INFORMATION                            !'
   IF(MSR)WRITE(IPT,*)'!                                                                !'
   IF(MSR)WRITE(IPT,*)'!  # OF PROCESSORS       :',NPROCS
!
!----------------Determine Number of Elements----------------------------------!
!
   LM1 = 1
   DO WHILE(.TRUE.)
     READ(INGRD,*,IOSTAT=IOS)J
     IF(IOS < 0)THEN
       IF(MSR)WRITE(IPT,*)'ERROR DETERMINING GRID DIMENSIONS'
       CALL PSTOP
     END IF
     IF(J == 1 .AND. LM1 /= 1)THEN
       NGL = LM1
       EXIT
     END IF
     LM1 = J
   END DO
   REWIND(INGRD)

!
!----------------Determine Number of Nodes-------------------------------------!
!
   ALLOCATE(NVTMP(NGL,3))
   DO I=1,NGL
     READ(INGRD,*)J,NVTMP(I,1),NVTMP(I,2),NVTMP(I,3)
   END DO
   MGL = MAX(MAXVAL(NVTMP(:,1)) , MAXVAL(NVTMP(:,2)) , MAXVAL(NVTMP(:,3)))
   DEALLOCATE(NVTMP)
   REWIND(INGRD)

!
!----------------Determine Number of Groundwater Flux Points ------------------!
!
   READ(INBFW,*)
   READ(INBFW,*) IBFW_GL
   REWIND(INBFW)
!
!----------------Determine Number of Freshwater Flux Points--------------------!
!
   READ(INRIV,*)
   READ(INRIV,*) NUMQBC_GL
   REWIND(INRIV)
   
!
!------Determine Number of Frictional Geostrophic Correction Inflow Nodes------!
!
   NOBCGEO_GL = 0
   IF(JMPOBC) READ(INJMP,*)NOBCGEO_GL
   REWIND(INJMP)
!
!----------------Report--------------------------------------------------------!
!
   IF(MSR)WRITE(IPT,*)'!  # OF NODES            :',MGL
   IF(MSR)WRITE(IPT,*)'!  # OF ELEMENTS         :',NGL
   IF(MSR)WRITE(IPT,*)'!  # OF OPEN BNDRY NODES :',IOBCN_GL
   IF(MSR)WRITE(IPT,*)'!  # OF BOTTOM FLUX PTS  :',IBFW_GL
   IF(MSR)WRITE(IPT,*)'!  # OF FRESH WATER PTS  :',NUMQBC_GL
   IF(MSR)WRITE(IPT,*)'!  # OF JMP OBC NODES    :',NOBCGEO_GL

!
!-------------Set Local Number of Elements/Nodes for Serial Case---------------!
!
   N = NGL
   M = MGL
   NT = N
   MT = M

!
!-------------Read In Global Node Numbering------------------------------------!
!
   ALLOCATE(NVG(0:NGL,4))
   DO I=1,NGL
     READ(INGRD,*) J,NVG(I,1),NVG(I,3),NVG(I,2)
     NVG(I,4) = NVG(I,1)
   END DO
   IF(MSR)WRITE(IPT,*)  '!  GLOBAL MESH READING   :    FINISHED'

   !==============================================================================|
!  GENERATE LOCAL NODE CONNECTIVITY (NV) FROM GLOBAL NODE CONNECTIVITY (NVG)   |
!  USING LOCAL TO GLOBAL MAPPING FOR INTERIOR ELEMENTS (EGID)                  |
!  AND LOCAL TO GLOBAL MAPPING FOR HALO ELEMENTS (HE_LST)                      |
!==============================================================================|

   IF(SERIAL) NV = NVG

!==============================================================================|
!   SET UP LOCAL MESH (HORIZONTAL COORDINATES)                                 |
!==============================================================================|


!--------------READ IN X AND Y GLOBAL COORDINATES AT NODES---------------------!

   ALLOCATE(XG(0:MGL),YG(0:MGL)) ; XG = 0.0_SP ; YG = 0.0_SP
   ALLOCATE(HG(0:MGL))  ; HG = 0.0_SP
   
   DO I=1,MGL
     READ(INGRD,*)J,XG(I),YG(I),HG(I)
   END DO
   CLOSE(INGRD)
   
   RETURN
   END SUBROUTINE GETDIM
!==============================================================================|
