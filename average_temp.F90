IMPLICIT NONE
INTEGER :: i,j,n,dummy
REAL*8,ALLOCATABLE :: v1(:),v2(:),KE(:),PE(:),TE(:),inst_T(:),T(:)
INTEGER,PARAMETER :: window=220
REAL*8 :: Tt, dummy1

OPEN(1,file='ener.txt')
OPEN(2,file='ENERGIES.txt')
CALL count(n)
ALLOCATE(v1(n),v2(n))
ALLOCATE(KE(n),PE(n),TE(n))
ALLOCATE(inst_T(n),T(n))
!------------------------------------------------------------------------
99 FORMAT(I6,2X,2F10.6,1X,3F10.6,2X,2F10.2)
WRITE(2,*)' MD_STEPS    vx       vy         KE        PE        TE       inst_Temp   Temp'
WRITE(2,*)'-------------------------------------------------------------------------------'
DO i = 1,n
  READ(1,*)dummy1,v1(i),v2(i),KE(i),PE(i),TE(i),inst_T(i)
END DO
DO i =1,n
     Tt = 0.0d0
     DO j = i,window+(i-1)
        Tt = inst_T(j) + Tt
     END DO
  T(i) = Tt/window
IF (i .le. (n-window))WRITE(2,99)i,v1(i),v2(i),KE(i),PE(i),TE(i),inst_T(i),T(i)
END DO
PRINT'(A,F10.2,A)',"SYSTEM TEMPERATURE   =",T(n-window),'K'
CLOSE(1,STATUS='DELETE') ; CLOSE(2)
END
!------------------------------------------------------------------------
SUBROUTINE count(nline)
INTEGER :: nline
CHARACTER*10 :: dummy
   DO
      READ(1,*,END=22) dummy
        nline = nline + 1
   END DO
22  CONTINUE
REWIND(1)
END SUBROUTINE
!------------------------------------------------------------------------

