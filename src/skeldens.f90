 SUBROUTINE Skel_Dens(dMesh,nP,IndCav,nAtm,AtmSym,SkV,SkM,SkD)

   IMPLICIT NONE

   INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

   INTEGER, INTENT(IN) :: nP, nAtm
   INTEGER, DIMENSION(:), INTENT(IN) :: IndCav
   INTEGER :: iP, iAtm

   REAL(DP), INTENT(IN) :: dMesh
   REAL(DP), INTENT(OUT) :: SkV, SkM, SkD
   REAL(DP), PARAMETER :: Fact=1.6606d0 ! Conversion from g mol-1 A-3 to g cm-3 = 1/0.6022
   REAL(DP) :: v

   CHARACTER(2), DIMENSION(:), INTENT(IN) :: AtmSym

   v = dMesh*dMesh*dMesh
   SkV = 0.d0
   SkM = 0.d0

! Include in the skeletal volume both filled and closed volume elements
   do iP = 1, nP
     if(IndCav(iP).eq.1.OR.IndCav(iP).eq.2) SkV = SkV + v
   end do
   
   do iAtm = 1, nAtm
     if(AtmSym(iAtm).eq."H")  SkM = SkM + 1.0d0
     if(AtmSym(iAtm).eq."C")  SkM = SkM + 12.0d0
     if(AtmSym(iAtm).eq."N")  SkM = SkM + 14.0d0
     if(AtmSym(iAtm).eq."O")  SkM = SkM + 16.0d0
     if(AtmSym(iAtm).eq."CH") SkM = SkM + 14.0d0
     if(AtmSym(iAtm).eq."Si") SkM = SkM + 28.0d0
     if(AtmSym(iAtm).eq."Ar")  SkM = SkM + 12.0d0  !for testing 
   end do

   SkD = SkM / SkV
   SkD = SkD * Fact

   return
 end
       
