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
     if(AtmSym(iAtm).eq."H") then
       SkM = SkM + 1.0d0
     elseif(AtmSym(iAtm).eq."C") then
       SkM = SkM + 12.0d0
     elseif(AtmSym(iAtm).eq."N") then
       SkM = SkM + 14.0d0
     elseif(AtmSym(iAtm).eq."O") then
       SkM = SkM + 16.0d0
     elseif(AtmSym(iAtm).eq."CH") then
       SkM = SkM + 14.0d0
     elseif(AtmSym(iAtm).eq."Si") then
       SkM = SkM + 28.0d0
     elseif(AtmSym(iAtm).eq."Zn") then
       SkM = SkM + 65.4d0
     elseif(AtmSym(iAtm).eq."P") then
       SkM = SkM + 30.9d0
     elseif(AtmSym(iAtm).eq."Ca") then
       SkM = SkM + 40.0d0
     elseif(AtmSym(iAtm).eq."S") then
       SkM = SkM + 32.1d0
     elseif(AtmSym(iAtm).eq."F") then
       SkM = SkM + 19.0d0
     elseif(AtmSym(iAtm).eq."Na") then
       SkM = SkM + 23.0d0
     elseif(AtmSym(iAtm).eq."Cl") then
       SkM = SkM + 35.5d0
     elseif(AtmSym(iAtm).eq."Br") then
       SkM = SkM + 79.9d0
     elseif(AtmSym(iAtm).eq."Ga") then
       SkM = SkM + 69.7d0
     elseif(AtmSym(iAtm).eq."Fe") then
       SkM = SkM + 55.8d0
     elseif(AtmSym(iAtm).eq."Al") then
       SkM = SkM + 27.0d0
     elseif(AtmSym(iAtm).eq."B") then
       SkM = SkM + 10.8d0
     elseif(AtmSym(iAtm).eq."In") then
       SkM = SkM + 114.8d0
     elseif(AtmSym(iAtm).eq."Ti") then
       SkM = SkM + 47.9d0
     elseif(AtmSym(iAtm).eq."Sc") then
       SkM = SkM + 45.0d0
     elseif(AtmSym(iAtm).eq."Sn") then
       SkM = SkM + 118.7d0
     elseif(AtmSym(iAtm).eq."Ge") then
       SkM = SkM + 72.6d0
     elseif(AtmSym(iAtm).eq."Te") then
       SkM = SkM + 127.6d0
     elseif(AtmSym(iAtm).eq."ZZ") then
       SkM = SkM + 0.5d0
     elseif(AtmSym(iAtm).eq."AZ") then
       SkM = SkM + 105.0d0
     else
       write(6,'("Atom weight not defined for ",a2)') AtmSym(iAtm)
       stop
     end if
   end do

   SkD = SkM / SkV
   SkD = SkD * Fact

   return
 end
       
