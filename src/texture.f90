 SUBROUTINE Texture(nP,dMesh,nVol,DiamStep,VMinD,Cumulative_VMinD,UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV,nR,iM1,iM2, &     
                   surf_computation,DistMin,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,Rad,IndCav,Surf,Accessible,IndSurf)

 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
 INTEGER, INTENT(IN) :: nVol, surf_computation, nP, iM1, iM2
 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCav
 INTEGER, DIMENSION(:), INTENT(IN) :: IndSurf, nR
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: VMinD
 REAL(DP), INTENT(IN) :: DiamStep, dMesh, Rad
 REAL(DP), INTENT(OUT) :: UltraV, MicroV, SmallMesoV, LargeMesoV, MacroV, TotPorV
 REAL(DP), INTENT(OUT) :: UltraS, MicroS, SmallMesoS, LargeMesoS, MacroS
 REAL(DP), DIMENSION(:), INTENT(IN) :: DistMin
 REAL(DP), DIMENSION(:), INTENT(OUT) :: Cumulative_VMinD
 REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: Surf
 LOGICAL, INTENT(IN) :: Accessible
 INTEGER :: iVol, iP, MinD, iC, iPNew, iX, iY, iZ, lCube
 REAL(DP) :: NAV, v_block, D, Dist_X, Dist_Y, Dist_Z, Rad1, XP, YP, ZP, XP1, YP1, ZP1
 REAL(DP), PARAMETER :: Zero=0.0d0, Two=2.0d0
 REAL(DP), PARAMETER :: UltraMax=7.0d0, MicroMax=20.0d0, SmallMesoMax=35.0d0, LargeMesoMax=50.0d0
 LOGICAL :: Overlap
 INTERFACE
   FUNCTION MoveCub(iP,iX,iY,iZ,nR)
     INTEGER, INTENT(IN) :: iP,iX,iY,iZ
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER :: MoveCub
   END FUNCTION MoveCub
 END INTERFACE

 NAV = Zero
 v_block = dMesh*dMesh*dMesh
 Rad1 = Rad+(dMesh*0.1)

! Assign the block to the suitable pore set
! Volumes associated to RMin are not added to VMinD
! At this point is the total volume (accessible and not accessible)
 do iP=1,nP
   if(IndCav(iP).ne.0) cycle
   MinD = INT(DistMin(iP)/DiamStep) + 1
   VMinD(MinD) = VMinD(MinD) + v_block 

   !Compute the Accessible volume, if required
   if(Accessible) then
     ! Find the block Cartesian coordinates
     iC = iP - 1
     ZP = INT(iC/(iM2))*dMesh + dMesh/2.0d0
     YP = INT(MOD(iC,iM2)/nR(1))*dMesh + dMesh/2.0d0
     XP = MOD(MOD(iC,iM2),nR(1))*dMesh + dMesh/2.0d0
     ! Check all the blocks that could fall inside Rad or Rad1. 
     Overlap = .False.  !If TRUE this probe overlaps to the wall (then it will be considered Not Accessible)

     lCube = nint(Rad1 / dMesh) + 1
     do iX = -lCube,lCube
       if(Overlap) exit                                                                                                                               
       do iY = -lCube,lCube
         if(Overlap) exit
         do iZ = -lCube,lCube
           if(Overlap) exit
           iPNew = MoveCub(iP, iX, iY, iZ, nR)
           ! Skip if the new block is void (IndSurf=0) or has already been classified as occupiable (IndSurf=3)
           if (IndCav(iPNew).eq.0.or.IndCav(iPNew).eq.3) cycle
           ! Find the coordinates of the new block
           iC = iPNew - 1
           ZP1 = INT(iC/(iM2))*dMesh + dMesh/2.0d0
           YP1 = INT(MOD(iC,iM2)/nR(1))*dMesh + dMesh/2.0d0
           XP1 = MOD(MOD(iC,iM2),nR(1))*dMesh + dMesh/2.0d0

           ! Compute the distance between blocks along the coordinates 
           Dist_X = abs(XP - XP1) 
           Dist_Y = abs(YP - YP1) 
           Dist_Z = abs(ZP - ZP1) 
           ! If the checked block fell outside the cell, it was converted to its periodic image inside the cell
           ! In this case, the distance results unexpectedly large, and it is recomputed correctly
           if (Dist_X.gt.2.d0*Rad1) Dist_X = nR(1)*dMesh - Dist_X
           if (Dist_Y.gt.2.d0*Rad1) Dist_Y = nR(2)*dMesh - Dist_Y
           if (Dist_Z.gt.2.d0*Rad1) Dist_Z = nR(3)*dMesh - Dist_Z
           ! Compute the actual distance
           D = sqrt(Dist_X*Dist_X + Dist_Y*Dist_Y + Dist_Z*Dist_Z)
         
           ! If D is lower than Rad the block overlaps with the probe, and it will be discarded later
           if (D.le.Rad1) then
              Overlap = .True.
              IndCav(iP) = 3
              VMinD(MinD) = VMinD(MinD) - v_block 
           end if

         end do
       end do
     end do

   end if
 end do

! Compute simplified, total cumulative volumes and surface
 TotPorV = Zero
 UltraV = Zero
 MicroV = Zero
 SmallMesoV = Zero
 LargeMesoV = Zero
 MacroV = Zero

 UltraS = Zero
 MicroS = Zero
 SmallMesoS = Zero
 LargeMesoS = Zero
 MacroS = Zero
 
 Allocate(Surf(nVol))
 Surf = 0.0d0

! Save the nature of the block, through the vector IndCav(iP)
! IndCav , Nature of block
!   1        filled
!   2        excluded because too small
!   3        Non Accessible Volume 
!   4        surf ultramicro (Rmin < 7 A)
!   5        surf micro (7 < Rmin < 20 A)
!   6        surf small meso (20 < Rmin < 35 A)
!   7        surf large meso (35 < Rmin < 50 A)
!   8        surf macro (Rmin > 50 A)
!   9        vol ultramicro (Rmin < 7 A)
!   10       vol micro (7 < Rmin < 20 A)
!   11       vol small meso (20 < Rmin < 35 A)
!   12       vol large meso (35 < Rmin < 50 A)
!   13       vol macro (Rmin > 50 A)           

 do iP= 1, nP
    if (IndCav(iP).eq.1.OR.IndCav(iP).eq.2) cycle     !Filled blocks
   
    if(surf_computation.eq.1) then
       if (IndSurf(iP).eq.3.OR.IndSurf(iP).eq.4) then
          if (DistMin(iP).le.UltraMax) then             !Ultramicro surf block
             IndCav(iP) = 4
          else if (DistMin(iP).le.MicroMax) then        !Micro surf block
             IndCav(iP) = 5
          else if (DistMin(iP).le.SmallMesoMax) then    !SmallMeso surf block
             IndCav(iP) = 6
          else if (DistMin(iP).le.LargeMesoMax) then    !LargeMeso surf block
             IndCav(iP) = 7
          else                                          !Macro surf block
             IndCav(iP) = 8
          end if                                                                
       end if

    elseif (IndCav(iP).ne.3) then

          if (DistMin(iP).le.UltraMax) then      !Ultramicro bulk block
             IndCav(iP) = 9
          else if (DistMin(iP).le.MicroMax) then       !Micro bulk block
             IndCav(iP) = 10
          else if (DistMin(iP).le.SmallMesoMax) then   !SmallMeso bulk block
             IndCav(iP) = 11
          else if (DistMin(iP).le.LargeMesoMax) then   !LargeMeso bulk block
             IndCav(iP) = 12
          else                                         !Macro bulk block   
             IndCav(iP) = 13
          end if
    end if
 end do

!  Assign the block to the suitable surface set
  do iP = 1,nP
   if(IndCav(iP).ge.4.AND.IndCav(iP).le.8) then
     MinD = INT(DistMin(iP)/DiamStep) + 1
     Surf(MinD) = Surf(MinD) + v_block
   end if
 end do

  ! compute total distribution of VMinD for volume and surface
 do iVol = 1, nVol
   MinD = iVol * DiamStep
   if(MinD.le.UltraMax) then
     UltraV = UltraV + VMinD(iVol)
     UltraS = UltraS + Surf(iVol)
   elseif(MinD.le.MicroMax) then
     MicroV = MicroV + VMinD(iVol)
     MicroS = MicroS + Surf(iVol)
   elseif(MinD.le.SmallMesoMax) then
     SmallMesoV = SmallMesoV + VMinD(iVol)
     SmallMesoS = SmallMesoS + Surf(iVol)
   elseif(MinD.le.LargeMesoMax) then
     LargeMesoV = LargeMesoV + VMinD(iVol)
     LargeMesoS = LargeMesoS + Surf(iVol)
   else
     MacroV = MacroV + VMinD(iVol)
     MacroS = MacroS + Surf(iVol)
   end if
   TotPorV = TotPorV + VMinD(iVol)
   Cumulative_VMinD(iVol) = TotPorV
 end do

 return
 end
