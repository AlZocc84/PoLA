 SUBROUTINE Texture(nP,dMesh,nVol,DiamStep,VMinD,Cumulative_VMinD,UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV,nR,iM1,iM2, &     
                   surf_computation,DistMin,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,RadAV,IndCav,Surf,Accessible,IndSurf,IndCon,Connect)

 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
 INTEGER, INTENT(IN) :: nVol, surf_computation, nP, iM1, iM2, Connect
 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCav, IndCon
 INTEGER, DIMENSION(:), INTENT(IN) :: IndSurf, nR
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: VMinD
 REAL(DP), INTENT(IN) :: DiamStep, dMesh, RadAV
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
 LOGICAL, DIMENSION(:), ALLOCATABLE :: NAccVol
 INTERFACE
   FUNCTION MoveCub(iP,iX,iY,iZ,nR)
     INTEGER, INTENT(IN) :: iP,iX,iY,iZ
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER :: MoveCub
   END FUNCTION MoveCub
 END INTERFACE

 NAV = Zero
 v_block = dMesh*dMesh*dMesh
 Rad1 = RadAV+(dMesh*0.1)

! Allocate the array for the NOT Accessible Volume. 
! NAccVal= False : ISN'T NOT Accessible Volume 
! (Can be a filled block, or Accessible)

 allocate(NAccVol(nP))
 NAccVol = .False.

!Compute the Accessible volume, if required
 if(Accessible) then
   do iP=1,nP
     if(IndCav(iP).ne.0) cycle

     ! Find the block Cartesian coordinates
     iC = iP - 1
     ZP = INT(iC/(iM2))*dMesh + dMesh/2.0d0
     YP = INT(MOD(iC,iM2)/nR(1))*dMesh + dMesh/2.0d0
     XP = MOD(MOD(iC,iM2),nR(1))*dMesh + dMesh/2.0d0
     ! Check all the blocks that are closer than Rad1. 
     Overlap = .False.  !If TRUE this probe overlaps to the wall (then it will be considered Not Accessible)

     lCube = nint(Rad1 / dMesh) + 1
     do iX = -lCube,lCube
       if(Overlap) exit                                                                                                                               
       do iY = -lCube,lCube
         if(Overlap) exit
         do iZ = -lCube,lCube
           if(Overlap) exit
           iPNew = MoveCub(iP, iX, iY, iZ, nR)
           ! Skip if the new block is void (IndCav=0) or has already been classified as Not Accessible (IndCav=3)
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
         
           ! If D is lower than Rad1 the block overlaps with the probe, and it is part of the Not Accessible volume
           if (D.le.Rad1) then
              Overlap = .True.
              IndCav(iP) = 3
              NAccVol(iP) = .True. 
           end if

         end do
       end do
     end do

   end do
 end if

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
    
    if (IndSurf(iP).eq.1.OR.IndSurf(iP).eq.2) then   !Are part of the surface
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

!  Fill VMinD array with the accessible volume  
 do iP = 1,nP
   if(IndCav(iP).eq.0) then     ! At this point all blocks should have IndCav different from 0, otherwise there is an error and the program is stopped
     write(6,'("IndCav of block ",i20," is 0")') iP
     write(6,'("At this point this should not happen.")')
     write(6,'("Program stopped.")')
     stop
   else if (IndCav(iP).eq.1.OR.IndCav(iP).eq.2) then ! Skip the filled blocks (1,2) 
     cycle
   end if
    
   if(.not.(NAccVol(iP))) then           ! VMinD is filled only with Accessible volume
     MinD = DistMin(iP) 
     iVol = INT(MinD/DiamStep) + 1
     VMinD(iVol) = VMinD(iVol) + v_block 
   end if
 end do

 ! Compute the volume of each classification (UltraMicro, Micro etc.)
 ! The UltraMicro volume is given by the ultramicro bulk (IndCav=9) and by the ACCESSIBLE volume of the ultamicro surface (IndCav=4 and not NaccVol)
 ! (and the same for all the other classifications)
 do iP= 1, nP
   if((IndCav(iP).eq.9).OR.(IndCav(iP).eq.4.AND.(.not.NaccVol(iP)))) then
     UltraV = UltraV + v_block    
   elseif((IndCav(iP).eq.10).OR.(IndCav(iP).eq.5.AND.(.not.NaccVol(iP)))) then
     MicroV = MicroV + v_block    
   elseif((IndCav(iP).eq.11).OR.(IndCav(iP).eq.6.AND.(.not.NaccVol(iP)))) then
     SmallMesoV = SmallMesoV + v_block    
   elseif((IndCav(iP).eq.12).OR.(IndCav(iP).eq.7.AND.(.not.NaccVol(iP)))) then
     LargeMesoV = LargeMesoV + v_block    
   elseif((IndCav(iP).eq.13).OR.(IndCav(iP).eq.8.AND.(.not.NaccVol(iP)))) then
     MacroV = MacroV + v_block    
   end if
 end do

 do iP= 1, nP
 ! Compute the surface volume of each classification (UltraMicro, Micro etc.)
 !  And Surf(iVol) is filled
   if(IndCav(iP).eq.4) then
     UltraS = UltraS + v_block    
     Surf(iVol) = Surf(iVol) + v_block
   elseif(IndCav(iP).eq.5) then
     MicroS = MicroS + v_block    
     Surf(iVol) = Surf(iVol) + v_block
   elseif(IndCav(iP).eq.6) then
     SmallMesoS = SmallMesoS + v_block    
     Surf(iVol) = Surf(iVol) + v_block
   elseif(IndCav(iP).eq.7) then
     LargeMesoS = LargeMesoS + v_block    
     Surf(iVol) = Surf(iVol) + v_block
   elseif(IndCav(iP).eq.8) then
     MacroS = MacroS + v_block    
     Surf(iVol) = Surf(iVol) + v_block
   end if
 end do

  ! compute the total volume and the Cumulative VMinD
 do iVol = 1, nVol
   TotPorV = TotPorV + VMinD(iVol)
   Cumulative_VMinD(iVol) = TotPorV
 end do

 return
 end
