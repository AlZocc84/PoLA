 SUBROUTINE Texture(nP,dMesh,nVol,DiamStep,VMinD,Cumulative_VMinD,UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV, &     
                   surf_computation,DistMin,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,Rad,IndCav,Surf,Accessible,IndCon)

 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
 INTEGER, INTENT(IN) :: nVol, surf_computation, nP
 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCav, IndCon
 REAL(DP), DIMENSION(:), INTENT(INOUT) :: VMinD
 REAL(DP), INTENT(IN) :: DiamStep, dMesh, Rad
 REAL(DP), INTENT(OUT) :: UltraV, MicroV, SmallMesoV, LargeMesoV, MacroV, TotPorV
 REAL(DP), INTENT(OUT) :: UltraS, MicroS, SmallMesoS, LargeMesoS, MacroS
 REAL(DP), DIMENSION(:), INTENT(IN) :: DistMin
 REAL(DP), DIMENSION(:), INTENT(OUT) :: Cumulative_VMinD
 REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: Surf
 LOGICAL, INTENT(IN) :: Accessible
 INTEGER :: iVol, iP, MinD
 REAL(DP) :: NAV, v_block
 REAL(DP), PARAMETER :: Zero=0.0d0, Two=2.0d0
 REAL(DP), PARAMETER :: UltraMax=7.0d0, MicroMax=20.0d0, SmallMesoMax=35.0d0, LargeMesoMax=50.0d0


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
 
 NAV = Zero

 v_block = dMesh*dMesh*dMesh
 
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
 if (IndCav(iP).eq.1.OR.IndCav(iP).eq.2) then     !Filled blocks
     cycle

   else if (IndCav(iP).eq.0) then
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

   else if (IndCav(iP).eq.3.OR.IndCav(iP).eq.4) then
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
 end do

 if(Accessible) then
    do iP=1,nP
       if((IndCav(iP).eq.9).AND.(DistMin(iP).le.(Rad*2.0))) then !this works only if 2Rad le UltraMax TODO:generalize 
         IndCav(iP) = 3
         IndCon(iP) = -2
         NAV = NAV + v_block
       end if
    end do
 end if
    
!  Assign the block to the suitable pore set
! At this point is only the accessible volume
  do iP = 1,nP
   if(IndCav(iP).eq.3) then
     MinD = INT(DistMin(iP)/DiamStep) + 1 
     VMinD(MinD) = VMinD(MinD) - v_block
   else if(IndCav(iP).ge.4.AND.IndCav(iP).le.8) then
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
