program pore_local_analysis

!****************************************************************************************************************
! Analyze the Local Porous Volume Distribution of a given porous material, through detailed geometrical analysis
!****************************************************************************************************************

 IMPLICIT NONE

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
 INTEGER, DIMENSION(:), ALLOCATABLE :: IndCav, SurfCenter
 INTEGER :: i, iP, nP, nAtm, iPN, iPNopp, nMono_dens, n_angles, MinD
 INTEGER ::  iM1, iM2, XCub, YCub, ZCub, iC, Print_xyz, nVol, nS, surf_computation
 INTEGER, DIMENSION(3) :: nR
 REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0, Six=6.0d0
 REAL(DP), PARAMETER :: UltraMax=7.0d0, MicroMax=20.0d0, SmallMesoMax=35.0d0, LargeMesoMax=50.0d0
 REAL(DP), PARAMETER :: A3_gMol_to_cm3_g = 0.6022d0
 REAL(DP), DIMENSION(:), ALLOCATABLE :: Cumulative_VMinD, VMinD, DistMin, Surf
 REAL(DP), DIMENSION(:,:), ALLOCATABLE:: versors
 REAL(DP) :: UltraV, MicroV, SmallMesoV, LargeMesoV, MacroV, TotPorV 
 REAL(DP) :: UltraS, MicroS, SmallMesoS, LargeMesoS, MacroS
 REAL(DP) :: dX, dY, dZ, R, RMin, RMin_temp, RMax, RN, Ropp
 REAL(DP) :: dMesh, v_block, DiamStep, Closed_Thresh, MaxDiameter
 REAL(DP) :: SkV, SkM, SkD, Rad

 CHARACTER(2), DIMENSION(:), ALLOCATABLE :: AtmSym

 LOGICAL:: Found, Accessible
 REAL(DP) :: Fact

 INTERFACE
   SUBROUTINE Input(dMesh,nR,nP,nAtm,AtmSym,IndCav,DiamStep,RMin,MaxDiameter,Closed_Thresh,Print_xyz, &
                 surf_computation,Rad,n_angles,versors,Accessible)
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: IndCav
     INTEGER, INTENT(OUT) :: nP, nAtm, Print_xyz, surf_computation, n_angles
     INTEGER, DIMENSION(3), INTENT(OUT) :: nR
     REAL(DP), INTENT(OUT) :: dMesh, DiamStep, RMin, Closed_Thresh, MaxDiameter, Rad
     REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: versors
     CHARACTER(2), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: AtmSym
     LOGICAL, INTENT(OUT) :: Accessible
   END SUBROUTINE Input
   SUBROUTINE Find_Wall(iP,iPN,nR,dMesh,XCub,YCub,ZCub,dX,dY,dZ,IndCav,Found,R)
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     INTEGER, INTENT(IN) :: iP, XCub, YCub, ZCub
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER, DIMENSION(:), INTENT(IN) :: IndCav
     INTEGER, INTENT(OUT) :: iPN
     REAL(DP), INTENT(IN) :: dMesh, dX, dY, dZ
     LOGICAL, INTENT(OUT) :: Found
     REAL(DP) :: R
   END SUBROUTINE Find_Wall
   SUBROUTINE Skel_Dens(dMesh,nP,IndCav,nAtm,AtmSym,SkV,SkM,SkD)
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     INTEGER, INTENT(IN) :: nP, nAtm
     INTEGER, DIMENSION(:), INTENT(IN) :: IndCav
     REAL(DP), INTENT(IN) :: dMesh
     REAL(DP), INTENT(OUT) :: SkV, SkM, SkD
     CHARACTER(2), DIMENSION(:), INTENT(IN) :: AtmSym
   END SUBROUTINE Skel_Dens
     SUBROUTINE Surface(Surf_Computation,dMesh,nR,nP,IndCav,Rad,SurfCenter,nS,nMono_dens)
     IMPLICIT NONE

     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

     INTEGER, INTENT(IN) :: Surf_Computation, nP                                                                                                      
     INTEGER, INTENT(OUT) :: nS, nMono_dens
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCav
     INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: SurfCenter
     REAL(DP), INTENT(IN) :: dMesh, Rad

   END SUBROUTINE Surface                                             
   SUBROUTINE Texture(nP,dMesh,nVol,DiamStep,VMinD,Cumulative_VMinD,UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV, &                                   
                   surf_computation,DistMin,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,Rad,IndCav,Surf,Accessible)
     IMPLICIT NONE
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     INTEGER, INTENT(IN) :: nVol, surf_computation, nP
     INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCav
     REAL(DP), DIMENSION(:), INTENT(INOUT) :: VMinD
     REAL(DP), INTENT(IN) :: DiamStep, dMesh, Rad
     REAL(DP), INTENT(OUT) :: UltraV, MicroV, SmallMesoV, LargeMesoV, MacroV, TotPorV
     REAL(DP), INTENT(OUT) :: UltraS, MicroS, SmallMesoS, LargeMesoS, MacroS
     REAL(DP), DIMENSION(:), INTENT(IN) :: DistMin
     REAL(DP), DIMENSION(:), INTENT(OUT) :: Cumulative_VMinD
     REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: Surf
     LOGICAL, INTENT(IN) :: Accessible
   END SUBROUTINE Texture
   SUBROUTINE Output(dMesh,nR,nP,nVol,IndCav,DiamStep,Cumulative_VMinD,VMinD,Accessible, &
             UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV,SkV,SkM,SkD,Print_xyz, &
             Rad,surf_computation,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,nMono_dens)
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     INTEGER, DIMENSION(:), INTENT(IN) :: IndCav
     INTEGER, INTENT(IN) :: nP, nVol, Print_xyz, surf_computation, nMono_dens
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     REAL(DP), INTENT(IN) :: dMesh, DiamStep, SkV, SkM, SkD
     REAL(DP), INTENT(IN) :: UltraV, MicroV, SmallMesoV, LargeMesoV, MacroV, TotPorV
     REAL(DP), INTENT(IN) :: UltraS, MicroS, SmallMesoS, LargeMesoS, MacroS, Rad
     REAL(DP), DIMENSION(:), INTENT(IN) :: Cumulative_VMinD, VMinD
     LOGICAL, INTENT(IN) :: Accessible
   END SUBROUTINE Output
 END INTERFACE

! Read in the mesh coordinates of the material: input, xyz file; output, 
! IndCav vector starts with 1 for "occupied" and 0 for void blocks
! and will contain a detailed classification of the void blocks (Volume and surface: <7A(ultramicro), >7 and <20 (micro), meso etc.)

 call Input(dMesh,nR,nP,nAtm,AtmSym,IndCav,DiamStep,RMin,MaxDiameter,Closed_Thresh,Print_xyz, &
                 surf_computation,Rad,n_angles,versors,Accessible)

 iM1 = nR(1)
 iM2 = nR(1)*nR(2)

 allocate(DistMin(nP))
 DistMin = 0.0d0

!nVol = INT(MaxDiameter / DiamStep)
! To obtain comparable files of local volumes and VMinD, the maximum diameter in the output
! is set to 120 A, way larger than the affordable porosities...
 nVol = INT(120/DiamStep)
 allocate(Cumulative_VMinD(nVol),VMinD(nVol))
 VMinD = Zero
 Cumulative_vMinD = Zero
 v_block = dMesh*dMesh*dMesh

! Loop on the void points
 open(1,file='distance.txt',status='unknown',form='formatted')
 write(1,'("# Index    XCub  YCub  ZCub   MinD")')
 do iP = 1, nP

   if(IndCav(iP).eq.1) then
     iC = iP - 1
     ZCub = INT(iC/iM2) + 1
     YCub = INT(MOD(iC,iM2)/nR(1)) + 1
     XCub = MOD(MOD(iC,iM2),nR(1)) + 1
     write(1,'(" ",i7,"   ",i3,"   ",i3,"   ",i3,"   -1")') iP, XCub, YCub, ZCub

     cycle
   end if
! RMin_temp/RMax will be the shortest/longest distances between opposite walls for this point
   RMin_temp = RMin
   RMax = Zero
! Find the "coordinates" of the block on the mesh
   iC = iP - 1
   ZCub = INT(iC/iM2) + 1
   YCub = INT(MOD(iC,iM2)/iM1) + 1
   XCub = MOD(MOD(iC,iM2),iM1) + 1

! Scan the spherical angles around iP
   do i = 1, n_angles
! Find the unit vector 
       dX = versors(1,i)
       dY = versors(2,i)
       dZ = versors(3,i)

! Find the wall in this direction
       Call Find_Wall(iP,iPN,nR,dMesh,XCub,YCub,ZCub,dX,dY,dZ,IndCav,Found,RN)
       if(.not.Found) cycle
! Search the same line in the opposite direction
       dX = -dX       
       dY = -dY       
       dZ = -dZ       
! Find the wall in this direction also
       Call Find_Wall(iP,iPNopp,nR,dMesh,XCub,YCub,ZCub,dX,dY,dZ,IndCav,Found,Ropp)
       if(.not.Found) cycle
! Find the distance between points iPN and iPNopp
       R = RN + Ropp
! Check the minimum/maximum distances
       if(R.le.RMin_temp) RMin_temp = R
       if(R.ge.RMax) RMax = R
   end do  

! Save each block and its minimum distance
   DistMin(iP) = RMin_temp
   iC = iP - 1
   ZCub = INT(iC/iM2) + 1
   YCub = INT(MOD(iC,iM2)/nR(1)) + 1
   XCub = MOD(MOD(iC,iM2),nR(1)) + 1
   write(1,'(" ",i7,"   ",i3,"   ",i3,"   ",i3,"   ",f8.4)') iP, XCub, YCub, ZCub, DistMin(iP)

! If RMax falls below a given threshold this is considered a "closed" pore, inaccessible
! by adsorbates, then this block becomes "filled" and is not considered in the porous volume analysis
   if(RMax.gt.Zero.AND.RMax.lt.Closed_Thresh) then
     IndCav(iP) = 2
     cycle
   end if
 end do
 close(1)
! Compute skeletal density
 call Skel_Dens(dMesh,nP,IndCav,nAtm,AtmSym,SkV,SkM,SkD)
 
! Assign the block to the suitable pore set
! Volumes associated to RMin are not added to VMinD
! At this point is the total volume (accessible and not accessible)
 do iP=1,nP
   if(IndCav(iP).eq.0) then
     MinD = INT(DistMin(iP)/DiamStep) + 1
     VMinD(MinD) = VMinD(MinD) + v_block 
   end if
 end do

 ! Compute specific surface
 call Surface(Surf_Computation,dMesh,nR,nP,IndCav,Rad,SurfCenter,nS,nMono_dens) 
 
! Compute porous volumes
 call Texture(nP,dMesh,nVol,DiamStep,VMinD,Cumulative_VMinD,UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV, &     
                   surf_computation,DistMin,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,Rad,IndCav,Surf,Accessible)

! Fact = A3_gMol_to_cm3_g / SkM
! open(5,file='Srf.txt', status='unknown', form='formatted')
! do MinD=1, nVol
!   write(5,'(i12,f14.1,e20.8)') MinD, Surf(MinD), Fact*Surf(MinD)
! end do
! close(5)

! Write the results
 call Output(dMesh,nR,nP,nVol,IndCav,DiamStep,Cumulative_VMinD,VMinD,Accessible, &
             UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV,SkV,SkM,SkD,Print_xyz, &
             Rad,surf_computation,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,nMono_dens)

 stop
 end

