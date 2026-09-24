 SUBROUTINE Surface(Surf_Computation,dMesh,nR,nP,Rad,SurfCenter,nS,nMono_dens,IndSurf)

 IMPLICIT NONE

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

 INTEGER, INTENT(IN) :: Surf_Computation, nP 
 INTEGER, INTENT(OUT) :: nS, nMono_dens
 INTEGER, DIMENSION(3), INTENT(IN) :: nR
 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndSurf
 INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: SurfCenter
 REAL(DP), INTENT(IN) :: dMesh, Rad

 REAL(DP) :: XP, YP, ZP, XP1, YP1, ZP1, Rad1
 REAL(DP) :: Dist_X, Dist_Y, Dist_Z, D, MonoVol
 REAL(DP), PARAMETER :: N2_Density = 0.808d0

INTEGER :: iP, iPNew, iC, iX, iY, iZ, lCube, iM2
Logical :: Overlap, Touch

 INTERFACE
   FUNCTION MoveCub(iP,iX,iY,iZ,nR)
     INTEGER, INTENT(IN) :: iP,iX,iY,iZ
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER :: MoveCub
   END FUNCTION MoveCub
 END INTERFACE

! Apply the spherical probe model to identify the occupiable void blocks (i.e. void blocks
! close to the surface and suitable to host the center of a spherical probe touching the surface, but
! not overlapping it) and to compute the specific surface area from a "complete monolayer" of spherical
! probes.

! Rad is the radius of the probe "molecule"; Rad1 is Rad plus dMesh/2.
! The probe overlaps to blocks falling inside Rad, and "touches" blocks falling between Rad and Rad1.
! Blocks touched by the probe contribute to the porous surface.

!SurfCenter is a VcT with all the boxes that are the center of the spherical probes
!so at the beginning the first element is -1
 allocate(SurfCenter(nP))
 SurfCenter(1) = -1

 nS = 0       ! Number of blocks that are the center of the spherical probes

 iM2 = nR(1) * nR(2)
 Rad1 = Rad + 0.6*dMesh

 do iP = 1, nP
! Loop only on void blocks
   if(IndSurf(iP).ne.0) cycle
! Find the block Cartesian coordinates
   iC = iP - 1
   ZP = INT(iC/(iM2))*dMesh + dMesh/2.0d0
   YP = INT(MOD(iC,iM2)/nR(1))*dMesh + dMesh/2.0d0
   XP = MOD(MOD(iC,iM2),nR(1))*dMesh + dMesh/2.0d0
! Check all the blocks that could fall inside Rad or Rad1. 
   Overlap = .False.  !If TRUE this probe overlaps to the wall (then it will be discarded)
   Touch = .False.  !If TRUE this probe is at right distance to the wall (then it will be accepted)

   lCube = nint(Rad1 / dMesh) + 1
   do iX = -lCube,lCube
     if(Overlap) exit
     do iY = -lCube,lCube
       if(Overlap) exit
       do iZ = -lCube,lCube
         if(Overlap) exit
         iPNew = MoveCub(iP, iX, iY, iZ, nR)
! Skip if the new block is void (IndSurf=0) or has already been classified as occupiable (IndSurf=3)
         if (IndSurf(iPNew).eq.0.or.IndSurf(iPNew).eq.3) cycle
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
         if (D.lt.Rad) Overlap = .True.
! If D is between Rad and Rad1 the probe touches this block 
         if (.not.Overlap .AND. D.le.Rad1) Touch = .True.
       end do
     end do
   end do

! If this probe touched some blocks, and it was not overlapped, change the IndSurf of the surface blocks
   if (Touch .AND. .not.Overlap) then
     nS = nS + 1
     IndSurf(iP) = 3
     SurfCenter(nS) = iP
     SurfCenter(nS+1) = -1
   end if
 end do 

! Compute the number of void blocks belonging to one of the spherical probes which form the surface monolayer
 do iP = 1, nP
! Loop only on blocks which could host the center of a probe touching the surface (occupiable blocks)
   if(IndSurf(iP).ne.3) cycle
! Find the block Cartesian coordinates
   iC = iP - 1
   ZP = INT(iC/(iM2))*dMesh + dMesh/2.0d0
   YP = INT(MOD(iC,iM2)/nR(1))*dMesh + dMesh/2.0d0
   XP = MOD(MOD(iC,iM2),nR(1))*dMesh + dMesh/2.0d0
! Check all the blocks that could fall inside Rad 
!   lCube = nint(Rad / dMesh) + 1
   lCube = nint(Rad / dMesh) 
   do iX = -lCube,lCube
     do iY = -lCube,lCube
       do iZ = -lCube,lCube
         iPNew = MoveCub(iP, iX, iY, iZ, nR)
! Skip if the block is filled or has been already assigned
         if (IndSurf(iPNew).ne.0) cycle
         ! IndSurf(iPNew) = 4                     uncomment this to get surface volume with cubical probe
 ! Find the coordinates of the new block                                                         comment from here
          iC = iPNew - 1                                                                          ! to get surface volume 
          ZP1 = INT(iC/(iM2))*dMesh + dMesh/2.0d0                                                  ! with cubical probe
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
 ! If the distance is lower than Rad, so that this block would fall inside one of the spherical probes
 ! forming the surface monolayer, its IndSurf is _temporarily_ set to 4
          if (D.le.Rad) IndSurf(iPNew) = 4                                                          ! comment to here
       end do
     end do
   end do
 end do

! Compute the volume of the surface monolayer
 MonoVol = 0.0d0
 do iP = 1, nP
! If the block belongs to the monolayer (IndSurf=3 is the center of a probe, 
! IndSurf=4 is falling inside some probe) increment the monolayer volume
   if (IndSurf(iP).eq.3 .OR. IndSurf(iP).eq.4) MonoVol = MonoVol + dMesh*dMesh*dMesh
 end do

! Find the number of spherical probes forming the monolayer
! This if we use the liquid N2 density to estimate the molecules in the monolayer
 nMono_dens = INT(N2_density * 0.6022d0 * MonoVol / 28.0d0)

 end
