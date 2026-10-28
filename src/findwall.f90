 SUBROUTINE Find_Wall(iP,iPN,nR,dMesh,XCub,YCub,ZCub,dX,dY,dZ,IndCav,Found,R)

 IMPLICIT NONE

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
 INTEGER, INTENT(IN) :: iP, XCub, YCub, ZCub
 INTEGER, DIMENSION(3), INTENT(IN) :: nR 
 INTEGER, DIMENSION(:), INTENT(IN) :: IndCav
 INTEGER, INTENT(OUT) :: iPN
 INTEGER :: iM1, iM2, iStep, dXcub, dYcub, dZcub, MaxR
 INTEGER ::  XCubN, YCubN, ZCubN, XCubTest, YCubTest, ZCubTest
 REAL(DP), INTENT(IN) :: dMesh, dX, dY, dZ
 REAL(DP) :: X, Y, Z, R2, R, DeltaX, DeltaY, DeltaZ
 LOGICAL, INTENT(OUT) :: Found
 LOGICAL Void

 iM1 = nR(1)  
 iM2 = nR(1)*nR(2)  
 MaxR= max(nR(1),nR(2),nR(3))

! Move along the selected direction from iP until a filled block is found
 iStep = 0
 iPN = 0
 Void = .True.
 Found = .True.
 do while (Void)
   iStep = iStep + 1

   ! Now getting the actual length of vector. steps are effecetively dMesh long (or we risk missing filled blocks.)
   ! If getting out of the cell on any direction mark not found.
   X = iStep*dX*dMesh
   Y = iStep*dY*dMesh
   Z = iStep*dZ*dMesh

   if(abs(X).ge.nR(1)*dMesh .OR. abs(Y).ge.nR(2)*dMesh .OR. abs(Z).ge.nR(3)*dMesh) then  
    Found = .False.
    exit
   end if 
  
! How many blocks we have to cross in X, Y, Z respectively?
   dXcub = NINT(X/dMesh) 
   dYcub = NINT(Y/dMesh) 
   dZcub = NINT(Z/dMesh) 

! The mesh coord. of the block: XCubN is the real index (that can be outside the cell, also) used
! to compute the distance from the wall; XCubeTest is the *cell* index (of which XCubN is possibly the periodic
! image) used to verify if the block is full or empty

   XCubN = XCub + dXcub
   XCubTest = XCubN
   if(XCubTest.lt.1) XCubTest = Mod(XCubTest,nR(1)) + nR(1)
   if(XCubTest.gt.nR(1)) then 
      XCubTest = Mod(XCubTest,nR(1))
      if (XCubTest.eq.0) XCubTest = nR(1)
   endif 

   YCubN = YCub + dYcub
   YCubTest = YCubN
   if(YCubTest.lt.1) YCubTest = Mod(YCubTest,nR(2)) + nR(2)
   if(YCubTest.gt.nR(2)) then
      YCubTest = Mod(YCubTest,nR(2)) 
      if (YCubTest.eq.0) YCubTest = nR(2)
   endif

   ZCubN = ZCub + dZcub
   ZCubTest = ZCubN
   if(ZCubTest.lt.1) ZCubTest = Mod(ZCubTest,nR(3)) + nR(3)
   if(ZCubTest.gt.nR(3)) then
      ZCubTest = Mod(ZCubTest,nR(3))
      if (ZCubTest.eq.0) ZCubTest = nR(3)
   endif

   iPN = iM2*(ZCubTest-1) + iM1*(YCubTest-1) + XCubTest 

   if(IndCav(iPN).eq.1) Void = .False.
   if(iPN.eq.iP) then
     Found = .False.
     exit 
   end if
 end do

! When the wall has been hit, compute the distance from point iP
 if (Found) then
!  Center_to_center distance
   DeltaX = (XCubN - XCub) * dMesh
   DeltaY = (YCubN - YCub) * dMesh
   DeltaZ = (ZCubN - ZCub) * dMesh
   R2 = DeltaX**2 + DeltaY**2 + DeltaZ**2
   R = sqrt(R2)
 end if
 
 return
 end

