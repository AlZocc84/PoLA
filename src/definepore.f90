SUBROUTINE DefinePore(iP,Npore,V,nR,IndCon)
 
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCon
 INTEGER, DIMENSION(3), INTENT(IN) :: nR
 INTEGER, INTENT(IN) :: iP, Npore
 INTEGER, INTENT(OUT) :: V

 INTEGER:: iPNew, l, Ncurrent, Nnext, iX, iY, iZ, iC, Nel, i, iM2, iCN
 INTEGER:: XCub, YCub, ZCub, XCubN, YCubN, ZCubN
 INTEGER, DIMENSION(:), ALLOCATABLE :: next
 INTEGER, DIMENSION(:), ALLOCATABLE :: current
 
 INTERFACE

  FUNCTION MoveCub(iP,iX,iY,iZ,nR)
   INTEGER, INTENT(IN) :: iP,iX,iY,iZ
   INTEGER, DIMENSION(3), INTENT(IN) :: nR
   INTEGER :: MoveCub
  END FUNCTION MoveCub
 
 END INTERFACE

 iM2 = nR(1)*nR(2)

 allocate(current(6))
 current = 0
 Ncurrent = 0
 V=1
! we start moving around the 26 blocks near the chosen one (iP)
 l = 1  
 do iX = -l,l
    do iY = -l,l 
       do iZ = -l,l
          if (iX == 0 .and. iY == 0 .and. iZ == 0) cycle
          if(.not.((iX == 0 .and. iY == 0 ).or.(iX == 0 .and. iZ == 0).or.(iY == 0 .and. iZ == 0))) cycle

          ! Find the "coordinates" of the block on the mesh
          iC = iP - 1
          ZCub = INT(iC/iM2) + 1
          YCub = INT(MOD(iC,iM2)/nR(1)) + 1
          XCub = MOD(MOD(iC,iM2),nR(1)) + 1

          XCubN = XCub + iX
          if(XCubN.lt.1) cycle
          if(XCubN.gt.nR(1)) cycle

          YCubN = YCub + iY                                                            
          if(YCubN.lt.1) cycle
          if(YCubN.gt.nR(2)) cycle

          ZCubN = ZCub + iZ
          if(ZCubN.lt.1) cycle
          if(ZCubN.gt.nR(3)) cycle

          iPNew = iM2*(ZCubN-1) + nR(1)*(YCubN-1) + XCubN

          if(IndCon(iPNew).ne.0) cycle
          IndCon(iPNew) = Npore
          V=V+1
          NCurrent = NCurrent +1
          current(Ncurrent) = iPNew
       end do
    end do
 end do

 Nel = Ncurrent

 do while(Nel.gt.0)
    
    allocate(next(Nel*6))
    next = 0
    Nnext = 0

    do i = 1, Nel
        iCN = current(i)
        l = 1  
        do iX = -l,l
           do iY = -l,l 
              do iZ = -l,l
                 if (iX == 0 .and. iY == 0 .and. iZ == 0) cycle
                 if(.not.((iX == 0 .and. iY == 0 ).or.(iX == 0 .and. iZ == 0).or.(iY == 0 .and. iZ == 0))) cycle

                 ! Find the "coordinates" of the block on the mesh
                 iC = iCN - 1
                 ZCub = INT(iC/iM2) + 1
                 YCub = INT(MOD(iC,iM2)/nR(1)) + 1
                 XCub = MOD(MOD(iC,iM2),nR(1)) + 1
                 
                 XCubN = XCub + iX
                 if(XCubN.lt.1) cycle
                 if(XCubN.gt.nR(1)) cycle
                 
                 YCubN = YCub + iY                                                            
                 if(YCubN.lt.1) cycle
                 if(YCubN.gt.nR(2)) cycle
                 
                 ZCubN = ZCub + iZ
                 if(ZCubN.lt.1) cycle
                 if(ZCubN.gt.nR(3)) cycle
                 
                 iPNew = iM2*(ZCubN-1) + nR(1)*(YCubN-1) + XCubN


                 if(IndCon(iPNew).ne.0) cycle
                 IndCon(iPNew) = Npore
                 V=V+1
                 Nnext = Nnext +1
                 next(Nnext) = iPNew
              end do
           end do
        end do

    end do 
    
    Nel = Nnext

    deallocate(current)
    allocate(current(Nnext))

    current(1:Nnext) = next(1:Nnext)

    deallocate(next)

 end do
END 
