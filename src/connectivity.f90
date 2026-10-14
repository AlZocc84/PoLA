subroutine Connectivity(IndCon,nR,dMesh,nP)

 IMPLICIT NONE
 
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCon
 INTEGER, DIMENSION(3), INTENT(IN) :: nR
 INTEGER, DIMENSION(6) :: Connection
 INTEGER, INTENT(IN) :: nP
 REAL(DP), INTENT(IN) :: dMesh
 
 INTEGER :: iM1, iM2, XCub, YCub, ZCub, V, iP, iC, Npore, iPNew, indexP, XCubN, YCubN, ZCubN
 REAL(DP) :: XX, YY, ZZ, XP, YP, ZP

 INTERFACE
   SUBROUTINE DefinePore(iP,Npore,V,nR,IndCon)
     INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCon
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER, INTENT(IN) :: iP, Npore
     INTEGER, INTENT(OUT) :: V
   END SUBROUTINE

 END INTERFACE

 iM1 = nR(1)
 iM2 = nR(1)*nR(2)

 ! Find connectivity
 ! At the beginning the Number of pore is set to Zero
 Npore = 0

 open(1,file='pores.txt',status='unknown',form='formatted')
 write(1,'("             Pore number           Volume (A^3 / cell)")')
 do iP =1, nP
    if(IndCon(iP).ne.0) cycle    ! Only the void blocks can be part of a pore, obviously, so we skip the filled blocks
    Npore = Npore + 1
    IndCon(iP) = Npore
       
    call DefinePore(iP,Npore,V,nR,IndCon)
    write(1,'("Found the pore ",i5," with volume ",i12)') Npore, V
   
 end do
 close(1)

 ! Find the pores passing through from opposite sides of the cell
 ! Loop over the pores:
 open(2,file='connectivity.txt',status='unknown',form='formatted')
 do indexP = 1, Npore

    ! Connection is an array with 6 value:
    !              1 2 3 4 5 6
    !  Connection (_|_|_|_|_|_)
    ! 
    ! For each pre, Connection is initialized to zero.
    !  Connection= (0,0,0,0,0,0)
    !
    ! Then, for each block, if it is located at one of the cell boundaries
    ! (i.e., if the X, Y, or Z coordinate is less than dMesh), 
    ! and if the adjacent cube on the other side of the wall belongs to the same pore, then
    ! Connection(1), Connection(2), or Connection(3) is incremented, respectively.
    !
    ! Otherwise, if it is on the other side of the cell boundaries
    ! (i.e., if the X, Y, or Z coordinate is greater than (nR(i)*dMesh)+dMesh/2),
    ! and if the adjacent cube on the other side of the wall belongs to the same pore, then
    ! Connection(4), Connection(5), or Connection(6) is incremented, respectively.
    !
    ! If there are edge or vertex blocks, all corresponding components are incremented.
    ! For example:
    ! Block with X = 0.5, Y = 0.5, Z = 13.5
    ! Connection: (0,0,0,0,0,0) --> (1,1,0,0,0,0)

    Connection = 0
 
    ! loop over the blocks
    do iP = 1, nP

       if(IndCon(iP).ne.indexP) cycle

       !find the coordinates of the block
       iC = iP - 1                                                                                         
       ! Now getting the cube index coordinates
       ZCub = INT(iC/iM2) + 1 
       YCub = INT(MOD(iC,iM2)/nR(1)) + 1
       XCub = MOD(MOD(iC,iM2),nR(1)) + 1

       ! The right component/s of Connection is/are incremented
       if (XCub.eq.1) then ! this means that block iP is in the X left border
         iPNew = iM2*(ZCub-1) + iM1*(YCub-1) + nR(1) ! this is the "image" of the
         if (IndCon(iPNew).eq.indexP) Connection(1) = Connection(1) + 1
       end if
       if (XCub.eq.nR(1)) then ! this means that block iP is in the X right borde
         iPNew = iM2*(ZCub-1) + iM1*(YCub-1) + 1 ! this is the "image" of the blo
         if (IndCon(iPNew).eq.indexP) Connection(4) = Connection(4) + 1
       end if

       if (YCub.eq.1) then ! this means that block iP is in the Y left border
         iPNew = iM2*(ZCub-1) + iM1*(nR(2)-1) + XCub ! this is the "image" of the
         if (IndCon(iPNew).eq.indexP) Connection(2) = Connection(2) + 1
       end if
       if (YCub.eq.nR(2)) then ! this means that block iP is in the Y right borde
         iPNew = iM2*(ZCub-1) + 0 + XCub ! this is the "image" of the block adjac
         if (IndCon(iPNew).eq.indexP) Connection(5) = Connection(5) + 1
       end if

       if (ZCub.eq.1) then ! this means that block iP is in the Z left border
         iPNew = iM2*(nR(3)-1) + iM1*(YCub-1) + XCub ! this is the "image" of the
         if (IndCon(iPNew).eq.indexP) Connection(3) = Connection(3) + 1
       end if
       if (ZCub.eq.nR(3)) then ! this means that block iP is in the Z right borde
         iPNew = 0 + iM1*(YCub-1) + XCub ! this is the "image" of the block adjac
         if (IndCon(iPNew).eq.indexP) Connection(6) = Connection(6) + 1
       end if
    end do

    ! Now all the blocks of the pore have been analyzed
    ! If the pore passes through from opposite sides of the cell both the connection's components 1-4 (or 2-5 or 3-6) are different from zero.
    
    if ((Connection(1).ne.0.and.Connection(4).ne.0).or.&
                                    (Connection(2).ne.0.and.Connection(5).ne.0).or.&
                                          (Connection(3).ne.0.and.Connection(6).ne.0)) then

       write(2,'("Pore ",i20," REACHABLE")') indexP

    else

       write(2,'("Pore ",i20," NOT REACHABLE")') indexP
       
    end if

 end do
 close(2) 

 ! Write the xyz output file
 open(3,file='pores.xyz',status='unknown',form='formatted')
 write(3,*)nP
   write(3,*)
   iM2 = nR(1)*nR(2)
   do iP = 1, nP
     iC = iP - 1
     ZCub = INT(iC/iM2) + 1
     YCub = INT(MOD(iC,iM2)/nR(1)) + 1
     XCub = MOD(MOD(iC,iM2),nR(1)) + 1
     XX = dMesh*(XCub - 0.5)
     YY = dMesh*((YCub + 1) - 0.5)
     ZZ = dMesh*((ZCub + 1) - 0.5)

     if(IndCon(iP).eq.(-1)) then
       write(3,'("SM",3f12.6)')XX,YY,ZZ   ! SM is the name assigned to the filled blocks (the Solid Material)
     else if (IndCon(iP).eq.(-2)) then
       write(3,'("NV",3f12.6)')XX,YY,ZZ   ! NV is the name assigned to the Not accessible Volume 
     else 
       write(3,'(i5,3f12.6)')IndCon(iP),XX,YY,ZZ ! the pores are identified with numbers: 1 will be rapresented by VMD with H, 2 with He, ecc..
     end if
 
   end do
 close(3)

end
