subroutine Connectivity(IndCon,nR,dMesh,nP)

 IMPLICIT NONE
 
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCon
 INTEGER, DIMENSION(3), INTENT(IN) :: nR
 INTEGER, DIMENSION(6) :: Connection
 LOGICAL, DIMENSION(:), ALLOCATABLE :: Reachable
 INTEGER, DIMENSION(:), ALLOCATABLE :: NewInd, Updated_Index
 REAL(DP), DIMENSION(:), ALLOCATABLE ::Volumes
 INTEGER, INTENT(IN) :: nP
 REAL(DP), INTENT(IN) :: dMesh
 
 INTEGER :: iM1, iM2, XCub, YCub, ZCub, V, iP, iC, Npore, iPNew, indexP, count_pore, higher, lower, i, iPi, iPj
 REAL(DP) :: XX, YY, ZZ

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

 do iP =1, nP
    if(IndCon(iP).ne.0) cycle    ! Only the void blocks can be part of a pore, obviously, so we skip the filled blocks
    Npore = Npore + 1
    IndCon(iP) = Npore
       
    call DefinePore(iP,Npore,V,nR,IndCon)
   
 end do

 ! Allocate the Reachable array
 ! This is an array of logical values: if the 1 value is True
 ! that means that the pore n.1 is Reachable
 ! If the 2 value is False, the pore n.2 is Not Reachable

 allocate(Reachable(Npore)) 
 allocate(Updated_Index(Npore)) 
 allocate(NewInd(Npore)) 
 allocate(Volumes(Npore)) 
 Reachable = .False.
 Volumes = -1
 Updated_Index = -1
 NewInd = -1
 do i = 1, NPore
   NewInd(i) = i
 end do

 ! Find the pores passing through from opposite sides of the cell
 ! Loop over the pores:
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

       Reachable(indexP) = .True.
       !write(2,'("Pore ",i20," REACHABLE")') indexP

    else

       !write(2,'("Pore ",i20," NOT REACHABLE")') indexP
       
    end if

 end do

 ! Now, we need to adjust the number of pores:
 ! The different pores that "touch" eachothers on the side of the cell, are the same pore.
 ! Analizing the X face:
 do ZCub = 1, nR(3)
    do YCub = 1, nR(2)

        iPi = iM2*(ZCub-1) + iM1*(YCub-1) + 1
        iPj = iM2*(ZCub-1) + iM1*(YCub-1) + nR(1)
        if((IndCon(iPi).lt.0).or.(IndCon(iPJ).lt.0)) cycle
        ! Check the value of IndCav
        ! if the two blocks have the same indcon, it's ok, it's already the same pore
        ! so in the new array the index of the pore will be unchanged
        if(IndCon(iPi).ne.IndCon(iPj)) then
           lower = min(IndCon(iPi),IndCon(iPj))
           higher = max(IndCon(iPi),IndCon(iPj))
           ! we want to check if the pore with the lower indCon has already been merged with another pore
           ! in thi case the NewIndCon of the higher indCon will be IndCon(min)
           ! for example: IndCon(iPi) = 3
           ! IndCon(iPj) = 8
           ! if NewInd(3) is 3 that means that the pore 3 is not merged with another pore
           ! So we can merge pore 8 and 3 in the pore n.3, so NewInd(8) = 3
           ! istead
           ! if NewInd(3) is different from 3 that means that the pore 3 is merged with another pore (e.g 1)
           ! So also the pore 8 will be merged with the pore e.g. 1
           
           if(NewInd(lower).eq.lower) then
             NewInd(higher) = lower
           else
             NewInd(higher) = NewInd(lower)
           end if

           ! Now we need to check if there were some pores already merged with the higher, 
           ! and in that case merge it with the lower
           do i = higher, Npore
              if (NewInd(i).eq.higher) NewInd(i) = lower
           end do
           
           ! We need to update the Reachable value:
           ! If we merge the pore 3 and 8, if one of the two is accessible, then the merged pore will also be accessible
           if(Reachable(lower).or.Reachable(higher)) Reachable(lower) = .True.
        end if

    end do
 end do


 ! Analizing the Y face:
 do ZCub = 1, nR(3)
    do XCub = 1, nR(1)

        iPi = iM2*(ZCub-1) + XCub
        iPj = iM2*(ZCub-1) + iM1*(nR(2)-1) + XCub
        if((IndCon(iPi).lt.0).or.(IndCon(iPJ).lt.0)) cycle
        ! Check the value of IndCav
        ! if the two blocks have the same indcon, it's ok, it's already the same pore
        ! so in the new array the index of the pore will be unchanged
        if(IndCon(iPi).ne.IndCon(iPj)) then
           lower = min(IndCon(iPi),IndCon(iPj))
           higher = max(IndCon(iPi),IndCon(iPj))
           ! we want to check if the pore with the lower indCon has already been merged with another pore
           ! in thi case the NewIndCon of the higher indCon will be IndCon(min)
           ! for example: IndCon(iPi) = 3
           ! IndCon(iPj) = 8
           ! if NewInd(3) is 3 that means that the pore 3 is not merged with another pore
           ! So we can merge pore 8 and 3 in the pore n.3, so NewInd(8) = 3
           ! istead
           ! if NewInd(3) is different from 3 that means that the pore 3 is merged with another pore (e.g 1)
           ! So also the pore 8 will be merged with the pore e.g. 1
           
           if(NewInd(lower).eq.lower) then
             NewInd(higher) = lower
           else
             NewInd(higher) = NewInd(lower)
           end if

           ! Now we need to check if there were some pores already merged with the higher, 
           ! and in that case merge it with the lower
           do i = higher, Npore
              if (NewInd(i).eq.higher) NewInd(i) = lower
           end do
           ! We need to update the Reachable value:
           ! If we merge the pore 3 and 8, if one of the two is accessible, then the merged pore will also be accessible
           if(Reachable(lower).or.Reachable(higher)) Reachable(lower) = .True.

        end if

    end do
 end do

   
 ! Analizing the Z face:
 do YCub = 1, nR(2)
    do XCub = 1, nR(1)

        iPi = iM1*(YCub-1) + XCub
        iPj = iM2*(nR(3)-1) + iM1*(YCub-1) + XCub
        if((IndCon(iPi).lt.0).or.(IndCon(iPJ).lt.0)) cycle
        ! Check the value of IndCav
        ! if the two blocks have the same indcon, it's ok, it's already the same pore
        ! so in the new array the index of the pore will be unchanged
        if(IndCon(iPi).ne.IndCon(iPj)) then
           lower = min(IndCon(iPi),IndCon(iPj))
           higher = max(IndCon(iPi),IndCon(iPj))
           ! we want to check if the pore with the lower indCon has already been merged with another pore
           ! in thi case the NewIndCon of the higher indCon will be IndCon(min)
           ! for example: IndCon(iPi) = 3
           ! IndCon(iPj) = 8
           ! if NewInd(3) is 3 that means that the pore 3 is not merged with another pore
           ! So we can merge pore 8 and 3 in the pore n.3, so NewInd(8) = 3
           ! istead
           ! if NewInd(3) is different from 3 that means that the pore 3 is merged with another pore (e.g 1)
           ! So also the pore 8 will be merged with the pore e.g. 1
           
           if(NewInd(lower).eq.lower) then
             NewInd(higher) = lower
           else
             NewInd(higher) = NewInd(lower)
           end if

           ! Now we need to check if there were some pores already merged with the higher, 
           ! and in that case merge it with the lower
           do i = higher, Npore
              if (NewInd(i).eq.higher) NewInd(i) = lower
           end do
           ! We need to update the Reachable value:
           ! If we merge the pore 3 and 8, if one of the two is accessible, then the merged pore will also be accessible
           if(Reachable(lower).or.Reachable(higher)) Reachable(lower) = .True.

        end if

    end do
 end do


 !Now change the value of IndCon and compute the volume of the pores
 do iP = 1, nP
   if(IndCon(iP).lt.0) cycle
   !DEBUG----------------------------------
   if(NewInd(IndCon(iP)).le.0) then

     write(6,'("ERROR")')
     write(6,'("iP: ",i20)') iP
     write(6,'("IndCon(iP): ",i20)') IndCon(iP)
     write(6,'("NewInd(IndCon(iP)): ",i20)') NewInd(IndCon(iP))

   end if
   !DEBUG----------------------------------
   Volumes(NewInd(IndCon(iP))) = Volumes(NewInd(IndCon(iP))) + dMesh*dMesh*dMesh
   IndCon(iP) = NewInd(IndCon(iP)) 
 end do
   
 ! Write the xyz output file
 count_pore = 0
 open(1,file='pores.txt',status='unknown',form='formatted')
 write(1,'("             Pore number           Volume (A^3 / cell)")')
 open(2,file='connectivity.txt',status='unknown',form='formatted')
 do indexP = 1, Npore
    if (Volumes(indexP).eq.-1) cycle
    count_pore = count_pore + 1
    Updated_Index(indexP) = count_pore
    if(Reachable(indexP)) then
      write(1,'("Found the channel ",i0," with volume ",f12.6)') count_pore, Volumes(IndexP)
      write(2,'("Channel ",i0," REACHABLE")') count_pore
    else
      write(1,'("Found the pore ",i0," with volume ",f12.6)') count_pore, Volumes(IndexP)
      write(2,'("Pore ",i0," NOT REACHABLE")') count_pore
    end if
 end do

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
       if(Reachable(IndCon(iP))) then
          write(3,'("CHANNEL",i0,3f12.6)')Updated_Index(IndCon(iP)),XX,YY,ZZ ! the pores are identified with numbers: 1 will be rapresented by VMD with H, 2 with He, ecc..
       else
          write(3,'("PORE",i0,3f12.6)')Updated_Index(IndCon(iP)),XX,YY,ZZ ! the pores are identified with numbers: 1 will be rapresented by VMD with H, 2 with He, ecc..
       end if
     end if
 
   end do
 close(3)

end
