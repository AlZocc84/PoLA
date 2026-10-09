subroutine Connectivity(IndCon,nR,dMesh,nP)

 IMPLICIT NONE
 
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

 INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCon
 INTEGER, DIMENSION(3), INTENT(IN) :: nR
 INTEGER, INTENT(IN) :: nP
 REAL(DP), INTENT(IN) :: dMesh
 
 INTEGER :: iM2, XCub, YCub, ZCub, V, iP, iC, Npore
 REAL(DP) :: XX, YY, ZZ

 INTERFACE
   SUBROUTINE DefinePore(iP,Npore,V,nR,IndCon)
     INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCon
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER, INTENT(IN) :: iP, Npore
     INTEGER, INTENT(OUT) :: V
   END SUBROUTINE

 END INTERFACE

 ! Find connectivity
 ! At the beginning the Number of pore is set to Zero
 Npore = 0

 open(3,file='connectivity.txt',status='unknown',form='formatted')
 do iP =1, nP
    if(IndCon(iP).ne.0) cycle    ! Only the void blocks can be part of a pore, obviously, so we skip the filled blocks
    Npore = Npore + 1
    IndCon(iP) = Npore
       
    call DefinePore(iP,Npore,V,nR,IndCon)
    write(3,'("Found the pore ",i5," with volume ",i12)') Npore, V
   
 end do
 close(3)
   
 ! Write the output file
 open(4,file='connectivity.xyz',status='unknown',form='formatted')
 write(4,*)nP
   write(4,*)
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
       write(4,'("XX",3f12.6)')XX,YY,ZZ
     else if (IndCon(iP).eq.(-2)) then
       write(4,'("YY",3f12.6)')XX,YY,ZZ
     else 
       write(4,'(i5,3f12.6)')IndCon(iP),XX,YY,ZZ
     end if
 
   end do
 close(4)

end
