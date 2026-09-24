subroutine Connectivity

 IMPLICIT NONE
 
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

 INTERFACE
   SUBROUTINE DefinePore(iP,Npore,V,nR,IndCon)
     INTEGER, DIMENSION(:), INTENT(INOUT) :: IndCon
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER, INTENT(IN) :: iP, Npore
     INTEGER, INTENT(OUT) :: V
   END SUBROUTINE

 END INTERFACE

 ! Find connectivity
 Npore = 0
!open(7,file='indcon_pre.xyz',status='unknown',form='formatted')                                                            
! write(7,*)nP
!   write(7,*)
!   do iP = 1, nP
!     write(7,'(i13)')IndCon(iP)
!   end do
! close(7)

if (Connect.gt.0) then
    open(3,file='connectivity.txt',status='unknown',form='formatted')
    do iP =1, nP
       if(IndCon(iP).ne.0) cycle
       !write(3,'("IndCon(iP) at beginnig ",i5)') IndCon(iP)
       Npore = Npore + 1
       IndCon(iP) = Npore
       
       !debug 
   !    write(3,'("iP= ",i13)') iP
   !    write(3,'("IndCon after change",i5)') IndCon(iP)
   !    write(3,'("Npore ",i5)') Npore
       
       call DefinePore(iP,Npore,V,nR,IndCon)
       write(3,'("Found the pore ",i5," with volume ",i12)') Npore, V
   
       !debug 
   !    write(3,'("IndCon after connect",i5)') IndCon(iP)
   !    write(3,'("Npore ",i5)') Npore
       !write(3,'("Found the pore ",i5," with volume ",i12)') IndCon(iP), V
   
    end do
    close(3)
   
   !open(7,file='indcon_post_conn.xyz',status='unknown',form='formatted')
   ! write(7,*)nP
   !   write(7,*)
   !   do iP = 1, nP
   !     write(7,'(i13)')IndCon(iP)
   !   end do
   ! close(7)
   
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

