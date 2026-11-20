subroutine checkface(iPi,iPj,nP,NPore,IndCon,NewInd,Reachable)                                     

 IMPLICIT NONE

 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)

 INTEGER, DIMENSION(:), INTENT(INOUT) :: NewInd
 LOGICAL, DIMENSION(:), INTENT(INOUT) :: Reachable
 INTEGER, DIMENSION(:), INTENT(IN) :: IndCon
 INTEGER, INTENT(IN) :: nP, iPi, iPj, NPore

 INTEGER :: higher, lower, i
 INTEGER :: lowerNewInd, change

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
    
    if(NewInd(lower).lt.NewInd(higher)) then
      lowerNewInd = NewInd(lower) 
      change = higher
    elseif(NewInd(higher).lt.NewInd(lower)) then
      lowerNewInd = NewInd(higher) 
      change = lower
    end if
    NewInd(change) = lowerNewInd

    ! Now we need to check if there were some pores already merged with the higher, 
    ! and in that case merge it with the lower
    do i = change, Npore
       if (NewInd(i).eq.change) then
          NewInd(i) = lowerNewInd
       end if
    end do
    
    ! We need to update the Reachable value:
    ! If we merge the pore 3 and 8, if one of the two is accessible, then the merged pore will also be access
    if(Reachable(lower).or.Reachable(higher)) Reachable(lower) = .True.
 end if

end
