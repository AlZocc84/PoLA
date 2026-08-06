   FUNCTION FindRvdW(Symb)

     IMPLICIT NONE

     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     REAL(DP) :: FindRvdW
     CHARACTER(2), INTENT(IN) :: Symb

     FindRvdW = 0.9d0
     if (Symb.eq."C")   FindRvdW = 1.7d0
     if (Symb.eq."H")   FindRvdW = 1.2d0
     if (Symb.eq."N")   FindRvdW = 1.55d0
     if (Symb.eq."O")   FindRvdW = 1.52d0
     if (Symb.eq."CH")  FindRvdW = 1.9d0
     if (Symb.eq."Si")  FindRvdW = 1.91d0  ! taken from PoreBlazer UFF Atoms
     if (Symb.eq."ZZ")  FindRvdW = 0.1d0   ! for debug
     if (Symb.eq."Ar")  FindRvdW = 0.5d0   ! for debug
   end

