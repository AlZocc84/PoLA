MODULE gas_parameters

 IMPLICIT NONE

  TYPE GasData
     CHARACTER(len=10) :: name
     REAL(8) :: radius
     REAL(8) :: density
     REAL(8) :: cross_section
     REAL(8) :: molecular_mass
 END TYPE GasData

 TYPE(GasData), PARAMETER :: gases(2) = (/ &
     GasData('N2',  2.0d0, 0.808d0, 16.2d0, 28.0d0), &
     GasData('H2',  1.3d0, 0.202d0, 7.1d0, 2.0d0) &
    /)

 contains

 SUBROUTINE get_gas_parameters(gas, radius, density, cross_section,molecular_mass)

 IMPLICIT NONE

 CHARACTER(len=*), INTENT(IN) :: gas
 REAL(8), INTENT(OUT) :: radius, density, cross_section, molecular_mass

 INTEGER :: i
 LOGICAL :: found

 found = .false.

 do i = 1, size(gases)

    ! adjust() Removes the spaces at the beginning by shifting characters to the left
    ! trim() removes spaces at the end
    ! trim(adjustl()) removes both spaces, at the beginnig and at the end

    if (trim(adjustl(gas)) == trim(gases(i)%name)) then 
       radius = gases(i)%radius
       density = gases(i)%density
       cross_section = gases(i)%cross_section
       molecular_mass = gases(i)%molecular_mass
       
       found = .true.
    end if

 end do

 if(.not.found) then
     write(6,'("ERROR: Unknown gas: ",a," in the input file")') trim(gas)
     stop
 end if

 END SUBROUTINE get_gas_parameters

END MODULE gas_parameters
