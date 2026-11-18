 SUBROUTINE Output(dMesh,nR,nP,nVol,IndCav,DiamStep,Cumulative_VMinD,VMinD,Accessible, &
             UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV,SkV,SkM,SkD,Print_xyz, &
             RadS,surf_computation,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,nMono_dens,density,molecular_mass,cross_section)
 
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
 INTEGER, DIMENSION(:), INTENT(IN) :: IndCav
 INTEGER, INTENT(IN) :: nP, nVol, Print_xyz, Surf_Computation
 INTEGER, DIMENSION(3), INTENT(IN) :: nR
 INTEGER :: iP, iC, iM2, XCub, YCub, ZCub, iVol
 REAL(DP), INTENT(IN) :: dMesh, DiamStep, SkV, SkM, SkD, nMono_dens
 REAL(DP), INTENT(IN) :: UltraV, MicroV, SmallMesoV, LargeMesoV, MacroV, TotPorV
 REAL(DP), INTENT(IN) :: UltraS, MicroS, SmallMesoS, LargeMesoS, MacroS, RadS, density, molecular_mass, cross_section
 REAL(DP), DIMENSION(:), INTENT(IN) :: Cumulative_VMinD, VMinD
 REAL(DP) :: XX, YY, ZZ, MinD, vMesh, vTot, Fact, Surf_convert, gas_cross
 REAL(DP), PARAMETER :: gMol_A3_to_g_cm3 = 1.6606d0 ! Conversion from g mol-1 A-3 to g cm-3 = 1/0.6022
 REAL(DP), PARAMETER :: A3_gMol_to_cm3_g = 0.6022d0 ! Conversion from A^3/(g mol-1) to cm^3 / g
 REAL(DP), PARAMETER :: A2_to_m2_mol = 6.022d3 ! Conversion from A2/cell to m2/mol
 LOGICAL, INTENT(IN) :: Accessible

 gas_cross = 6.022 * (cross_section/10.0d0) * 1.0d4 ! Conversion from monolayer mol/g to m2/g
 
! Write the results
 if (Print_xyz.ge.1) then
    open(1,file='porous.xyz',status='unknown',form='formatted')
    open(2,file='porous_legend.txt',status='unknown',form='formatted')
 end if
 open(3, file='Cumulative_volume.txt', status='unknown', form='formatted')
 open(4, file='VMinD.txt', status='unknown', form='formatted')
 open(7, file='Texture.txt', status='unknown', form='formatted')
 open(8, file='Simplified_vol.txt', status='unknown', form='formatted')

! If required, print the porous structure, indicating filled, void, "closed" and surface blocks 
! IndCav, Nature of block
!   1               filled                                     --> Material
!   2               excluded because too small                 --> Excluded
!   3               Non Accessible volume                      --> NotAccessible
!   4-5-6-7-8       surf                                       --> Surface
!   9-10-11-12-13   vol                                        --> Bulk
 if (Print_xyz.eq.1) then
   write(2,'("Filled                       --> Material  ")')
   write(2,'("Excluded because too small   --> Excluded  ")')
   write(2,'("Non Accessible volume        --> NotAccessible  ")')
   write(2,'("Surf                         --> Surface  ")')
   write(2,'("Vol                          --> Bulk  ")')
   write(1,*)nP
   write(1,*)
   iM2 = nR(1)*nR(2)
   do iP = 1, nP
     iC = iP - 1
     ZCub = INT(iC/iM2) + 1
     YCub = INT(MOD(iC,iM2)/nR(1)) + 1
     XCub = MOD(MOD(iC,iM2),nR(1)) + 1
     XX = dMesh*(XCub - 0.5)
     YY = dMesh*((YCub + 1) - 0.5)
     ZZ = dMesh*((ZCub + 1) - 0.5)
     if(IndCav(iP).eq.1) write(1,'("MATERIAL",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.2) write(1,'("EXCLUDED",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.3) write(1,'("NOTACC",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).ge.4.AND.IndCav(iP).le.8) write(1,'("SURFACE",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).ge.9.AND.IndCav(iP).le.13) write(1,'("BULK",3f12.6)')XX,YY,ZZ
   end do
 end if

! If a more detailed structure is required, indicate with different elements the nature of the voids
! (in this case, don't signal the surface blocks)
! IndCav, Nature of block
!   1         filled                                     --> MATERIAL
!   2         excluded because too small                 --> MATERIAL
!   3         Non Accessible volume                      --> NOTACC
!   4-9    surf and vol ultramicro (Rmin < 7 A)       --> UMICRO
!   5-10   surf and vol micro (7 < Rmin < 20 A)       --> MICRO
!   6-11   surf and vol small meso (20 < Rmin < 35 A) --> SMESO
!   7-12   surf and vol large meso (35 < Rmin < 50 A) --> LMESO
!   8-13   surf and vol macro (Rmin > 50 A)           --> MACRO

 if (Print_xyz.eq.2) then
   write(2,'("Filled                                     --> MATERIAL   ")')
   write(2,'("Excluded because too small                 --> MATERIAL   ")')
   write(2,'("Non Accessible volume                      --> NOTACC   ")')
   write(2,'("Surf and vol ultramicro (Rmin < 7 A)       --> UMICRO   ")')
   write(2,'("Surf and vol micro (7 < Rmin < 20 A)       --> MICRO   ")')
   write(2,'("Surf and vol small meso (20 < Rmin < 35 A) --> SMESO   ")')
   write(2,'("Surf and vol large meso (35 < Rmin < 50 A) --> LMESO   ")')
   write(2,'("Surf and vol macro (Rmin > 50 A)           --> MACRO   ")')
   write(1,*)nP
   write(1,*)
   iM2 = nR(1)*nR(2)
   do iP = 1, nP
     iC = iP - 1
     ZCub = INT(iC/iM2) + 1
     YCub = INT(MOD(iC,iM2)/nR(1)) + 1
     XCub = MOD(MOD(iC,iM2),nR(1)) + 1
     XX = dMesh*(XCub - 0.5)
     YY = dMesh*((YCub + 1) - 0.5)
     ZZ = dMesh*((ZCub + 1) - 0.5)
     if(IndCav(iP).eq.1) write(1,'("MATERIAL",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.2) write(1,'("MATERIAL",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.3) write(1,'("NOTACC",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.4.OR.IndCav(iP).eq.9.OR.IndCav(iP).eq.14) write(1,'("UMICRO",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.5.OR.IndCav(iP).eq.10.OR.IndCav(iP).eq.15) write(1,'("MICRO",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.6.OR.IndCav(iP).eq.11.OR.IndCav(iP).eq.16) write(1,'("SMESO",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.7.OR.IndCav(iP).eq.12.OR.IndCav(iP).eq.17) write(1,'("LMESO",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.8.OR.IndCav(iP).eq.13.OR.IndCav(iP).eq.18) write(1,'("MACRO",3f12.6)')XX,YY,ZZ
   end do
 end if

! If a complete and detailed structure is required, indicate with different elements the nature of the voids
! in this case, signal the surface blocks
! IndCav , Nature of block
!   1      filled                             --> MATERIAL
!   2      excluded because too small         --> MATERIAL
!   3      Non Accessible volume              --> NOTACC
!   4      surf ultramicro (Rmin < 7 A)       --> UMICROS
!   5      surf micro (7 < Rmin < 20 A)       --> MICROS
!   6      surf small meso (20 < Rmin < 35 A) --> SMESOS 
!   7      surf large meso (35 < Rmin < 50 A) --> SMESOS
!   8      surf macro (Rmin > 50 A)           --> MACROS
!   9      vol ultramicro (Rmin < 7 A)        --> UMICROV
!   10     vol micro (7 < Rmin < 20 A)        --> MICROV
!   11     vol small meso (20 < Rmin < 35 A)  --> CSMESOS
!   12     vol large meso (35 < Rmin < 50 A)  --> SMESOS
!   13     vol macro (Rmin > 50 A)            --> Ba

 if (Print_xyz.eq.3) then
   write(2,'("Filled                             --> MATERIAL   ")')
   write(2,'("Excluded because too small         --> MATERIAL   ")')
   write(2,'("Non Accessible volume              --> NOTACC   ")')
   write(2,'("Surf ultramicro (Rmin < 7 A)       --> UMICROS   ")')
   write(2,'("Surf micro (7 < Rmin < 20 A)       --> MICROS   ")')
   write(2,'("Surf small meso (20 < Rmin < 35 A) --> SMESOS    ")')
   write(2,'("Surf large meso (35 < Rmin < 50 A) --> LMESOS   ")')
   write(2,'("Surf macro (Rmin > 50 A)           --> MACROS   ")')
   write(2,'("Vol ultramicro (Rmin < 7 A)        --> UMICROV   ")')
   write(2,'("Vol micro (7 < Rmin < 20 A)        --> MICROV   ")')
   write(2,'("Vol small meso (20 < Rmin < 35 A)  --> SMESOV   ")')
   write(2,'("Vol large meso (35 < Rmin < 50 A)  --> LMESOV   ")')
   write(2,'("Vol macro (Rmin > 50 A)            --> MACROV   ")')
   write(1,*)nP
   write(1,*)
   iM2 = nR(1)*nR(2)
   do iP = 1, nP
     iC = iP - 1
     ZCub = INT(iC/iM2) + 1
     YCub = INT(MOD(iC,iM2)/nR(1)) + 1
     XCub = MOD(MOD(iC,iM2),nR(1)) + 1
     XX = dMesh*(XCub - 0.5)
     YY = dMesh*((YCub + 1) - 0.5)
     ZZ = dMesh*((ZCub + 1) - 0.5)
     if(IndCav(iP).eq.1) then
        write(1,'("MATERIAL",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.2) then
        write(1,'("MATERIAL",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.3) then
        write(1,'("NOTACC",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.4) then
        write(1,'("UMICROS",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.5) then
        write(1,'("MICROS",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.6) then
        write(1,'("SMESOS",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.7) then
        write(1,'("LMESOS",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.8) then
        write(1,'("MACROS",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.9) then
        write(1,'("UMICROV",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.10) then
        write(1,'("MICROV",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.11) then
        write(1,'("SMESOV",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.12) then
        write(1,'("LMESOV",3f12.6)')XX,YY,ZZ
     elseif(IndCav(iP).eq.13) then
        write(1,'("MACROV",3f12.6)')XX,YY,ZZ
     else
        write(1,'("UK",3f12.6)')XX,YY,ZZ
     end if
   end do
 end if

! The VMinD analysis; Fact converts A3 to cm3/g
 Fact = A3_gMol_to_cm3_g / SkM
 Surf_convert = (A3_gMol_to_cm3_g * density)/molecular_mass

 MinD = 0.0d0
 do iVol = 1, nVol
   MinD = MinD + DiamStep
   write(3,'(f12.4,f14.1,e20.8)') MinD, Cumulative_VMinD(iVol), Fact*Cumulative_VMinD(iVol)
   write(4,'(f12.4,f14.1,e20.8)') MinD, VMinD(iVol), Fact*VMinD(iVol)
 end do

! Summary of the results
 vMesh = dMesh*dMesh*dMesh
 Vtot = nP * vMesh
 write(7,'("======================================================================================")')
 write(7,'("====================== Summary of textural results ===================================")')
 write(7,'("======================================================================================")')
 write(7,'(" Simulation performed with dMesh =   ",f3.1," A")') dMesh
 write(7,'("--------------------------------------------------------------------------------------")')

 write(7,*)
 write(7,'("Solid weight (g/mol)        ",f10.1)') SkM
 write(7,'("Unit Cell volume (A^3)          ",f10.1)') Vtot
 write(7,'("Apparent density (g/cm^3)   ",f10.3)') gMol_A3_to_g_cm3 * SkM/Vtot

 write(7,*)
 write(7,'("Skeletal volume (A^3)       ",f10.1)') SkV
 write(7,'("Skeletal density (g/cm^3)   ",f10.3)') gMol_A3_to_g_cm3 * SkM/SkV

 write(7,*)
 write(7,'("====================== Porous volume analysis ========================================")')

 write(7,*)
 if (Accessible) then
    write(7,'("Porous accessible volume fraction      ",f10.3)') TotPorV/Vtot
    write(7,'("Total accessible porous volume (A^3)   ",f10.1)') TotPorV
    write(7,'("Total accessible porous volume (cm^3/g)",f10.3)') A3_gMol_to_cm3_g * TotPorV/SkM
 else 
    write(7,'("Porous volume fraction      ",f10.3)') TotPorV/Vtot
    write(7,'("Total porous volume (A^3)   ",f10.1)') TotPorV
    write(7,'("Total porous volume (cm^3/g)",f10.3)') A3_gMol_to_cm3_g * TotPorV/SkM
 end if

 write(7,*)
 write(7,'("                                                        with MinD < 7 A    ",f8.3)') &
                                                                    A3_gMol_to_cm3_g * UltraV/SkM
 write(7,'("                                                      /")')
 write(7,'("V(MinD) [cm^3/g], MinD < 20 A   ",f10.3," ---- of which ")') A3_gMol_to_cm3_g * (UltraV+MicroV)/SkM
 write(7,'("                                                      \")')
 write(7,'("                                                        with 7 A < MinD < 20 A      ",f8.3)') &
                                                                    A3_gMol_to_cm3_g * MicroV/SkM
  
 write(7,*)
 write(7,'("                                                        with 20 A < MinD < 35 A     ",f8.3)') &
                                                                    A3_gMol_to_cm3_g * SmallMesoV/SkM
 write(7,'("                                                      /")')
 write(7,'("V(MinD) [cm^3/g], 20 A < MinD < 50 A ",f10.3," ---- of which ")') A3_gMol_to_cm3_g * (SmallMesoV+LargeMesoV)/SkM
 write(7,'("                                                      \")')
 write(7,'("                                                        with 35 A < MinD < 50 A",f8.3)') &
                                                                    A3_gMol_to_cm3_g * LargeMesoV/SkM

 write(7,*)
 write(7,'("V(MinD) [cm^3/g],  MinD > 50 A   ",f10.3)') A3_gMol_to_cm3_g * MacroV/SkM

 write(7,*)                                                                                                                     

 if(Surf_Computation.ge.1) then
   write(7,'("==================== Accessible surface area analysis ================================")')
 end if
 
 if(Surf_Computation.eq.1) then
   write(7,*)

   write(7,'("Total surface (m^2/g) (liquid nitrogen monolayer)     ",f10.1)') &
                                                                           nMono_dens * gas_cross / SkM
   write(7,*)
   write(7,'("                                                        with MinD < 7 A    ",f17.1)') &
                                                                      Surf_convert * UltraS * gas_cross / SkM
   write(7,'("                                                      /")')
   write(7,'("Surface [m^2/g] with MinD < 20 A   ",f17.1,"  ---- of which ")') Surf_convert * (UltraS+MicroS) * gas_cross / SkM 
   write(7,'("                                                      \")')
   write(7,'("                                                        with 7 A < MinD < 20 A      ",f17.1)') &
                                                                      Surf_convert * MicroS * gas_cross / SkM
    
   write(7,*)
   write(7,'("                                                        with 20 A < MinD < 35 A   ",f17.1)') &
                                                                      Surf_convert * SmallMesoS * gas_cross / SkM
   write(7,'("                                                      /")')
   write(7,'("Surface [m^2/g] with 20 < MinD < 50 A ",f17.1," -- of which ")') Surf_convert*(SmallMesoS+LargeMesoS)*gas_cross/SkM
   write(7,'("                                                      \")')
   write(7,'("                                                        with 35 < MinD < 50 A   ",f17.1)') &
                                                                      Surf_convert * LargeMesoS * gas_cross / SkM
                                                                                                                                
   write(7,*)
   write(7,'("Surface [m^2/g] with MinD > 50 A   ",f17.1)') Surf_convert * MacroS * gas_cross / SkM
                                                                                                                                
   write(7,*)                                                                                                                   


  ! ! Write a simple file to plot a bar graph with micro and meso surface
  ! TotS= UltraS + MicroS + SmallMesoS + LargeMesoS + MacroS
  ! write(10,'("# Total, Ultramicro, micro, small meso, large meso and macroporous volumes")')
  ! write(10,'("1 ",f8.3,"  2 ",f8.3,"  3 ",f8.3,"  4 ",f8.3,"  5 ",f8.3,"  6 ",f8.3)') &
  !       A3_gMol_to_cm3_g * TotS/SkM, &
  !       A3_gMol_to_cm3_g * UltraS/SkM, &
  !       A3_gMol_to_cm3_g * MicroS/SkM, &
  !       A3_gMol_to_cm3_g * SmallMesoS/SkM, &
  !       A3_gMol_to_cm3_g * LargeMesoS/SkM, &
  !       A3_gMol_to_cm3_g * MacroS/SkM
 end if 

! Write a simple file to plot a bar graph with micro and meso volumes
 write(8,'("#                    V(MinD) [cm^3/g] with :                              ")')
 write(8,'("# Total | MinD < 7 A | 7 A < MinD < 20 A | 20 A < MinD < 35 A | 35 < MinD < 50 A | MinD > 50 A")')
 write(8,'("1 ",f8.3,"  2 ",f8.3,"  3 ",f8.3,"  4 ",f8.3,"  5 ",f8.3,"  6 ",f8.3)') &
       A3_gMol_to_cm3_g * TotPorV/SkM, &
       A3_gMol_to_cm3_g * UltraV/SkM, &
       A3_gMol_to_cm3_g * MicroV/SkM, &
       A3_gMol_to_cm3_g * SmallMesoV/SkM, &
       A3_gMol_to_cm3_g * LargeMesoV/SkM, &
       A3_gMol_to_cm3_g * MacroV/SkM  

 close(1)
 close(2)
 close(3)
 close(4)
 close(7)
 close(8)

 return
 end

