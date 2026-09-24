 SUBROUTINE Output(dMesh,nR,nP,nVol,IndCav,DiamStep,Cumulative_VMinD,VMinD,Accessible, &
             UltraV,MicroV,SmallMesoV,LargeMesoV,MacroV,TotPorV,SkV,SkM,SkD,Print_xyz, &
             Rad,surf_computation,UltraS,MicroS,SmallMesoS,LargeMesoS,MacroS,nMono_dens)
 
 IMPLICIT NONE
 INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
 INTEGER, DIMENSION(:), INTENT(IN) :: IndCav
 INTEGER, INTENT(IN) :: nP, nVol, Print_xyz, Surf_Computation, nMono_dens
 INTEGER, DIMENSION(3), INTENT(IN) :: nR
 INTEGER :: iP, iC, iM2, XCub, YCub, ZCub, iVol
 REAL(DP), INTENT(IN) :: dMesh, DiamStep, SkV, SkM, SkD
 REAL(DP), INTENT(IN) :: UltraV, MicroV, SmallMesoV, LargeMesoV, MacroV, TotPorV
 REAL(DP), INTENT(IN) :: UltraS, MicroS, SmallMesoS, LargeMesoS, MacroS, Rad
 REAL(DP), DIMENSION(:), INTENT(IN) :: Cumulative_VMinD, VMinD
 REAL(DP) :: XX, YY, ZZ, MinD, vMesh, vTot, Fact, Fact_sup, Surf_convert
 REAL(DP), PARAMETER :: gMol_A3_to_g_cm3 = 1.6606d0 ! Conversion from g mol-1 A-3 to g cm-3 = 1/0.6022
 REAL(DP), PARAMETER :: A3_gMol_to_cm3_g = 0.6022d0 ! Conversion from A^3/(g mol-1) to cm^3 / g
 REAL(DP), PARAMETER :: A2_to_m2_mol = 6.022d3 ! Conversion from A2/cell to m2/mol
 REAL(DP), PARAMETER :: N2_cross = 6.022 * 1.62 * 1.0d4 ! Conversion from monolayer mol/g to m2/g
 LOGICAL, INTENT(IN) :: Accessible

! Write the results
 if (Print_xyz.ge.1) open(1,file='porous.xyz',status='unknown',form='formatted')
 open(2, file='Cumulative_volume.txt', status='unknown', form='formatted')
 open(4, file='VMinD.txt', status='unknown', form='formatted')
 open(7, file='Texture.txt', status='unknown', form='formatted')
 open(8, file='Simplified_vol.txt', status='unknown', form='formatted')
 !open(10, file='Simplified_surface.txt', status='unknown', form='formatted')

! If required, print the porous structure, indicating filled, void, "closed" and surface blocks 
! IndCav, Nature of block
!   1               filled                                     --> Ar
!   2               excluded because too small                 --> Kr
!   3               Non Accessible volume                      --> He
!   4-5-6-7-8       surf                                       --> Ne
!   9-10-11-12-13   vol                                        --> Xe
 if (Print_xyz.eq.1) then
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
     if(IndCav(iP).eq.1) write(1,'("Ar",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.2) write(1,'("Kr",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.3) write(1,'("He",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).ge.4.AND.IndCav(iP).le.8) write(1,'("Ne",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).ge.9.AND.IndCav(iP).le.13) write(1,'("Xe",3f12.6)')XX,YY,ZZ
   end do
 end if

! If a more detailed structure is required, indicate with different elements the nature of the voids
! (in this case, don't signal the surface blocks)
! IndCav, Nature of block
!   1         filled                                     --> Ar
!   2         excluded because too small                 --> Ar
!   3         Non Accessible volume                      --> He
!   4-9    surf and vol ultramicro (Rmin < 7 A)       --> Ne
!   5-10   surf and vol micro (7 < Rmin < 20 A)       --> Kr
!   6-11   surf and vol small meso (20 < Rmin < 35 A) --> Xe
!   7-12   surf and vol large meso (35 < Rmin < 50 A) --> Rn
!   8-13   surf and vol macro (Rmin > 50 A)           --> Ca

 if (Print_xyz.eq.2) then
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
     if(IndCav(iP).eq.1) write(1,'("Ar",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.2) write(1,'("Ar",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.3) write(1,'("He",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.4.OR.IndCav(iP).eq.9.OR.IndCav(iP).eq.14) write(1,'("Ne",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.5.OR.IndCav(iP).eq.10.OR.IndCav(iP).eq.15) write(1,'("Kr",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.6.OR.IndCav(iP).eq.11.OR.IndCav(iP).eq.16) write(1,'("Xe",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.7.OR.IndCav(iP).eq.12.OR.IndCav(iP).eq.17) write(1,'("Rn",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.8.OR.IndCav(iP).eq.13.OR.IndCav(iP).eq.18) write(1,'("Ca",3f12.6)')XX,YY,ZZ
   end do
 end if

! If a complete and detailed structure is required, indicate with different elements the nature of the voids
! in this case, signal the surface blocks
! IndCav , Nature of block
!   1      filled                             --> Ar
!   2      excluded because too small         --> Ar
!   3      Non Accessible volume              --> Kr
!   4      surf ultramicro (Rmin < 7 A)       --> Li
!   5      surf micro (7 < Rmin < 20 A)       --> Na
!   6      surf small meso (20 < Rmin < 35 A) --> K 
!   7      surf large meso (35 < Rmin < 50 A) --> Rb
!   8      surf macro (Rmin > 50 A)           --> Cs
!   9      vol ultramicro (Rmin < 7 A)        --> Be
!   10     vol micro (7 < Rmin < 20 A)        --> Mg
!   11     vol small meso (20 < Rmin < 35 A)  --> Ca
!   12     vol large meso (35 < Rmin < 50 A)  --> Sr
!   13     vol macro (Rmin > 50 A)            --> Ba

 if (Print_xyz.eq.3) then
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
     if(IndCav(iP).eq.1) write(1,'("Ar",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.2) write(1,'("Ar",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.3) write(1,'("Kr",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.4) write(1,'("Li",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.5) write(1,'("Na",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.6) write(1,'("K",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.7) write(1,'("Rb",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.8) write(1,'("Cs",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.9) write(1,'("Be",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.10) write(1,'("Mg",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.11) write(1,'("Ca",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.12) write(1,'("Sr",3f12.6)')XX,YY,ZZ
     if(IndCav(iP).eq.13) write(1,'("Ba",3f12.6)')XX,YY,ZZ
   end do
 end if

! The VMinD analysis; Fact converts A3 to cm3/g
 Fact = A3_gMol_to_cm3_g / SkM
 Fact_sup= 3.570795*Rad*Rad
 Surf_convert = (A3_gMol_to_cm3_g * 0.808)/28

 MinD = 0.0d0
 do iVol = 1, nVol
   MinD = MinD + DiamStep
   write(2,'(f12.4,f14.1,e20.8)') MinD, Cumulative_VMinD(iVol), Fact*Cumulative_VMinD(iVol)
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
   write(7,'("Total surface (m^2/g) (liquid nitrogen monolayer)   ",f10.1)') nMono_dens * N2_cross / SkM
   write(7,*)
   write(7,'("                                                        with MinD < 7 A    ",f8.3)') &
                                                                      Surf_convert * UltraS * N2_cross / SkM
   write(7,'("                                                      /")')
   write(7,'("Surface [m^2/g] with MinD < 20 A   ",f10.3," ---- of which ")') Surf_convert * (UltraS+MicroS) * N2_cross / SkM
   write(7,'("                                                      \")')
   write(7,'("                                                        with 7 A < MinD < 20 A      ",f8.3)') &
                                                                      Surf_convert * MicroS * N2_cross / SkM
    
   write(7,*)
   write(7,'("                                                        with 20 A < MinD < 35 A   ",f8.3)') &
                                                                      Surf_convert * SmallMesoS * N2_cross / SkM
   write(7,'("                                                      /")')
   write(7,'("Surface [m^2/g] with 20 < MinD < 50 A ",f10.3," ---- of which ")') Surf_convert*(SmallMesoS+LargeMesoS)*N2_cross/SkM
   write(7,'("                                                      \")')
   write(7,'("                                                        with 35 < MinD < 50 A   ",f8.3)') &
                                                                      Surf_convert * LargeMesoS * N2_cross / SkM
                                                                                                                                
   write(7,*)
   write(7,'("Surface [m^2/g] with MinD > 50 A   ",f10.3)') Surf_convert * MacroS * N2_cross / SkM
                                                                                                                                
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
 close(10)

 return
 end

