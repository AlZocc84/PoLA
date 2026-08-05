 SUBROUTINE Input(dMesh,nR,nP,nAtm,AtmSym,IndCav,IndCon,DiamStep,RMin,MaxDiameter,Closed_Thresh,Print_xyz, &
                 Surf_computation,Rad,n_angles,versors,Accessible)
 
 USE angle_scan

 IMPLICIT NONE
 INTEGER, ALLOCATABLE, INTENT(OUT) :: IndCav(:), IndCon(:)
 INTEGER, INTENT(OUT) :: nP, nAtm, Print_xyz, Surf_computation, n_angles
 INTEGER, DIMENSION(3), INTENT(OUT) :: nR
 INTEGER :: iAtm, lCube, ios
 INTEGER :: iCub, XCub, YCub, ZCub, iCubNew, iC, iM2, iX, iY, iZ, XP, YP, ZP
 INTEGER, PARAMETER :: nTh=10, Scan_Phi(nTh)=[1,1,3,5,6,7,6,5,3,1]
 REAL(DP), PARAMETER :: One=1.0d0, Two=2.0d0, Three=3.0d0, Four=4.0d0
 REAL(DP), PARAMETER :: Pi=Four*Atan(One)
 REAL(DP), PARAMETER :: dRad=Pi/nTh
 REAL(DP), INTENT(OUT) :: dMesh, DiamStep, RMin, Closed_Thresh, MaxDiameter, Rad
 REAL(DP), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: versors
 REAL(DP), ALLOCATABLE :: XAtm(:), YAtm(:), ZAtm(:)
 REAL(DP) :: Xmin, Ymin, Zmin, DX, DY, DZ, dist_x, dist_y, dist_z, D
 CHARACTER(2), ALLOCATABLE, INTENT(OUT) :: AtmSym(:)
 CHARACTER(100) :: File_name, line, print_por, Surf_com, Fix, do_ac_vol
 LOGICAL, INTENT(OUT) :: Accessible
 LOGICAL :: FixCoord

 INTERFACE
   FUNCTION FindRvdW(Symb)
     INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14)
     REAL(DP) :: FindRvdW
     CHARACTER(2), INTENT(IN) :: Symb
   END FUNCTION FindRvdW
   FUNCTION MoveCub(iCub,iX,iY,iZ,nR)
     INTEGER, INTENT(IN) :: iCub,iX,iY,iZ
     INTEGER, DIMENSION(3), INTENT(IN) :: nR
     INTEGER :: MoveCub
   END FUNCTION MoveCub
 END INTERFACE

! Assign default values 
 dMesh = 1.0
 DiamStep = 1.0
 Closed_Thresh = 7.0
 Print_xyz = 0
 Surf_computation = 1
 Rad = 2.0
 FixCoord = .False.
 n_angles = 120
 Accessible = .True.

! Read values from input file 
 open(1, file='input.dat', status='old', form='formatted')
 read(1,'(A)') line
 
 do while(line(1:3).ne.'END')
   if(line(1:1).eq.'#') then

     if(line(3:15).eq.'Molecule_File') then
        read(1,*,iostat=ios) File_name            ! file of coordinates
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if

     else if (line(3:16).eq.'Cube_Mesh_Size') then
        read(1,*,iostat=ios) dMesh                ! size of the cubic grid
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if

     else if(line(3:15).eq.'Out_Mesh_Size') then
        read(1,*,iostat=ios) DiamStep             ! step to output the volumes (usually = dMesh)
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if

     else if(line(3:10).eq.'Box_info') then
        read(1,*,iostat=ios) DX, DY, DZ         ! Cell periodic constants (only orthogonal cells, presently)
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if

     else if(line(3:17).eq.'Threshold_small') then
        read(1,*,iostat=ios) Closed_Thresh ! Threshold on the maximum wall distance to discard small "pores"
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if

     else if(line(3:14).eq.'XYZ_out_file') then
        read(1,*,iostat=ios) Print_por   ! Index to print the xyz file of the cavity (0=don't, 1=does: full/void-surface/void-volume, 2=does w. colors)
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if
        if(Print_por(1:2).eq.'no') then
           Print_xyz = 0
        else if(Print_por(1:5).eq.'basic') then                            ! 1=does: full/void-surface/void-volume
           Print_xyz = 1
        else if(Print_por(1:14).eq.'classification') then
           Print_xyz = 2                                          ! 2= does w. colors : full/ultramicro/.../macro 
        else if(Print_por(1:5).eq.'total') then
           Print_xyz = 3                                          ! 3= does w. color : full/surf-ultramicro/.../surf-macro/vol-ultramicro/.../vol-macro
        else
           write(6,'("Error in file input: value for "A20" is invalid")') line
           stop
        end if

     else if(line(3:21).eq.'Surface_computation') then
             read(1,*,iostat=ios) Surf_com   ! Index to compute the Surface of the cavity (0=doesn't, 1=does)
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if
        if (surf_com(1:2).eq.'no') then
           surf_computation = 0
        else if(surf_com(1:3).eq.'yes') then
           surf_computation = 1
        else
           write(6,'("Error in file input: value for "A20" is invalid")') line
           stop
        end if

     else if(line(3:9).eq.'FixBox') then
       read(1,*,iostat=ios) Fix ! Are the material's coordinates to be centered and wrapped?
       if(Fix(1:3).eq.'yes') FixCoord = .True.
       if(ios.ne.0) then
         write(6,'("Error in file input: value for "A20" is invalid")') line
         stop
        end if
        
     else if(line(3:11).eq.'Probe_Rad') then
        read(1,*,iostat=ios) Rad            ! radius in A of the probe used for the surface computation
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if

     else if(line(3:12).eq.'Angle_scan') then
        read(1,*,iostat=ios) n_angles   ! Number of different direction for angle scan
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if
     
     else if(line(3:15).eq.'Do_Accessible') then
        read(1,*,iostat=ios) do_ac_vol   ! Index to consider only the accessible volumee (no=don't, yes=does)
        if(ios.ne.0) then
                write(6,'("Error in file input: value for "A20" is invalid")') line
                stop
        end if
        if (do_ac_vol(1:2).eq.'no') then
           Accessible = .False.
        else if(do_ac_vol(1:3).eq.'yes') then
           Accessible = .True.
        end if

     else 
        write(6,'("Error: "A20" is invalid")') line 
        stop
     end if
   end if
   read(1,'(A)') line
 end do

 close(1)

 ! chose the versors to map the material:
 allocate(versors(3,n_angles))
 if(n_angles.eq.120) then
    versors = directions120
 elseif(n_angles.eq.180) then
    versors = directions180
 end if

 open(2, file=File_name, status='old', form='formatted')
 read(2,*) nAtm
 read(2,*)
 
 allocate(AtmSym(nAtm),XAtm(nAtm),YAtm(nAtm),ZAtm(nAtm))

 do iAtm = 1, nAtm
   read(2,*) AtmSym(iAtm), XAtm(iAtm), YAtm(iAtm), ZAtm(iAtm)
 end do
 close(2)

 Xmin = minval(XAtm)
 Ymin = minval(YAtm)
 Zmin = minval(ZAtm)

! RMin will be the shortest distance between opposite walls for this solid
 RMin = max(DX, DY, DZ)
! MaxDiameter is the maximum possible size for a "pore" in this solid
 MaxDiameter = max(DX, DY, DZ)

! The following operation is performed by default
! Only in special cases (e.g. when this is part of the analysis of the surface
! coverage during MC) the material's coordinates have to be left unchanged
 if (.not.FixCoord) then
!   Translate the cubic cell so the origin is in the first corner 
!   This has changed (27/09/24): to avoid out-of-range errors in some unfortunate cases, 
!   the first atoms are no longer put in the center of their block
   do iAtm = 1, nAtm
     XAtm(iAtm) = XAtm(iAtm) - Xmin 
     YAtm(iAtm) = YAtm(iAtm) - Ymin 
     ZAtm(iAtm) = ZAtm(iAtm) - Zmin 
!   Check if some coordinate falls outside the given box, meaning that the atom needs to
!   be wrapped inside the box
     if (XAtm(iAtm).gt.DX) XAtm(iAtm) = Mod(XAtm(iAtm),DX)
     if (YAtm(iAtm).gt.DY) YAtm(iAtm) = Mod(YAtm(iAtm),DY)
     if (ZAtm(iAtm).gt.DZ) ZAtm(iAtm) = Mod(ZAtm(iAtm),DZ)
   end do

!   For debugging, if needed print the wrapped coordinates
   open(1,file='wrapped_coord.xyz', status='unknown', form='formatted')
   write(1,*)nAtm
   write(1,*)
   do iAtm = 1, nAtm
     write(1,'(a2,3f12.6)') AtmSym(iAtm),XAtm(iAtm),YAtm(iAtm),ZAtm(iAtm)
   end do
   close(1)
 end if

 nR = (/ NINT(DX/dMesh), NINT(DY/dMesh), NINT(DZ/dMesh) /)
 nP = nR(1)*nR(2)*nR(3)
 iM2 = nR(1)*nR(2)
 allocate(IndCav(nP))
 allocate(IndCon(nP))
! The volume is divided in blocks with edge dMesh. For each block, IdnCav = 0 if void, = 1 if occupied by the
! material skeleton.
 IndCav = 0
 IndCon = 0
! Loop on atoms and find to which "mesh coordinate" each one belongs, then find the block index and put
! the corresponding IndCav to 1 (i.e. this block is filled)
 do iAtm = 1, nAtm
   XCub = INT(XAtm(iAtm)/dMesh) + 1
     if(XCub.eq.(nR(1)+1)) XCub = nR(1) ! Take into account possible rounding errors in INT function 
   YCub = INT(YAtm(iAtm)/dMesh) + 1
     if(YCub.eq.(nR(2)+1)) YCub = nR(2) ! Take into account possible rounding errors in INT function
   ZCub = INT(ZAtm(iAtm)/dMesh) + 1
     if(ZCub.eq.(nR(3)+1)) ZCub = nR(3) ! Take into account possible rounding errors in INT function
   iCub = XCub + nR(1)*(YCub-1) + iM2*(ZCub-1)
   IndCav(iCub) = 1
   IndCon(iCub) = -1

! Consider the vdW radius of this atom and fill the blocks whose center falls inside the vdW sphere 
   lCube = int(FindRvdW(AtmSym(iAtm)) / dMesh)+1
   ! Loop through blocks which could fall inside the vdW sphere
   do iX = -lcube,lcube
    do iY = -lcube,lcube
     do iZ = -lcube,lcube
        iCubNew = MoveCub(iCub, iX, iY, iZ, nR)
        iC = iCubNew - 1
        ! Now getting cell center coordinates
        ZP = INT(iC/(iM2))*dMesh + dMesh/2.0d0
        YP = INT(MOD(iC,iM2)/nR(1))*dMesh + dMesh/2.0d0
        XP = MOD(MOD(iC,iM2),nR(1))*dMesh + dMesh/2.0d0

        ! Now computing distance between centers of blocks and atom taking PBC into account 
        dist_x = XP-XAtm(iAtm)
        dist_y = YP-YAtm(iAtm)
        dist_z = ZP-ZAtm(iAtm)
        dist_x=dist_x-Dx*anint(dist_x/Dx) 
        dist_y=dist_y-Dy*anint(dist_y/Dy)
        dist_z=dist_z-Dz*anint(dist_z/Dz)
        ! Check if center is within radius from atom.

        D = sqrt(dist_x**2+dist_y**2+dist_z**2)
        if (D.le.FindRvdW(AtmSym(iAtm))) then
           IndCav(iCubNew) = 1
           IndCon(iCubNew) = -1
        end if
     end do
    end do
   end do

 end do

 return
 end

