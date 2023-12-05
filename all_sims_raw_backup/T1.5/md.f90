! ----------------------------------------------------------------
! Sequential MD program with L-J potential
! By Cheng Zhang, CACS, USC, 2004
! ----------------------------------------------------------------
! The main program
! ----------------------------------------------------------------
program test
   use md_vars
   implicit none
   integer :: i

! Initialize lattice and table

   call initialize()
   call make_tables()
   call buildLlist()
   call cal_intra_force(epot)

! Main body of dynamics loops

   do i = 1, istep
      call time_verlet(i)
      call cal_properties(i)
      if (i .eq. istable) call init_properties()
   enddo

   if (istable .lt. istep) call final_properties()

! Save configuation of atoms for future runs

   call writefiles()

end

! ----------------------------------------------------------------
! This subroutine sets up the lattice and parameters once for all
! ----------------------------------------------------------------
subroutine initialize()
   use md_vars
   implicit none
   real(8) :: currentr(3), boxread(3)
   integer :: i, j, k, l, lcn3, lcn3_p3, lcn3_p6, lcn3_p9
   real(8) :: xlcread, xratio

! First read in all the parameters needed

   open (unit=10, file='md.in', form='formatted', status='old')
   read (10, *) istep, istable, iprint
   read (10, *) ini_new, iscale_x, iscale_v
   read (10, *) req_temp, xlcd
   read (10, *) (ncell(i), i=1, 3)
   read (10, *) isquench, iquench, coeff_quench
   read (10, *) lfreeze
   close (10)


! Examine if the parameter given are too large. If they exceed preset values, reset them

   if (ini_new .eq. 1) then
      do i = 1, 3
         if (ncell(i) .gt. NCMAX) then
            ncell(i) = NCMAX
            print *, 'ncell(', i, ') exceeds NCMAX =', &
               NCMAX, '. Resetting it to NCMAX.'
         endif
      enddo
   endif

   if (istep .gt. NSMAX) then
      istep = NSMAX
      print *, 'istep exceeds NSMAX =', NSMAX, &
         '. Resetting istep to NSMAX.'
   endif

! If we are gonna generate the whole system from lattice positions

   if (ini_new .eq. 1) then

! Set up the dimension parameters of the system
! boxmd contains the sizes of the system in x, y, z respectively
! mc stores the numbers of reduced cells (not unit cell!) in each direction
! and cellsize has the sizes of a single cell

      do i = 1, 3
         boxmd(i) = dble(ncell(i))*xlc(i)*xlcd
      end do

! Set up Argon lattice

      localN = 0
      do i = 1, ncell(3)
         currentr(3) = dble(i - 0.75d0)*xlc(3)*xlcd + 1d-12
         do j = 1, ncell(2)
            currentr(2) = dble(j - 0.75d0)*xlc(2)*xlcd + 1d-12
            do k = 1, ncell(1)
               currentr(1) = dble(k - 0.75d0)*xlc(1)*xlcd + 1d-12
               lcn3 = 3*localN
               lcn3_p3 = lcn3 + 3
               lcn3_p6 = lcn3 + 6
               lcn3_p9 = lcn3 + 9
               do l = 1, 3
                  x(lcn3 + l) = currentr(l)
                  x(lcn3_p3 + l) = currentr(l) + v1(l)*xlcd
                  x(lcn3_p6 + l) = currentr(l) + v2(l)*xlcd
                  x(lcn3_p9 + l) = currentr(l) + v3(l)*xlcd
               enddo
               localN = localN + 4
            enddo
         enddo
      enddo

! Set up random velocities using Maxwell distribution

      call Rand_V()

! Otherwise, we will read in the parameters and coordinates from the file

   else

      open (unit=18, file='xvconf', form='formatted', status='old')
      read (18, *) localN, (boxread(i), i=1, 3), xlcread
      read (18, *) (x(i), i=1, localN*3)
      read (18, *) (v(i), i=1, localN*3)
      close (18)

! Then we calculate those system quantities similar to the other case

      if (iscale_x .eq. 1) then
         xratio = xlcd/xlcread
         do i = 1, 3
            boxmd(i) = boxread(i)*xratio
         enddo
         do i = 1, localN*3
            x(i) = x(i)*xratio
         enddo
      else
         do i = 1, 3
            boxmd(i) = boxread(i)
         enddo
      endif

   endif

   do i = 1, 3
      halfboxmd(i) = 0.5d0*boxmd(i)
      mc(i) = int(boxmd(i)/cutoff)
      cellsize(i) = boxmd(i)/dble(mc(i))
   end do

   volume = boxmd(1)*boxmd(2)*boxmd(3)
   density = dble(localN)/volume
   print *, '# of atoms is', localN, ' Density is ', density
   if ((iscale_x .eq. 1) .and. (ini_new .eq. 0)) then
      print *, ' Lattice constant scaled from', xlcread, ' to', xlcd
   endif

   if ((iscale_v .eq. 1) .or. (ini_new .eq. 1)) then
      call scale_temp(1)
   endif

   if(lfreeze > 0.d0) then
      msd_num_layers = int(boxmd(3)/msd_layer_thickness)
      allocate(msd_by_layer(msd_num_layers,0:1))
      msd_by_layer(:,:)=0.d0
   endif

   call write_xyzfiles(1)
end

!----------------------------------------------------------------------
! This subroutine makes the table, which is built just once and used
! many times in the force calculation to replace real-time calculation
! ---------------------------------------------------------------------
subroutine make_tables()
   use md_vars
   implicit none
   integer :: i
   real(8) :: currentr, ri, ri12, ri2, ri6

! First calculate the offset in potential table and force table

   ri = 1d0/cutoff
   ri2 = ri*ri
   ri6 = ri2**3
   ri12 = ri6*ri6
   voffset = 4d0*(ri12-ri6)
   foffset = (48d0*ri12-24d0*ri6)*ri

! For each sample point, calculate the value of the tables

   do i = 1, NTMAX1 - 1
      currentr = dble(i)*dr
      ri = 1d0/currentr
      ri2 = ri*ri
      ri6 = ri2**3
      ri12 = ri6*ri6
      vtable(i) = 4d0*(ri12-ri6) - voffset + foffset*(currentr - cutoff)
      ftable(i) = (48d0*ri12-24d0*ri6)*ri2 - foffset*ri
   enddo

! Take care of the last point

   vtable(NTMAX1) = 0.d0
   ftable(NTMAX1) = 0.d0

   vtable(NTMAX) = 0.d0
   ftable(NTMAX) = 0.d0

end

! --------------------------------------------------------------
! This subroutine implements the velocity-verlet algorithm
! --------------------------------------------------------------
subroutine time_verlet(is)
   use md_vars
   implicit none
   integer :: is
   integer :: i, i3, k

! Update the velocities and coordinates

   do i = 1, localN
      i3 = i + i + i - 3
      do k = 1, 3
         v(i3 + k) = v(i3 + k) + dthalf*f(i3 + k)
         x(i3 + k) = x(i3 + k) + v(i3 + k)*dt

! Take care of PBC in the case atoms move out of the box.
! Need to wrap them back and also update the x reference array for
! calculating mean square displacement purpose.

         if (x(i3 + k) .lt. 0d0) then
            x(i3 + k) = x(i3 + k) + boxmd(k)
            originx(i3 + k) = originx(i3 + k) + boxmd(k)
         elseif (x(i3 + k) .gt. boxmd(k)) then
            x(i3 + k) = x(i3 + k) - boxmd(k)
            originx(i3 + k) = originx(i3 + k) - boxmd(k)
         endif
      enddo
   enddo

! Update the linklist

   if (mod(is, 1) .eq. 0) then
      call buildLlist()
   endif

! Calculate the forces and update velocities again to complete the step

   call cal_intra_force(epot)

   do i = 1, localN*3
      v(i) = v(i) + dthalf*f(i)
   enddo

   if ((isquench .eq. 1) .and. (mod(is, iquench) .eq. 1)) then
      call scale_temp(2)
   endif

! Freeze atoms lfreeze unit from bottom 
   do i=1, localN
      if(x(i*3)<lfreeze) v(i*3-2:i*3)=0.d0
   enddo

   return
end

!---------------------------------------------------------------------
! Set up the random velocities
! --------------------------------------------------------------------
subroutine Rand_V()
   use md_vars
   implicit none
   integer :: i, i3, ia
   real(8) :: vav(3), vbuf(3), facv, rnd1, rnd2, temp_new, twopi

! Set up constants and initialize

   twopi = 2d0*Pi
   temp_new = 0.3d0

! Prepare Maxwellian velocities

   facv = dsqrt(3d0*temp_new)
   dseed = 13579d0

! Physical velocity. The Maxwell-Boltzmann distribution is
! f(v) = c * v * v * exp(- v * v / T) in normalized units. So we can
! work the opposite way to generate velocities with M-B distribution using
! random numbers from 0 to 1.
! Here facv is the factor to scale velocities to the desired temperature.
! dcos() is the angular part (assuming even angular distribution) and
! the dsqrt() is the magnitude part that follows Gaussian distribution.

   do i3 = 1, 3*localN
      call drnd(rnd2)
      call drnd(rnd1)
      v(i3) = facv*dsqrt(-dlog(rnd1))*dcos(twopi*rnd2)
   enddo

! Make total momentum zero
! Calculate the global average momentum & center-of-mass

   do ia = 1, 3
      vav(ia) = 0d0
      vbuf(ia) = 0d0
   enddo

   do i = 1, localN
      do ia = 1, 3
         vav(ia) = vav(ia) + v(3*i - 3 + ia)
      enddo
   enddo

   do ia = 1, 3
      vav(ia) = vav(ia)/dble(LocalN)
   enddo

! Shift velocities to make the total momentum zero

   do i = 1, localN
      do ia = 1, 3
         v(3*i - 3 + ia) = v(3*i - 3 + ia) - vav(ia)
      enddo
   enddo

   return
end

!----------------------------------------------------------------
! Random-number generator
!----------------------------------------------------------------
subroutine drnd(rnd)
   use md_vars
   implicit none
   real(8) :: rnd
   real(8), parameter :: d2p31m = 2147483647d0
   real(8), parameter :: d2p31 = 2147483648d0

! Generate a new random number

   dseed = dmod(16807d0*dseed, d2p31m)
   rnd = dseed/d2p31
   dseed = dseed + 1d0

   return
end

!----------------------------------------------------------------
! Scale the temperature to the preset value
!----------------------------------------------------------------
subroutine scale_temp(itype)
   use md_vars
   implicit none
   integer :: i, i3, itype
   real(8) :: scale, temp

! Calculate the current temperature
! Find out the ratio and scale the velocities

   if (itype .eq. 1) then
      ekin = 0d0
      do i = 1, localN*3
         ekin = ekin + v(i)*v(i)
      enddo
      temp = ekin/dble(3*localN)
      scale = dsqrt(req_temp/temp)
      print *, 'Temperature scaled to', req_temp
   elseif (itype .eq. 2) then
      scale = dsqrt(coeff_quench)
      print *, 'Temerature scaled by', coeff_quench
   endif

   do i = 1, localN*3
      v(i) = v(i)*scale
   enddo

   return
end

! -------------------------------------------------------------------
! This subroutine builds up the linklist for future force calculation
! -------------------------------------------------------------------
subroutine buildLlist()
   use md_vars
   implicit none
   integer :: i, j, k, i3, lcoor(3)

! Initialize the linklist

   do i = 1, mc(1)
      do j = 1, mc(2)
         do k = 1, mc(3)
            nheader(i, j, k) = 0
            lsize(i, j, k) = 0
         enddo
      enddo
   enddo

   do i = 1, NMAX
      linklist(i) = 0
   enddo

! Decide which cell an atom belongs to and put it in the linklist

   do i = 1, localN
      i3 = (i - 1)*3
      do j = 1, 3
         lcoor(j) = int(x(i3 + j)/cellsize(j)) + 1
      enddo
      linklist(i) = nheader(lcoor(1), lcoor(2), lcoor(3))
      nheader(lcoor(1), lcoor(2), lcoor(3)) = i
      lsize(lcoor(1), lcoor(2), lcoor(3)) = lsize(lcoor(1), lcoor(2), lcoor(3)) + 1

   enddo

end

! ---------------------------------------------------------------------
! This subroutine calculated the forces on atoms and update the f array
! ---------------------------------------------------------------------
subroutine cal_intra_force(vsum)
   use md_vars
   implicit none
   real(8) :: vsum
   real(8) :: dx(3), r, fr, vr, rsq, rdiv
   integer :: i, j, k, m, n, ni, nj, ixm, iym, izm, ip, i3, j3, idiv
   integer :: ix, iy, iz, jx, jy, jz

! Set accumulators to zero

   do i = 1, NMAX3
      f(i) = 0.0
   enddo
   vsum = 0d0

! Loop over cells

!$OMP parallel do collapse(3) default(shared) &
!$OMP private(i,j,k,m,n,ni,nj,ixm,iym,izm,ip,i3,j3,idiv) &
!$OMP private(ix,iy,iz,jx,jy,jz,r,fr,vr,rsq,rdiv,dx) &
!$OMP reduction(+:vsum,pressure)
   do ix = 1, mc(1)
      do iy = 1, mc(2)
         do iz = 1, mc(3)

! First find particles in cell (ix, iy, iz)

            i = nheader(ix, iy, iz)

            do ni = 1, lsize(ix, iy, iz)
               i3 = i + i + i - 3

! Find particles in neighbouring cells

               do jx = -1, 1
               do jy = -1, 1
               do jz = -1, 1

                  ixm = ix + jx
                  iym = iy + jy
                  izm = iz + jz

! Periodic boundary conditions in x, y and z

                  if (ixm .eq. mc(1) + 1) then
                     ixm = 1
                  else if (ixm .eq. 0) then
                     ixm = mc(1)
                  endif
                  if (iym .eq. mc(2) + 1) then
                     iym = 1
                  else if (iym .eq. 0) then
                     iym = mc(2)
                  endif
                  if (izm .eq. mc(3) + 1) then
                     izm = 1
                  else if (izm .eq. 0) then
                     izm = mc(3)
                  endif

! For each cell, save local atoms and neighboring atoms in an array

                  j = nheader(ixm, iym, izm)
                  do nj = 1, lsize(ixm, iym, izm)
                     if (i < j) then
                        j3 = j + j + j - 3

! Calculate distances between particles in the array

                        rsq = 0d0

! Calculate distances after applying periodic boundary conditions

                        do k = 1, 3
                           dx(k) = x(i3 + k) - x(j3 + k)
                           if (dx(k) .gt. halfboxmd(k)) dx(k) = dx(k) - boxmd(k)
                           if (dx(k) .lt. -halfboxmd(k)) dx(k) = dx(k) + boxmd(k)
                           rsq = rsq + dx(k)*dx(k)
                        enddo

! If that distance is smaller than the cutoff, proceed to force calculation
! using that distance

                        if (rsq < cutoffsq) then

                           r = dsqrt(rsq)
                           rdiv = r*dri
                           idiv = int(rdiv)
                           rdiv = rdiv - idiv
                           vr = vtable(idiv) + rdiv*(vtable(idiv + 1) - vtable(idiv))
                           fr = ftable(idiv) + rdiv*(ftable(idiv + 1) - ftable(idiv))

                           vsum = vsum + vr
                           pressure = pressure + rsq*fr

! Store force in the carrier array temporarily

                           do k = 1, 3
!$OMP atomic
                              f(i3 + k) = f(i3 + k) + dx(k)*fr
!$OMP atomic
                              f(j3 + k) = f(j3 + k) - dx(k)*fr
                           enddo

                        endif

                     endif

                     j = linklist(j)
                  enddo
               enddo; enddo; enddo

               i = linklist(i)
            enddo

         enddo; enddo; enddo
!$OMP end parallel do

   return
end

! -------------------------------------------------------------------
! This subroutine initializes the variables needed for calculating properties
! -------------------------------------------------------------------
subroutine init_properties()
   use md_vars
   implicit none
   integer :: i

! Initializing and calculate quantities used later

   do i = 1, localN*3
      originx(i) = x(i)
   enddo

   akinsq = 0d0
   akin = 0d0
   pressure = 0d0

! Determine the cutoff and grid size of g(r)

   grcutoff = xlcd*NUMXLCDFORGR
   if (grcutoff .gt. halfboxmd(2)) grcutoff = halfboxmd(2)
   if (grcutoff .gt. halfboxmd(3)) grcutoff = halfboxmd(3)
   grcutoffsq = grcutoff*grcutoff
   gdr = grcutoff/dble(NGMAX)

   print'(a,f10.5)', 'grcutoff : ', grcutoff

! Initialize variables

   do i = 1, NGMAX
      g(i) = 0d0
      cnum(i) = 0d0
   enddo

   return
end

! ----------------------------------------------------------------------
! Calculate the properties of the system
! ----------------------------------------------------------------------
subroutine cal_properties(is)
   use md_vars
   implicit none
   integer :: i, is, isp
   real(8) :: total_diff

   integer :: idx
   real(8) :: rr(3)
   logical :: isFileOpen

! Calculate kinetic energy

   ekin = 0d0
   do i = 1, localN*3
      ekin = ekin + v(i)*v(i)
   enddo

   ekin = 0.5d0*ekin/dble(localN)
   epot = epot/dble(localN)

! Calculate total energy

   etot = ekin + epot
   if (mod(is, iprint) .eq. 0) then
      isp = is/iprint
      energy(isp, 1) = etot
      energy(isp, 2) = epot
      energy(isp, 3) = ekin
   endif

! If the system is considered to be thermolized, then start calculating Cv
! and diffusion constant. First calculate diffusion for the current step,
! then accumulate akin and akinsq for Cv calculation. At last call
! cal_gofr() subroutine to accumulate g(r) and coordination number array.

   if (is .gt. istable) then
      total_diff = 0d0
      do i = 1, localN*3
         total_diff = total_diff + (x(i) - originx(i))**2
      enddo
      total_diff = total_diff/dble(localN)
      if (mod(is - istable, 20) .eq. 0) then
         diffusion((is - istable)/20) = total_diff
      endif
      akin = akin + ekin
      akinsq = akinsq + ekin*ekin
      call cal_gofr()
      pressure = pressure + 2d0*ekin*dble(localN)

      if(lfreeze > 0.d0) then
         msd_by_layer(:,:)=0.d0
         do i=1, localN
            idx = int(x(i*3)/msd_layer_thickness)+1
            rr(1:3)=x(i*3-2:i*3)-originx(i*3-2:i*3)
            msd_by_layer(idx,0)=msd_by_layer(idx,0)+1
            msd_by_layer(idx,1)=msd_by_layer(idx,1)+sum(rr(1:3)*rr(1:3))
         enddo
      endif

   endif

   if (mod(is, iprint) .eq. 0) then
      if (is .le. istable) then
         write (*, '(I6, A, F11.7, A, F11.6, A, F11.6)') &
            is, ' E=', etot, ' PE=', epot, ' KE=', ekin
      else
         write (*, '(I6, A, F11.7, A, F11.6, A, F11.6, A, F11.6)') &
            is, ' E=', etot, ' PE=', epot, ' KE=', ekin, &
            ' Meansq displacement=', total_diff

         if(lfreeze > 0.d0) then
            inquire(file="msd.out", opened=isFileOpen)
            if(.not. isFileOpen)  open(42,file="msd.out")
            write(*,'(a5, i6 $)') 'msd: ', is
            write(42,'(a5, i6 $)') 'msd: ', is
            do idx = 2, msd_num_layers
               if ( msd_by_layer(idx,0) > 30.d0 ) then
                  write(*,'(f10.5 $)') msd_by_layer(idx,1)/msd_by_layer(idx,0)
                  write(42,'(f10.5 $)') msd_by_layer(idx,1)/msd_by_layer(idx,0)
               else
                  write(*,'(f10.5 $)') 0.d0
                  write(42,'(f10.5 $)') 0.d0
               endif
            enddo 
            write(*,'(a3 $)') ' - samples '
            write(42,'(a3 $)') ' - samples '
            do idx = 2, msd_num_layers
               write(*,'(i5 $)') int(msd_by_layer(idx,0))
               write(42,'(i5 $)') int(msd_by_layer(idx,0))
            enddo 
            write(*,*)
            write(42,*)
         endif
      endif

! Save atom config into xyz file
      call write_xyzfiles(0)
   endif

   return
end

! -------------------------------------------------------------------
! This subroutine calculates g(r)
! -------------------------------------------------------------------
subroutine cal_gofr()
   use md_vars
   implicit none
   real(8) :: dx(3), gt(NGMAX), r, rsq, radius
   integer :: i, j, k, i3, j3, idiv, ix, iy, iz, jx, jy, jz, ixm, iym, izm, m, ni, nj

! Initialize the array for storage

   do i = 1, NGMAX
      gt(i) = 0d0
   enddo

! Loop through pairs in the system

   do i = 1, localN - 1
      do j = i + 1, localN
         i3 = i * 3 - 3
         j3 = j * 3 - 3
         rsq = 0d0

! Calculate distances after applying periodic boundary conditions

         do k = 1, 3
            dx(k) = x(i3 + k) - x(j3 + k)
            if(dx(k) .gt. halfboxmd(k)) dx(k) = dx(k) - boxmd(k)
            if(dx(k) .lt. -halfboxmd(k)) dx(k) = dx(k) + boxmd(k)
            rsq = rsq + dx(k) * dx(k)
         enddo

! Find out which slot this pair distance belongs to

         if (rsq .lt. grcutoffsq) then
            r = dsqrt(rsq)
            idiv = int(r/gdr) + 1
            gt(idiv) = gt(idiv) + 1
         endif

      enddo
   enddo

! Normalizing

   do i = 1, NGMAX
      radius = dble(i)*gdr
      gt(i) = gt(i)/(2d0*Pi*radius*radius*gdr*dble(localN)*density)
      g(i) = g(i) + gt(i)
   enddo

   return
end

! -------------------------------------------------------------------
! This subroutine calculates the final properties at the end of run
! -------------------------------------------------------------------
subroutine final_properties()
   use md_vars
   implicit none
   integer :: i, j, isp
   real(8) :: atemp, radius, Cv, silo

! Saving energy data

   open (unit=12, file='epotke', form='formatted', status='unknown')
   do i = 1, istep/iprint
      write (12, '(3F12.7)') (energy(i, j), j=1, 3)
   enddo
   close (12)

!3000    format (3F12.7)

! Calculating Cv using fluctuation of kinetic energy

   akin = akin/dble(istep - istable)
   atemp = akin/1.5d0
   akinsq = akinsq/dble(istep - istable)
   akin = akin*akin
   silo = 2d0/(3d0*localN) - (akinsq - akin)/akin
   Cv = 1d0/(silo*dble(localN))
   print *, 'T = ', atemp, ' Cv = ', Cv

! Saving the diffusion data

   open (unit=14, file='meansqdisp', form='formatted', status='unknown')
   do i = 1, (istep - istable)/20
      write (14, *) dble(i*20)*dt, diffusion(i)
   enddo
   close (14)

! Calculate coordination number by integrating g(r)
! Saving g(r) and coor num results in files called gr and nr

   open (unit=15, file='gr', form='formatted', status='unknown')
   open (unit=16, file='nr', form='formatted', status='unknown')
   do i = 1, NGMAX
      radius = dble(i)*gdr
      g(i) = g(i)/dble(istep - istable)
      if (i .gt. 1) cnum(i) = cnum(i - 1) + g(i) &
                              *4d0*Pi*radius*radius*gdr*density
      write (15, '(2F12.6)') radius, g(i)
      write (16, '(2F12.6)') radius, cnum(i)
   enddo
   close (15)
   close (16)

! Calculating pressure using average

   pressure = pressure/(3d0*(istep - istable)*volume)
   print *, 'Pressure = ', pressure

! Saving properties in a summary file

   open (unit=17, file='properties', form='formatted', status='unknown')
   write (17, *) 'Cv = ', Cv, ' T = ', atemp
   write (17, *) 'Pressure = ', pressure, ' Energy = ', etot
   do i = 1, 5
      isp = i*(istep - istable)/100
      write (17, *) 'Mean square displacement at', isp, &
         'th step is', diffusion(isp)
   enddo
   close (17)

   return
end

! -------------------------------------------------------------------
! This subroutine writes the configuration data to disk
! -------------------------------------------------------------------
subroutine writefiles()
   use md_vars
   implicit none
   integer :: i

   open (unit=18, file='xvconf', form='formatted', status='unknown')
   write (18, *) localN, (boxmd(i), i=1, 3), xlcd
   write (18, *) (x(i), i=1, localN*3)
   write (18, *) (v(i), i=1, localN*3)
   close (18)

   return
end

! -------------------------------------------------------------------
! This subroutine writes atom data in XYZ format to disk
! -------------------------------------------------------------------
subroutine write_xyzfiles(inew)
   use md_vars
   implicit none
   integer,intent(in) :: inew
   integer :: i

   if(inew==1) then
      open (unit=18, file='md.xyz', form='formatted', status='unknown')
   else
      open (unit=18, file='md.xyz', form='formatted', position='append',status='old')
   endif

   write (18, '(i6)') localN
   write (18, '(4f12.5)') (boxmd(i), i=1, 3), xlcd
   do i=1, localN
      write (18, '(a3,7f8.3)') 'Ar', x(i*3-2:i*3), originx(i*3-2:i*3), &
            sum(v(i*3-2:i*3)*v(i*3-2:i*3))
   enddo
   close (18)

   return
end
