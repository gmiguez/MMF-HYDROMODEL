MODULE module_forcings

   use module_parallel
   use interp_lib

   implicit none
  integer, parameter :: dimerax=250,dimeray=290,xdim=1440,ydim=721,xdimeraland=3600,ydimeraland=1801,dimeralandx=620,dimeralandy=720
   real, parameter :: eraswlat = -89.46282, eraswlon = 0., dxera = 0.703125

CONTAINS

!********************************************************************************

   SUBROUTINE READFORCINGS(n2, js, je, filename, irec, hour, topo, lats, lons &
                           , vels, temp, pres, qair, landmask, topoera5, varpack, request, istart)
      use netcdf

      integer :: n2, js, je, i, j, ii, jj, n, istart, hour
      integer :: domblocke2obs, request, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2)
      character*50 :: filename
      character*200 :: filee2obs
      integer, dimension(n2, js:je) :: landmask
      real, dimension(n2, js:je) :: topo, topoera5, lats, lons, vels, temp, pres, qair
      real :: tempint, elev, relh
      integer*2, allocatable, dimension(:, :, :) :: varread
      real, allocatable, dimension(:,:,:) :: dtemp,press,pressesat,pressesattd,rh,tempera,wind,rain,swdown,rnetera,st1,st2,st3,st4
      real, dimension(dimerax, dimeray, 0:23, 11) :: varpack
      integer :: irec
      real :: dx, dy, grx, gry, xswlat, xswlon, xlon, wt1, wt2, wt3, wt4
      double precision :: scale_factor, add_offset

      integer :: ncid, varid, statusnc
      integer :: start(3), count(3), startp(3), countp(3)

      call MPI_Type_CONTIGUOUS(24*dimerax*dimeray, MPI_REAL, domblocke2obs, ierr)
      call MPI_Type_commit(domblocke2obs, ierr)

      if (hour .eq. 0) then
      IF (pid .eq. 0) THEN

!if(hour.eq.0)then

         allocate (varread(xdim, ydim, 0:23))
         allocate (wind(xdim, ydim, 0:23))
         allocate (dtemp(xdim, ydim, 0:23))
         allocate (press(xdim, ydim, 0:23))
         allocate (pressesat(xdim, ydim, 0:23))
         allocate (pressesattd(xdim, ydim, 0:23))
         allocate (tempera(xdim, ydim, 0:23))
         allocate (rh(xdim, ydim, 0:23))

         count(1) = xdim
         count(2) = ydim
         count(3) = 24
         start(1) = 1
         start(2) = 1
         start(3) = irec

         countp(1) = xdim
         countp(2) = ydim
         countp(3) = 1
         startp(1) = 1
         startp(2) = 1
         startp(3) = irec

!Wind velocity

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/WIND/ERA5_wind_speed_'//filename(1:len_trim(filename))

!write(6,*)'reading wind forcing from ',filee2obs

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         call readvar(ncid, xdim, ydim, 'ws10', start, count, varread, scale_factor, add_offset)

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         wind = dble(varread)*scale_factor + add_offset
!stop

!Temperature

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/TEMP/ERA5_2m_temperature_'//filename(1:len_trim(filename))

!write(6,*)'reading forcing from ',filee2obs

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         call readvar(ncid, xdim, ydim, 't2m', start, count, varread, scale_factor, add_offset)

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         tempera = dble(varread)*scale_factor + add_offset
         pressesat = exp(17.67*(tempera - 273.15)/(243.5 + (tempera - 273.15)))!*611.2 !in Pa

!write(6,*)'finished with Tair',pid

!Pressure

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/SFCPRESS/ERA5_surface_pressure_'//filename(1:len_trim(filename))

!write(6,*)'reading forcing from ',filee2obs

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         call readvar(ncid, xdim, ydim, 'sp', start, count, varread, scale_factor, add_offset)

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         press = dble(varread)*scale_factor + add_offset !in Pa

!write(6,*)'finished with PSurf',pid

!Mixing ratio to rh

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/DEWTEMP/ERA5_2m_dewpoint_temperature_'//filename(1:len_trim(filename))

!write(6,*)'reading forcing from ',filee2obs

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         call readvar(ncid, xdim, ydim, 'd2m', start, count, varread, scale_factor, add_offset)

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         dtemp = dble(varread)*scale_factor + add_offset
         pressesattd = exp(17.67*(dtemp - 273.15)/(243.5 + (dtemp - 273.15)))!*611.2
         rh = min(max(pressesattd/pressesat, 0.), 1.)

         deallocate (dtemp, pressesat, pressesattd)

         allocate (rain(xdim, ydim, 0:23))
         allocate (swdown(xdim, ydim, 0:23))
         allocate (rnetera(xdim, ydim, 0:23))

!write(6,*)'finished with RH',pid

!Precipitation

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/PRECIP/ERA5_precip_'//filename(1:len_trim(filename))

!write(6,*)'reading from ',filee2obs,irec

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         call readvar(ncid, xdim, ydim, 'tp', start, count, varread, scale_factor, add_offset)

         rain = dble(varread)*scale_factor + add_offset

         call readvar(ncid, xdim, ydim, 'sf', start, count, varread, scale_factor, add_offset)

         rain = max(rain - (dble(varread)*scale_factor + add_offset), 0.)*1000. !in mm

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!write(6,*)'read rain'

!Shortwave radiation

    filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/RADIATION/ERA5_surface_solar_radiation_downwards_'//filename(1:len_trim(filename))

!write(6,*)'reading from ',filee2obs,irec

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         call readvar(ncid, xdim, ydim, 'ssrd', start, count, varread, scale_factor, add_offset)

         swdown = dble(varread)*scale_factor + add_offset

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!write(6,*)'read SW'

!net radiation

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/RADIATION/ERA5_surface_net_radiation_'//filename(1:len_trim(filename))

!write(6,*)'reading from ',filee2obs,irec

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         call readvar(ncid, xdim, ydim, 'sr', start, count, varread, scale_factor, add_offset)

         rnetera = dble(varread)*scale_factor + add_offset

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!soil levels

         allocate (st1(xdim, ydim, 0:23))
         allocate (st2(xdim, ydim, 0:23))
         allocate (st3(xdim, ydim, 0:23))
         allocate (st4(xdim, ydim, 0:23))

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/SOILTEMP/ERA5_soiltemp_'//filename(1:len_trim(filename))

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!write(6,*)'reading from ',filee2obs,irec

         call readvar(ncid, xdim, ydim, 'stl1', start, count, varread, scale_factor, add_offset)
         st1 = dble(varread)*scale_factor + add_offset
         call readvar(ncid, xdim, ydim, 'stl2', start, count, varread, scale_factor, add_offset)
         st2 = dble(varread)*scale_factor + add_offset
         call readvar(ncid, xdim, ydim, 'stl3', start, count, varread, scale_factor, add_offset)
         st3 = dble(varread)*scale_factor + add_offset
         call readvar(ncid, xdim, ydim, 'stl4', start, count, varread, scale_factor, add_offset)
         st4 = dble(varread)*scale_factor + add_offset

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!if(numtasks.gt.1.and.istart.eq.0) call MPI_waitall(numtasks-2,req,stats,ierr)
         if (numtasks .gt. 1 .and. istart .eq. 0) call MPI_wait(request, status, ierr)
         istart = 0

!298:587 is equivalent to 135:424

!now send to nodes
         varpack(1:250, 1:290, 0:23, 1) = wind(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 2) = tempera(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 3) = press(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 4) = rh(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 5) = rain(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 6) = swdown(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 7) = rnetera(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 8) = st1(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 9) = st2(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 10) = st3(1067:1316, 298:587, 0:23)
         varpack(1:250, 1:290, 0:23, 11) = st4(1067:1316, 298:587, 0:23)

!write(6,*)'mirar',varpack(100,100,0,1),varpack(100,100,0,2),varpack(100,100,0,3),varpack(100,100,0,4),varpack(100,100,0,5),varpack(100,100,0,6),varpack(100,100,0,7)

!      call MPI_bcast(varpack(1,1,0,1),1,domblocke2obs,0,MPI_COMM_WORLD,ierr)

!        do n=1,numtasks-2
!           call MPI_isend(varpack(1,1,0,1),1,domblocke2obs,n,548,MPI_COMM_WORLD,req(n),ierr)
!        enddo

         deallocate (varread, press, rh, tempera, wind, rain, swdown, rnetera, st1, st2, st3, st4)

      END IF

      if (numtasks .gt. 1) call MPI_ibcast(varpack, 11, domblocke2obs, 0, MPI_COMM_WORLD, request, ierr)

      if (pid .ne. 0) call MPI_wait(request, status, ierr)

      end if

      IF ((pid .eq. 0 .and. numtasks .eq. 1) .or. (pid .gt. 0 .and. pid .lt. numtasks - 1)) then

!if(hour.eq.0.)then

!    if(numtasks.gt.1)then
!        call MPI_irecv(varpack(1,1,0,1),1,domblocke2obs,0,548,MPI_COMM_WORLD,request,ierr)
!        call  MPI_wait(request,status,ierr)
!    endif

!endif

!now interpolate

         dy = 0.25
         dx = 0.25
!    xswlat=-90.
!    xswlon=0.
         xswlat = -56.5
         xswlon = 266.5

         do j = js, je

!          gry = (lats(i,j) - xswlat) / dy + 1.
!          gry = (-lats(1,j) - xswlat) / dy + 1. !the grid as read from ERA5 has latitudes reversed
!          gry = ydim - 298 + 1 - (lats(1,j) - xswlat) / dy
            gry = dimeray + (xswlat - lats(1, j))/dy

            do i = 1, n2

               if (landmask(i, j) .gt. 0) then

                  if (lons(i, j) .lt. 0.) xlon = lons(i, j) + 360.
                  grx = (xlon - xswlon)/dx + 1.

                  call weights(grx, gry, ii, jj, wt1, wt2, wt3, wt4)

!wind
                  call gdtost3(varpack(1, 1, hour, 1), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, vels(i, j))
!temp
                  call gdtost3(varpack(1, 1, hour, 2), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, tempint)
!pres
                  call gdtost3(varpack(1, 1, hour, 3), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, pres(i, j))
!rh
                  call gdtost3(varpack(1, 1, hour, 4), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, relh)

                  vels(i, j) = max(vels(i, j), 0.)
                  temp(i, j) = tempint - 0.0065*(topo(i, j) - topoera5(i, j))

                  pres(i, j) = pres(i, j)*(temp(i, j)/tempint)**(9.81/(287.*0.0065))

                  relh = min(max(relh, 0.), 1.)

                  if (pid .gt. 0 .or. numtasks .eq. 1) then
                     qair(i, j) = relh*rslf(pres(i, j), temp(i, j))
                  end if

               end if

            end do
         end do

      END IF

      call MPI_TYPE_FREE(domblocke2obs, ierr)

   end subroutine readforcings

!********************************************************************************

   SUBROUTINE READFORCINGSACC(n2, js, je, filename, irec, hour, lats, lons &
                              , pcpgl, rshort, netrad, landmask, varpack, req, istart)
      use netcdf
      use interp_lib

      integer :: n2, js, je, i, j, ii, jj, n, istart, hour
      integer :: domblocke2obs, request, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2)
      character*50 :: filename
      character*200 :: filee2obs
      integer, dimension(n2, js:je) :: landmask
      real, dimension(n2, js:je) :: lats, lons, pcpgl, rshort, netrad
      real, allocatable, dimension(:, :, :) :: rain, swdown, rnetera
      integer*2, allocatable, dimension(:, :, :) :: varread
      real, dimension(dimerax, dimeray, 0:23, 11) :: varpack
      integer :: irec
      double precision :: scale_factor, add_offset
      real :: dx, dy, grx, gry, xswlat, xswlon, xlon, wt1, wt2, wt3, wt4

      integer :: ncid, varid, statusnc
      integer :: start(3), count(3)

      IF ((pid .eq. 0 .and. numtasks .eq. 1) .or. (pid .gt. 0 .and. pid .lt. numtasks - 1)) then

!now interpolate

         dy = 0.25
         dx = 0.25
!    xswlat=-90
!    xswlon=0.
         xswlat = -56.5
         xswlon = 266.5

         do j = js, je

!          gry = (lats(i,j) - xswlat) / dy + 1.
!          gry = (-lats(1,j) - xswlat) / dy + 1. !the grid as read from ERA5 has latitudes reversed
            gry = dimeray + (xswlat - lats(1, j))/dy

            do i = 1, n2

               if (landmask(i, j) .gt. 0) then

                  if (lons(i, j) .lt. 0.) xlon = lons(i, j) + 360.
                  grx = (xlon - xswlon)/dx + 1.

                  call weights(grx, gry, ii, jj, wt1, wt2, wt3, wt4)

!rain
                  call gdtost3(varpack(1, 1, hour, 5), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, pcpgl(i, j))
!rshort
                  call gdtost3(varpack(1, 1, hour, 6), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, rshort(i, j))
!rnet
                  call gdtost3(varpack(1, 1, hour, 7), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, netrad(i, j))

!convert from J/m2 to W/m2
                  rshort(i, j) = rshort(i, j)/3600.
                  netrad(i, j) = netrad(i, j)/3600.
               end if
            end do
         end do

      END IF

   end subroutine readforcingsacc

!********************************************************************************

   SUBROUTINE READFORCINGSRNET(n2, js, je, filename, irec, hour, lats, lons, rnet, landmask, netrad, req, istart)
      use netcdf
      use interp_lib

      integer :: n2, js, je, i, j, n, hour, istart
      integer ::domblocke2obs, request, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2)
      character*50 :: filename
      character*200 :: filee2obs
      integer, dimension(n2, js:je) :: landmask
      real, dimension(n2, js:je) :: lats, lons, rnet
      integer*2, dimension(dimerax, dimeray) :: varnc, varpastnc
      real, dimension(dimerax, dimeray) :: varpack, vartemp, netrad
      integer :: irec
      real :: dx, dy, grx, gry, xswlat, xswlon
      double precision :: scale_factor, add_offset

      integer :: ncid, varid, statusnc
      integer :: start(3), count(3)

      call MPI_Type_CONTIGUOUS(dimerax*dimeray, MPI_REAL, domblocke2obs, ierr)
      call MPI_Type_commit(domblocke2obs, ierr)

      IF (pid .eq. 0) THEN

         count(1) = dimerax
         count(2) = dimeray
         count(3) = 1
         start(1) = 1
         start(2) = 1
         start(3) = irec

         varpastnc = 0.

!Net short wave radiation

         filee2obs = '/mnt/netapp2/Store_uni/home/usc/fm/coc/Datos_GMM/Surface-Net-Solar-Radiation/'//filename(1:len_trim(filename))

!write(6,*)'reading from ',filee2obs,irec

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_inq_varid(ncid, 'ssr', varid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_var(ncid, varid, varnc, start, count)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'add_offset', add_offset)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         if (hour .ne. 3 .and. hour .ne. 15) then
            start(3) = irec - 1
            statusnc = nf90_get_var(ncid, varid, varpastnc, start, count)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)
            varpack = dble(varnc - varpastnc)*scale_factor !489.717729998627d0
         else
            varpack = dble(varnc)*scale_factor + add_offset !489.717729998627d0 + 16046091.141135d0
         end if

         varpastnc = 0
         start(3) = irec

!write(6,*)'read net SW'

!Net long wave radiation

         statusnc = nf90_inq_varid(ncid, 'str', varid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_var(ncid, varid, varnc, start, count)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'add_offset', add_offset)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         if (hour .ne. 3 .and. hour .ne. 15) then
            start(3) = irec - 1
            statusnc = nf90_get_var(ncid, varid, varpastnc, start, count)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)
            varpack = varpack + dble(varnc - varpastnc)*scale_factor !237.36992049807d0
         else
            varpack = varpack + dble(varnc)*scale_factor + add_offset !237.36992049807d0 - 4304363.18496025d0
         end if

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         varpastnc = 0

         varpack = varpack/(3.*3600.) !units from Jm-2 to Wm-2

         do j = 1, dimeray
            vartemp(:, dimeray - j + 1) = varpack(:, j)
         end do
         varpack = vartemp

!write(6,*)'read LW'

!now send to nodes
         if (numtasks .gt. 1 .and. istart .eq. 0) call MPI_waitall(numtasks - 2, req, stats, ierr)
         istart = 0

         netrad = varpack

         do n = 1, numtasks - 2
            call MPI_isend(netrad(1, 1), 1, domblocke2obs, n, 558, MPI_COMM_WORLD, req(n), ierr)
         end do

      END IF

      IF ((pid .eq. 0 .and. numtasks .eq. 1) .or. (pid .gt. 0 .and. pid .lt. numtasks - 1)) then

         if (numtasks .gt. 1) then
            call MPI_irecv(varpack(1, 1), 1, domblocke2obs, 0, 558, MPI_COMM_WORLD, request, ierr)
            call MPI_wait(request, status, ierr)
         end if

!now interpolate

         dy = 0.75
         dx = 0.75
         xswlat = -90.
         xswlon = 0.

         do j = js, je
            do i = 1, n2
               if (landmask(i, j) .gt. 0) then
                  gry = (lats(i, j) - xswlat)/dy + 1.
                  if (lons(i, j) .lt. 0.) then
                     grx = (lons(i, j) + 360.-xswlon)/dx + 1.
                  else
                     grx = (lons(i, j) - xswlon)/dx + 1.
                  end if
!netrad
                  call gdtost2(varpack(1, 1), dimerax, dimeray, grx, gry, rnet(i, j))
               end if
            end do
         end do

      END IF

      call MPI_TYPE_FREE(domblocke2obs, ierr)

   end subroutine readforcingsrnet
!********************************************************************************

   SUBROUTINE READFORCINGSSNOW(n2, js, je, filename, irec, hour, lats, lons, snow, landmask, snowera, request, istart)
      use netcdf
      use interp_lib

      integer :: n2, js, je, i, j, n, hour, istart, icount
      integer ::domblocke2obs, request, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2)
      character*50 :: filename
      character*200 :: filee2obs
      integer, dimension(n2, js:je) :: landmask
      real, dimension(n2, js:je) :: lats, lons, snow
      real, dimension(dimeralandx, dimeralandy, 8) :: snowera
      integer*2, allocatable, dimension(:, :, :) :: varread
      integer :: irec
      real :: dx, dy, grx, gry, xswlat, xswlon
      double precision :: scale_factor, add_offset

      integer :: ncid, varid, statusnc
      integer :: start(3), count(3)

      call MPI_Type_CONTIGUOUS(dimeralandx*dimeralandy*8, MPI_REAL, domblocke2obs, ierr)
      call MPI_Type_commit(domblocke2obs, ierr)

      if (hour .eq. 0) then
         IF (pid .eq. 0) THEN

            allocate (varread(xdimeraland, ydimeraland, 8))

            count(1) = xdimeraland
            count(2) = ydimeraland
            count(3) = 8
            start(1) = 1
            start(2) = 1
            start(3) = irec

!Net short wave radiation

            filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/SNOWMELT/ERA5LAND_snowmelt_'//filename(1:len_trim(filename))

!write(6,*)'reading from ',filee2obs,irec

            statusnc = nf90_open(filee2obs, 0, ncid)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

            statusnc = nf90_inq_varid(ncid, 'smlt', varid)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

            statusnc = nf90_get_var(ncid, varid, varread, start, count)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

            statusnc = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

            statusnc = nf90_get_att(ncid, varid, 'add_offset', add_offset)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

            statusnc = nf90_close(ncid)
            if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!      do j=1,ydimeraland
!         snowera(:,ydimeraland-j+1) = 1000.*max( dble( varread(:,j) )*scale_factor + add_offset , 0.)!in mm
!      enddo

!now send to nodes
!if(numtasks.gt.1.and.istart.eq.0) call MPI_waitall(numtasks-2,req,stats,ierr)

            if (numtasks .gt. 1 .and. istart .eq. 0) call MPI_wait(request, status, ierr)
            istart = 0

!747:1466 is equivalent to 336:1055

            snowera(1:620, 1:720, 1:8) = 1000*max(dble(varread(2666:3285, 747:1466, 1:8))*scale_factor + add_offset, 0.) !in mm

!        do n=1,numtasks-2
!           call MPI_isend(snowera(1,1,1),1,domblocke2obs,n,358,MPI_COMM_WORLD,req(n),ierr)
!        enddo

            deallocate (varread)

         END IF

         if (numtasks .gt. 1) call MPI_ibcast(snowera, 1, domblocke2obs, 0, MPI_COMM_WORLD, request, ierr)

         if (pid .ne. 0) call MPI_wait(request, status, ierr)

      end if

      IF ((pid .eq. 0 .and. numtasks .eq. 1) .or. (pid .gt. 0 .and. pid .lt. numtasks - 1)) then

!if(hour.eq.0)then
!    if(numtasks.gt.1)then
!        call MPI_irecv(snowera(1,1,1),1,domblocke2obs,0,358,MPI_COMM_WORLD,request,ierr)
!        call  MPI_wait(request,status,ierr)
!    endif
!endif

!now interpolate

         icount = hour/3 + 1

         dy = 0.1
         dx = 0.1
!    xswlat=-90
!    xswlon=0
         xswlat = -56.5
         xswlon = 266.5

         do j = js, je
!          gry = ydimeraland - 747 + 1 - (lats(1,j) - xswlat) / dy
            gry = dimeralandy + (xswlat - lats(1, j))/dy
            do i = 1, n2
               if (landmask(i, j) .gt. 0) then
!          gry = (lats(i,j) - xswlat) / dy + 1.
!          gry = (-lats(i,j) -xswlat ) / dy + 1. !the grid as read from ERA5 has latitudes reversed
!          gry = ydimeraland - 747 + 1 - (lats(1,j) - xswlat) / dy
                  if (lons(i, j) .lt. 0.) then
                     grx = (lons(i, j) + 360.-xswlon)/dx + 1.
                  else
                     grx = (lons(i, j) - xswlon)/dx + 1.
                  end if
!netrad
                  call gdtost2(snowera(1, 1, icount), dimeralandx, dimeralandy, grx, gry, snow(i, j))
               end if
            end do
         end do

      END IF

      call MPI_TYPE_FREE(domblocke2obs, ierr)

   end subroutine readforcingssnow

!********************************************************************************

   SUBROUTINE READFORCINGSSOILT(n2, js, je, nzg, filename, irec, hour, lats, lons, icefactor, landmask, topo, topoera &
                                , smoipack, req, istart)
      use netcdf
      use interp_lib

      integer :: n2, js, je, i, j, ii, jj, n, nlev, nzg, istart, hour
      integer ::domblocke2obs, request, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2)
      character*50 :: filename
      character*200 :: filee2obs
      character*6 :: directory
      character*15 :: prefix
      character*4 :: varname
      real :: tempint, elev, soilt, topofactor
      integer, dimension(n2, js:je) :: landmask
      real, dimension(n2, js:je) :: lats, lons, topo, topoera
      integer*1, dimension(n2, js:je, nzg - 14:nzg) :: icefactor
      real, allocatable, dimension(:, :, :) :: rain, swdown, rnetera
      integer*2, allocatable, dimension(:, :, :) :: varread
      real, dimension(dimerax, dimeray, 0:23, 11) :: smoipack
      integer :: irec
      real :: dx, dy, grx, gry, xswlat, xswlon, xlon, wt1, wt2, wt3, wt4
      double precision :: scale_factor, add_offset

      integer :: ncid, varid, statusnc
      integer :: start(3), count(3)

      IF ((pid .eq. 0 .and. numtasks .eq. 1) .or. (pid .gt. 0 .and. pid .lt. numtasks - 1)) then

!now interpolate

         dy = 0.25
         dx = 0.25
!    xswlat=-90.
!    xswlon=0.
         xswlat = -56.5
         xswlon = 266.5

         icefactor = 0

         do j = js, je

!          gry = (lats(i,j) - xswlat) / dy + 1.
!          gry = (-lats(1,j) - xswlat) / dy + 1. !the grid as read from ERA5 has latitudes reversed
            gry = dimeray + (xswlat - lats(1, j))/dy

            do i = 1, n2

               if (landmask(i, j) .gt. 0) then

                  if (lons(i, j) .lt. 0.) xlon = lons(i, j) + 360.
                  grx = (xlon - xswlon)/dx + 1.

                  topofactor = -0.0065*(topo(i, j) - topoera(i, j))

                  call weights(grx, gry, ii, jj, wt1, wt2, wt3, wt4)

                  do nlev = 1, 4
                     call gdtost3(smoipack(1, 1, hour, nlev + 7), dimerax, dimeray, ii, jj, wt1, wt2, wt3, wt4, tempint)

                     soilt = tempint + topofactor
                     select case (nlev)
                     case (1)
                        if (soilt .le. 273.15) icefactor(i, j, nzg) = 1
                     case (2)
                        if (soilt .le. 273.15) icefactor(i, j, nzg - 2:nzg - 1) = 1
                     case (3)
                        if (soilt .le. 273.15) icefactor(i, j, nzg - 7:nzg - 3) = 1
                     case (4)
                        if (soilt .le. 273.15) icefactor(i, j, nzg - 14:nzg - 8) = 1
                     end select

                  end do

               end if

            end do
         end do

      END IF

   end subroutine readforcingssoilt

!******************************************************************************************************
   subroutine readvar(ncid, xeraread, yeraread, varname, start, count, varread, scale_factor, add_offset)
      use netcdf

      integer :: xeraread, yeraread
      integer*2, dimension(xeraread, yeraread, 0:23) :: varread
      character(len=*) :: varname
      integer :: ncid, varid, statusnc
      double precision :: scale_factor, add_offset
      integer :: start(3), count(3)

      statusnc = nf90_inq_varid(ncid, varname(1:len_trim(varname)), varid)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_var(ncid, varid, varread, start, count)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

      statusnc = nf90_get_att(ncid, varid, 'add_offset', add_offset)
      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

   end subroutine readvar
!******************************************************************************************************
   subroutine READTOPOERA5(n2, n3, js, je, landmask, lats, lons, topoera5)
      use netcdf
      use interp_lib

      integer :: n2, n3, js, je, i, j, irec, iun, k, n
      real, dimension(n2, js:je) :: lats, lons, topoera5
      character*200 :: filee2obs
      integer, dimension(n2, js:je) :: landmask

      integer :: domblocke2obs, tag, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2), request

      integer*2, dimension(xdim, ydim) :: varread
      real, dimension(xdim, ydim) :: topoera

      real :: dx, dy, grx, gry, xswlat, xswlon, xlon
      double precision :: scale_factor, add_offset

      integer :: ncid, varid, statusnc
!integer :: start(3),count(3)

      call MPI_Type_CONTIGUOUS(xdim*ydim, MPI_REAL, domblocke2obs, ierr)
      call MPI_Type_commit(domblocke2obs, ierr)

      if (pid .eq. 0) then

         write (6, *) 'sending topera to nodes'

         filee2obs = '/mnt/netapp2/Store_uscfmgmm/ERA5/ERA5_static.nc4'

         statusnc = nf90_open(filee2obs, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_inq_varid(ncid, 'z', varid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_var(ncid, varid, varread)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'scale_factor', scale_factor)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'add_offset', add_offset)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         do j = 1, ydim
            topoera(:, ydim - j + 1) = dble(varread(:, j))*scale_factor + add_offset
         end do

         topoera = topoera/9.81

!        open(56,file='testera.dat' &
!           ,form='unformatted',convert='big_endian',access='direct',recl=4*xdim*ydim)

!                          write(56,rec=1) ((topoera(i,j),i=1,xdim),j=1,ydim)
!        close(56)

         do n = 1, numtasks - 2
            call MPI_isend(topoera(1, 1), 1, domblocke2obs, n, 8883, MPI_COMM_WORLD, req(n), ierr)
         end do

         if (numtasks .eq. 1) then

!now interpolate

            dy = 0.25
            dx = 0.25
            xswlat = -90
            xswlon = 0.

            do j = js, je
               do i = 1, n2
                  if (landmask(i, j) .gt. 0) then

                     if (lons(i, j) .lt. 0.) xlon = lons(i, j) + 360.
                     gry = (lats(i, j) - xswlat)/dy + 1.
                     grx = (xlon - xswlon)/dx + 1.
!elev
                     call gdtost2(topoera, xdim, ydim, grx, gry, topoera5(i, j))
                  end if
               end do
            end do

         else
            call MPI_waitall(numtasks - 2, req, stats, ierr)
         end if

      elseif (pid .lt. numtasks - 1) then
         call MPI_irecv(topoera(1, 1), 1, domblocke2obs, 0, 8883, MPI_COMM_WORLD, request, ierr)
         call MPI_wait(request, status, ierr)

!now interpolate

         dy = 0.25
         dx = 0.25
         xswlat = -90
         xswlon = 0.

         do j = js, je
            do i = 1, n2
               if (landmask(i, j) .gt. 0) then

                  if (lons(i, j) .lt. 0.) xlon = lons(i, j) + 360.
                  gry = (lats(i, j) - xswlat)/dy + 1.
                  grx = (xlon - xswlon)/dx + 1.

!elev
                  call gdtost2(topoera, xdim, ydim, grx, gry, topoera5(i, j))
               end if
            end do
         end do

      end if

      call MPI_TYPE_FREE(domblocke2obs, ierr)

   end subroutine readtopoera5

!******************************************************************************************************
   SUBROUTINE READLAI(n2, js, je, lats, lons, year, month, day, lai)
      use netcdf

      implicit none

      integer, parameter :: nxlai = 1200, nylai = 1200/9.81
      double precision, parameter :: radius = 6371007.181
      double precision, parameter :: pixel_size = 926.625433
      double precision, parameter :: d2r = 0.0174532925199
      double precision, parameter :: r2d = 57.2957795130823209
      double precision, parameter :: dx = 0.00833333, dy = 0.00833333

      integer :: n2, js, je, irec, year, month, day, currday, firstday
      integer, dimension(nxlai, nylai) :: modislai
      real, dimension(n2, js:je) :: lai, lats, lons
      real :: maxlat, minlat, maxlon, minlon, lat, lon
      integer :: i, j, ii, jj, minh, maxh, minv, maxv, hloc, vloc

      integer :: modid, varid, statusnc
      integer :: start(3), count(3)
      character(len=200) :: filemodis, modispath

      firstday = julday(month, 1, year)
      currday = julday(month, day, year)

!     irec = mod(day-1,4) + 1

      irec = (currday - firstday)/4 + 1
      start = (/1, 1, irec/)
      count = (/nxlai, nylai, 1/)

      write (modispath, '(a54,i4.4,a1,i4.4,i2.2,a1)') &
         '/home/usc/fm/coc/store/MODIS-LAI-NETCDF/South-America/', year, '/', year, month, '/'

!Tiles needed for this domain
      minlat = minval(lats(:, js:je))
      minlon = minval(lons(:, js:je)) - 360.
      maxlat = maxval(lats(:, js:je))
      maxlon = maxval(lons(:, js:je)) - 360.

!write(6,*)'esquinas',minlat,minlon,maxlat,maxlon,lats(1,js),lons(1,js),lats(1,je),lons(n2,js)

      minv = (nint(10800.-radius*maxlat*d2r/pixel_size) - 1)/nylai
      maxv = (nint(10800.-radius*minlat*d2r/pixel_size) - 1)/nylai

      minh = (nint(21600.+radius*minlon*d2r*cos(maxlat*d2r)/pixel_size) - 1)/nxlai
      maxh = (nint(21600.+radius*maxlon*d2r*cos(minlat*d2r)/pixel_size) - 1)/nxlai

!write(6,*)minh,maxh,minv,maxv

      do hloc = minh, maxh
         do vloc = minv, maxv

            if (month .eq. 12) then
               write (filemodis, '(a4,a1,i2.2a1,i2.2,a1,i4.4,i2.2,a3)') 'LAI-', 'h', hloc, 'v', vloc, '-', year + 1, 0, '.nc'
            else
               write (filemodis, '(a4,a1,i2.2a1,i2.2,a1,i4.4,i2.2,a3)') 'LAI-', 'h', hloc, 'v', vloc, '-', year, month, '.nc'
            end if
            write (6, *) 'reading LAI', filemodis(1:len_trim(filemodis)), irec, pid

            statusnc = nf90_open(modispath(1:len_trim(modispath))//filemodis(1:len_trim(filemodis)), 0, modid)
            if (statusnc .ne. nf90_noerr) then
               write (6, *) filemodis, 'is not there', statusnc !!call check(statusnc)
               cycle
            end if

            statusnc = nf90_inq_varid(modid, 'LAI', varid)
            if (statusnc .ne. nf90_noerr) call check(statusnc)

            statusnc = nf90_get_var(modid, varid, modislai, start, count)
            if (statusnc .eq. nf90_einvalcoords) then
               write (6, *) 'missing time in LAI file', filemodis(1:len_trim(filemodis))
               return
            end if
            if (statusnc .ne. nf90_noerr) call check(statusnc)

            statusnc = nf90_close(modid)
            if (statusnc .ne. nf90_noerr) call check(statusnc)

            where (modislai .gt. 100) modislai = 0

!now to process this tile

!           do jj=1,nylai
!                    lat = (10800. - float((vloc+1)*nylai-jj)) * pixel_size / radius
!                    j  = nint((lat*r2d - lats(1,js)) / dy )+ js
!              do ii=1,nxlai
!                    lon = (float(ii+hloc*nxlai) - 21600.) * pixel_size / (radius*cos(lat))
!                    i = nint((lon*r2d - lons(1,js)+360.) / dx )+ 1

!                       if(i.ge.1.and.i.le.n2.and.j.ge.js.and.j.le.je)then
!                            lai(i,j)=0.1*float(modislai(ii,jj))
!                       endif

!              enddo
!           enddo

            do j = js, je
               do i = 1, n2
                  jj = nint(radius*d2r*lats(i, j)/pixel_size) - 10800 + (vloc + 1)*nylai
                  ii = nint(d2r*(lons(i, j) - 360.)*radius*cos(lats(i, j)*d2r)/pixel_size) + 21600 - hloc*nxlai

                  if (ii .ge. 1 .and. ii .le. nxlai .and. jj .ge. 1 .and. jj .le. nylai) then
                     lai(i, j) = 0.1*float(modislai(ii, jj))
                  end if

               end do
            end do

         end do
      end do

   end subroutine readlai

!******************************************************************************************************
   subroutine READLAICLIM(n2, n3, js, je, month, lai)
      use netcdf

      integer :: n2, n3, js, je, i, j, irec, iun, k, n, month
      real, dimension(n2, js:je) :: lai

      integer :: tag, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2), request

      real, allocatable :: varread(:, :)
      integer*2, allocatable :: varreadbig(:, :)

      integer :: ncid, varid, statusnc
      integer :: start(3), count(3)
      double precision :: scalefactor, addoffset

      if (pid .eq. 0) then

         write (6, *) 'sending wtd to nodes'

         allocate (varreadbig(n2big, n3big))
         allocate (varread(n2, n3))

         start = (/1, 1, month/)
         count = (/n2big, n3big, 1/)

         statusnc = nf90_open('/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/lai_modis.nc', 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_inq_varid(ncid, 'LAI', varid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'scale_factor', scalefactor)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_att(ncid, varid, 'add_offset', addoffset)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_var(ncid, varid, varreadbig, start, count)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         varread(1:n2, 1:n3) = dble(varreadbig(nw:ne, ns:nn))*scalefactor + addoffset

!write(6,*)'mirar veg',varread(100,100),varreadbig(nw+100,n3big*3+100)

         do n = 1, numtasks - 2
            call MPI_isend(varread(1, nini(n)), 1, domblock(n), n, 53, MPI_COMM_WORLD, req(n), ierr)
         end do

         if (numtasks .eq. 1) then

            do j = js, je
               do i = 1, n2
                  lai(i, j) = varread(i, j)
               end do
            end do

         else
            call MPI_waitall(numtasks - 2, req, stats, ierr)
         end if

         deallocate (varread, varreadbig)

         close (iun)

      elseif (pid .lt. numtasks - 1) then
         call MPI_irecv(lai(1, js), 1, domblock(pid), 0, 53, MPI_COMM_WORLD, request, ierr)
         call MPI_wait(request, status, ierr)
      end if

   end subroutine readlaiclim
!******************************************************************************************************
   subroutine READLAICHINA(n2, n3, js, je, filelai, lai)
      use netcdf

      integer :: n2, n3, js, je, i, j, irec, iun, k, n, month
      real, dimension(n2, js:je) :: lai

      integer :: tag, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2), request

      real, allocatable :: varread(:, :)
      integer*1, allocatable :: varreadbig(:, :)

      integer :: ncid, varid, statusnc
      integer :: start(3), count(3)
      double precision :: scalefactor, addoffset
      character(len=*) :: filelai

      if (pid .eq. 0) then

         write (6, *) 'sending wtd to nodes'

         allocate (varreadbig(n2big, n3big))
         allocate (varread(n2, n3))

!     start=(/1,1,month/)
!     count=(/n2big,n3big,1/)

         statusnc = nf90_open(filelai, 0, ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_inq_varid(ncid, 'LAI', varid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!      statusnc = nf90_get_att(ncid,varid,'scale_factor',scalefactor)
!      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

!      statusnc = nf90_get_att(ncid,varid,'add_offset',addoffset)
!      if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_get_var(ncid, varid, varreadbig) !,start,count)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         statusnc = nf90_close(ncid)
         if (statusnc .ne. nf90_noerr) call handle_err(statusnc)

         varread(1:n2, 1:n3) = real(varreadbig(nw:ne, ns:nn))*0.1

!write(6,*)'mirar veg',varread(100,100),varreadbig(nw+100,n3big*3+100)

         do n = 1, numtasks - 2
            call MPI_isend(varread(1, nini(n)), 1, domblock(n), n, 53, MPI_COMM_WORLD, req(n), ierr)
         end do

         if (numtasks .eq. 1) then

            do j = js, je
               do i = 1, n2
                  lai(i, j) = varread(i, j)
               end do
            end do

         else
            call MPI_waitall(numtasks - 2, req, stats, ierr)
         end if

         deallocate (varread, varreadbig)

         close (iun)

      elseif (pid .lt. numtasks - 1) then
         call MPI_irecv(lai(1, js), 1, domblock(pid), 0, 53, MPI_COMM_WORLD, request, ierr)
         call MPI_wait(request, status, ierr)
      end if

   end subroutine readlaichina

!******************************************************************************************************

   integer function julday(imonth, iday, iyear)

      implicit none
      integer imonth, iday, iyear

! compute the julian day from a normal date

      julday = iday &
               + min(1, max(0, imonth - 1))*31 &
               + min(1, max(0, imonth - 2))*(28 + (1 - min(1, mod(iyear, 4)))) &
               + min(1, max(0, imonth - 3))*31 &
               + min(1, max(0, imonth - 4))*30 &
               + min(1, max(0, imonth - 5))*31 &
               + min(1, max(0, imonth - 6))*30 &
               + min(1, max(0, imonth - 7))*31 &
               + min(1, max(0, imonth - 8))*31 &
               + min(1, max(0, imonth - 9))*30 &
               + min(1, max(0, imonth - 10))*31 &
               + min(1, max(0, imonth - 11))*30 &
               + min(1, max(0, imonth - 12))*31

   end function julday

!******************************************************************************************************
   integer function daynumber(day, month, year)

      implicit none
      integer :: day, month, year, n

      daynumber = 0
      do n = 2003, year - 1
         daynumber = daynumber + 365 + (1 - min(1, mod(n, 4)))
      end do

      daynumber = daynumber + julday(month, day, year)

   end function daynumber
!     ******************************************************************

   subroutine handle_err(statusnc)
      use netcdf
      integer statusnc
      if (statusnc .ne. nf90_noerr) then
         print *, nf90_strerror(statusnc)
         stop 'Stopped'
      end if
   end subroutine handle_err

!********************************************************************************

   SUBROUTINE READFORCINGS6h(n2, js, je, lats, lons, topo, temp, press, ncount, FILE_NAME_NETCDF)
      use netcdf

      integer :: n2, js, je, i, j, ncount, n, ihour
      real :: grx, gry, elev, tempint

      integer :: request, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2), domblockera

      character(len=*) :: FILE_NAME_NETCDF
      character*80 :: VAR_NAME
      integer, parameter ::  NDIMS = 3
      integer :: start(NDIMS), count(NDIMS), ncid, varid
      integer*2, dimension(dimerax, dimeray) :: varread
      real, dimension(dimerax, dimeray) :: topo_era, temp_era, press_era
      real, dimension(dimerax, dimeray, 3) :: varpack
      real, dimension(n2, js:je) :: topo, lats, lons, temp, press

      call MPI_Type_CONTIGUOUS(3*dimerax*dimeray, MPI_REAL, domblockera, ierr)
      call MPI_Type_commit(domblockera, ierr)

      IF (pid .eq. 0) then

!read elevation

         VAR_NAME = 'z'

         call check(nf90_open('/mnt/netapp2/Store_uscfmgmm/ERAinterim/rootdepth_forcing/ERAInt_topography.nc', nf90_nowrite, ncid))
         call check(nf90_inq_varid(ncid, VAR_NAME, varid))
         count = (/dimerax, dimeray, 1/)
         start = (/1, 1, 1/)
         call check(nf90_get_var(ncid, varid, varread, start=start, count=count))
         call check(nf90_close(ncid))

         where (varread .eq. -32767) varread = 0

         topo_era = float(varread)*0.8265606640929 + 25942.8174325586
         write (6, *) 'mirar topo', topo_era(415, 132)

         varpack(:, :, 1) = topo_era(:, :)/9.81

!read temperature and calculate daily mean

         VAR_NAME = 't2m'

         call check(nf90_open(FILE_NAME_NETCDF, nf90_nowrite, ncid))
         call check(nf90_inq_varid(ncid, VAR_NAME, varid))
         count = (/dimerax, dimeray, 1/)

         temp_era = 0.
         n = ncount

         do ihour = 0, 18, 6

            n = n + 1
            start = (/1, 1, n/)
            call check(nf90_get_var(ncid, varid, varread, start=start, count=count))

            temp_era = temp_era + 0.25*(float(varread)*0.00197099929541104 + 259.350027683946)

         end do

         call check(nf90_close(ncid))

         write (6, *) 'mirar temp', temp_era(415, 132)

         varpack(:, :, 2) = temp_era(:, :)

!read pressures and calculate daily mean

         VAR_NAME = 'sp'

         call check(nf90_open(FILE_NAME_NETCDF, nf90_nowrite, ncid))
         call check(nf90_inq_varid(ncid, VAR_NAME, varid))
         count = (/dimerax, dimeray, 1/)

         press_era = 0.
         n = ncount

         do ihour = 0, 18, 6

            n = n + 1
            start = (/1, 1, n/)
            call check(nf90_get_var(ncid, varid, varread, start=start, count=count))

            press_era = press_era + 0.25*(float(varread)*0.846866273480537 + 78433.4203168633)

         end do

         call check(nf90_close(ncid))

         write (6, *) 'mirar press', press_era(415, 132)

         varpack(:, :, 3) = press_era(:, :)*1.e-2

!send to the nodes

         do n = 1, numtasks - 2
            call MPI_isend(varpack(1, 1, 1), 1, domblockera, n, 38, MPI_COMM_WORLD, req(n), ierr)
         end do

         if (numtasks .gt. 1) call MPI_waitall(numtasks - 2, req, stats, ierr)

      END IF

      IF ((pid .eq. 0 .and. numtasks .eq. 1) .or. (pid .gt. 0 .and. pid .lt. numtasks - 1)) then

         if (numtasks .gt. 1) then
            call MPI_irecv(varpack(1, 1, 1), 1, domblockera, 0, 38, MPI_COMM_WORLD, request, ierr)
            call MPI_wait(request, status, ierr)
         end if

!now interpolate to model grid

         do j = js, je
            do i = 1, n2

               gry = (lats(i, j) - eraswlat)/dxera + 1.
               grx = (lons(i, j) - eraswlon)/dxera + 1.

               call gdtost2(varpack(1, 1, 1), dimerax, dimeray, grx, gry, elev)
               call gdtost2(varpack(1, 1, 2), dimerax, dimeray, grx, gry, tempint)
               call gdtost2(varpack(1, 1, 3), dimerax, dimeray, grx, gry, press(i, j))

!adjust temp y press for elevation
               temp(i, j) = tempint - 0.0065*(topo(i, j) - elev)
               press(i, j) = press(i, j)*(temp(i, j)/tempint)**(9.81/(287.*0.0065))
            end do
         end do

      END IF

      call MPI_TYPE_FREE(domblockera, ierr)

   end subroutine readforcings6h

!********************************************************************************

   SUBROUTINE READFORCINGS3h(n2, js, je, lats, lons, rad, prep, ncount, FILE_NAME_NETCDF)
      use netcdf

      integer :: n2, js, je, i, j, ncount, n, ihour
      real :: grx, gry

      integer :: request, req(numtasks - 2), stats(MPI_STATUS_SIZE, numtasks - 2), domblockera

      character(len=*) :: FILE_NAME_NETCDF
      character*80 :: VAR_NAME
      integer, parameter ::  NDIMS = 3
      integer :: start(NDIMS), count(NDIMS), ncid, varid
      integer*2, dimension(dimerax, dimeray) :: varread
      real, dimension(dimerax, dimeray) :: rad_era, pp_era
      real, dimension(dimerax, dimeray, 2) :: varpack
      real, dimension(n2, js:je) :: lats, lons, rad, prep

      call MPI_Type_CONTIGUOUS(dimerax*dimeray*2, MPI_REAL, domblockera, ierr)
      call MPI_Type_commit(domblockera, ierr)

      IF (pid .eq. 0) then

!read shortwave radiation and calculate daily total

         VAR_NAME = 'ssr'

         call check(nf90_open(FILE_NAME_NETCDF, nf90_nowrite, ncid))
         call check(nf90_inq_varid(ncid, VAR_NAME, varid))
         count = (/dimerax, dimeray, 1/)

         rad_era = 0.
         n = ncount

         do ihour = 12, 24, 12

            n = n + 4
            start = (/1, 1, n/)
            call check(nf90_get_var(ncid, varid, varread, start=start, count=count))

            rad_era = rad_era + float(varread)*500.968245006333 + 16414725.5158775

         end do

         call check(nf90_close(ncid))

!read longwave radiation and calculate daily total

         VAR_NAME = 'str'

         call check(nf90_open(FILE_NAME_NETCDF, nf90_nowrite, ncid))
         call check(nf90_inq_varid(ncid, VAR_NAME, varid))
         count = (/dimerax, dimeray, 1/)

         n = ncount

         do ihour = 12, 24, 12

            n = n + 4
            start = (/1, 1, n/)
            call check(nf90_get_var(ncid, varid, varread, start=start, count=count))

            rad_era = rad_era + float(varread)*204.870248577053 - 4175573.43512429

         end do

         call check(nf90_close(ncid))

!change from J/m2 to W/m2

         varpack(:, :, 1) = rad_era(:, :)/(24.*3600.)

!now precipitation

         VAR_NAME = 'tp'

         call check(nf90_open(FILE_NAME_NETCDF, nf90_nowrite, ncid))
         call check(nf90_inq_varid(ncid, VAR_NAME, varid))
         count = (/dimerax, dimeray, 1/)

         pp_era = 0.
         n = ncount

         do ihour = 12, 24, 12

            n = n + 4
            start = (/1, 1, n/)
            call check(nf90_get_var(ncid, varid, varread, start=start, count=count))

            pp_era = pp_era + float(varread)*7.85746258683602e-06 + 0.257457619120269

         end do

         call check(nf90_close(ncid))

!change from m to mm

         varpack(:, :, 2) = pp_era(:, :)*1.e3

!    open(63,file= 'testforcingsera.dat'&
!      ,form='unformatted',convert='big_endian',access='direct',recl=4*dimerax*dimeray)
!          write(63,rec=1)((rad_era(i,j),i=1,dimerax),j=1,dimeray)
!          write(63,rec=2)((pp_era(i,j),i=1,dimerax),j=1,dimeray)
!    close(63)

!    stop

!send to the nodes

         do n = 1, numtasks - 2
            call MPI_isend(varpack(1, 1, 1), 1, domblockera, n, 39, MPI_COMM_WORLD, req(n), ierr)
         end do

         if (numtasks .gt. 1) call MPI_waitall(numtasks - 2, req, stats, ierr)

      END IF

      IF ((pid .eq. 0 .and. numtasks .eq. 1) .or. (pid .gt. 0 .and. pid .lt. numtasks - 1)) then

         if (numtasks .gt. 1) then
            call MPI_irecv(varpack(1, 1, 1), 1, domblockera, 0, 39, MPI_COMM_WORLD, request, ierr)
            call MPI_wait(request, status, ierr)
         end if

!now interpolate to model grid

         do j = js, je
            do i = 1, n2

               gry = (lats(i, j) - eraswlat)/dxera + 1.
               grx = (lons(i, j) - eraswlon)/dxera + 1.

               call gdtost2(varpack(1, 1, 1), dimerax, dimeray, grx, gry, rad(i, j))
               call gdtost2(varpack(1, 1, 2), dimerax, dimeray, grx, gry, prep(i, j))

            end do
         end do

!    open(63,file= 'testforcings.dat'&
!      ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*(je-js+1))
!          write(63,rec=1)((rad(i,j),i=1,n2),j=js,je)
!          write(63,rec=2)((prep(i,j),i=1,n2),j=js,je)
!    close(63)

!    stop

      END IF

      call MPI_TYPE_FREE(domblockera, ierr)

   end subroutine readforcings3h
!     ******************************************************************

   subroutine check(status)
      use netcdf
      integer, intent(in) :: status

      if (status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
   end subroutine check

!     ******************************************************************

   real function rslf(p, t)

!     This function calculates the liquid saturation vapor mixing ratio as
!     a function of pressure and Kelvin temperature

      implicit none
      real esl, x, t, p, c0, c1, c2, c3, c4, c5, c6, c7, c8
      parameter(c0=.6105851e+03, c1=.4440316e+02, c2=.1430341e+01)
      parameter(c3=.2641412e-01, c4=.2995057e-03, c5=.2031998e-05)
      parameter(c6=.6936113e-08, c7=.2564861e-11, c8=-.3704404e-13)

      x = max(-80., t - 273.16)

      esl = c0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + x*(c6 + x*(c7 + x*c8)))))))
      rslf = .622*esl/(p - esl)

   end function rslf

!     ******************************************************************

END MODULE MODULE_FORCINGS
