program driver

use module_parallel
use module_forcings
use module_rootdepth
use module_io
use module_wtable
use module_initial

implicit none

!integer, parameter :: n2=100,n3=100,nzg=40,restart=1,freedrain=0,riverswitch=1,nvar_out=15,writepar=0
integer, parameter :: n2=7320,n3=8520,nzg=40,restart=1,freedrain=0,riverswitch=1,nvar_out=15,writepar=1
!integer, parameter :: n2=500,n3=500,nzg=32,restart=0,freedrain=0,nvar_out=13 !change this when writing in parallel
!integer, parameter :: n2=1100,n3=300,nzg=32,restart=1,freedrain=0,nvar_out=14,writepar=0 !change this when writing in parallel
!integer, parameter :: n2=400,n3=400,nzg=32,restart=1,freedrain=0,nvar_out=14,writepar=0
!integer, parameter :: n2=100,n3=100,nzg=40,restart=0,freedrain=0,nvar_out=12,writepar=0
real, parameter :: deltat=1.*3600.,deltatwtd=24.*3600.,deltatriver=5.*60.,dxy=1./120.,steps=1.
integer, parameter :: maxinactivedays = 365 * 24 * nint(steps)

integer :: i,j,js,je,k,nmax,icount,n,nday,niter_river,iter_river

real, dimension(nzg+1) :: slz
real, dimension(nzg) :: dz
real :: inpair(2),outpair(2)
integer, allocatable, dimension(:,:,:) :: soiltxt,inactivedays
integer*1, allocatable, dimension(:,:,:) :: icefactor
integer, allocatable, dimension(:,:) :: landmask,fd,bfd
real, allocatable, dimension(:,:,:) :: smoi,watext,smoieq
real, allocatable, dimension(:,:) :: topo,topoera5,area,lats,lons,wind,temp,qair,netrad,rshort,rlon,press,precip,snow,snow_past,snow_fut,rain,lai &
                                    ,lai_past,lai_fut &
                                    ,press_past,press_fut,temp_past,temp_fut,qair_past,qair_fut,wind_past,wind_fut
real, allocatable, dimension(:,:) :: veg,hveg,wtd,smoiwtd,fdepth,rech,deeprech,qsrun,qsrunsum,qslat,qlat,qlatsum,qsprings &
                                     ,et_s,et_i,et_c,intercepstore,ppacum,waterdeficit,watextdeep,pppendepth
integer*1, allocatable, dimension(:,:) :: pppendepthold
real, allocatable, dimension(:,:) :: riverflow,qrf,qrfsum,delsfcwat,delsfcwatsum &
                  ,slope,riverdepth,riverwidth,riverlength,maxdepth,riverflowmean,floodheight,topoflood,riverarea,floodarea,riverchannel
real, allocatable, dimension(:,:,:) :: wtdflux,et_s_daily,et_c_daily,transptop,wtd_daily,smoi_daily
integer*1, allocatable, dimension(:,:,:) :: infilk
real, dimension(nzg,nstyp) :: fieldcp

real, allocatable, dimension(:,:,:) :: varoutput,varoutput2,varoutput3,varhis,varhis2,varhis3 
real, allocatable, dimension(:,:,:,:) :: varforcing,varforcing2,varforcing5
real, allocatable, dimension(:,:,:) :: varforcing4
real, allocatable, dimension(:,:) :: varforcing3

character*200 filename
integer :: hour,day,month,year,dayforc,monthforc,yearforc,lastday,daylai,monthlai,yearlai,jday,jdaylaipast,jdaylaifut,mdays,daypast,monthpast,yearpast,daysfromstart
integer :: daysforoutput(8)
integer :: hoursn,daysn,monthsn,yearsn
real :: tfact,dtlr,swlat,swlon
real*8 :: t1,t2,t3,t4,t5
integer :: mondays(12)
data mondays/31,28,31,30,31,30,31,31,30,31,30,31/
character*3 monthname(12)
data monthname/'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'/

character*100 date
integer :: irec,istart,request,requestsnow,req(nzg),req2(nzg),req3(nvar_out),rc,nsoil &
         ,istartforcing,istartforcing2,istartforcing3,istartforcing4,istartforcing5
integer :: istarthis,reqhis(nzg),reqhis2(0:nzg+1),reqhis3(3:7)
integer, allocatable, dimension(:) :: reqforcing,reqforcing2,reqforcing3,reqforcing4,reqforcing5

character (len = 300) :: filesoil,filetopo,filef,filewtd,fileveg,filehveg,filesmoieq,filerivers

swlat=-55.9958333333333
swlon=-92.9958333333333

filesoil='/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/sa_soil_modified.dat'
filef='/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/fdepth_sa.nc'
filetopo='/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/RIVERPARAMETERS/riverparameters_definitive_yama.dat'
filewtd='/mnt/netapp2/Store_uscfmgmm/GLOBALWTD/sciencepaper/ncfileshaibin_newrun/S_America_model_wtd_v2.nc'
fileveg='/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/vegsam_modis.nc'
filehveg='/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/hveg_sam.nc'
filesmoieq='/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/smoieq_SA_new.nc'
filerivers='/mnt/netapp2/Store_uscfmgmm/ROOTDEPTH/SHUWAL/PARALLEL/RIVERPARAMETERS/riverparameters_definitive_yama.dat'


!gmm intialize mpi stuff

   call MPI_INIT(ierr)
   if (ierr .ne. MPI_SUCCESS) then
      print *,'Error starting MPI program. Terminating.'
      call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
   end if

   call INITIALIZEDOMAIN(n2,n3,nzg,filetopo)

!Y dimension for each node

  NMAX=nend(pid)-nini(pid)+1
  js=nini(pid)
  je=nend(pid)

!gmm allocate variables

  allocate(landmask(n2,js:je))
  allocate(veg(n2,js:je))
  allocate(hveg(n2,js:je))
  allocate(soiltxt(2,n2,js:je))
  allocate(smoi(nzg,n2,js:je))
  allocate(watext(nzg,n2,js:je))
  allocate(smoieq(nzg,n2,js:je))

  allocate(wtd(n2,js:je))
  allocate(smoiwtd(n2,js:je))
  allocate(fdepth(n2,js:je))
  allocate(rech(n2,js:je))
  allocate(deeprech(n2,js:je))
  allocate(qsrun(n2,js:je))
  allocate(qsrunsum(n2,js:je))
  allocate(qslat(n2,js:je))
  allocate(qlat(n2,js:je))
  allocate(qlatsum(n2,js:je))
  allocate(qsprings(n2,js:je))
  allocate(et_s(n2,js:je))
  allocate(et_i(n2,js:je))
  allocate(et_c(n2,js:je))
  allocate(intercepstore(n2,js:je))
  allocate(ppacum(n2,js:je))
  allocate(waterdeficit(n2,js:je))
  allocate(watextdeep(n2,js:je))
  allocate(pppendepth(n2,js:je))
  allocate(pppendepthold(n2,js:je))
  allocate(inactivedays(0:nzg+1,n2,js:je))


  allocate(topo(n2,js:je))
  allocate(area(n2,js:je))
  allocate(lats(n2,js:je))
  allocate(lons(n2,js:je))

  allocate(riverflow(n2,js:je))
  allocate(qrf(n2,js:je))
  allocate(qrfsum(n2,js:je))
  allocate(delsfcwat(n2,js:je))
  allocate(delsfcwatsum(n2,js:je))
  allocate(slope(n2,js:je))
  allocate(riverdepth(n2,js:je))
  allocate(riverwidth(n2,js:je))
  allocate(riverlength(n2,js:je))
  allocate(maxdepth(n2,js:je))
  allocate(riverflowmean(n2,js:je))
  allocate(floodheight(n2,js:je))
  allocate(topoflood(n2,js:je))
  allocate(riverarea(n2,js:je))
  allocate(floodarea(n2,js:je))
  allocate(riverchannel(n2,js:je))
  allocate(fd(n2,js:je))
  allocate(bfd(n2,js:je))

  allocate(wind(n2,js:je))
  allocate(temp(n2,js:je))
  allocate(qair(n2,js:je))
  allocate(netrad(n2,js:je))
  allocate(rshort(n2,js:je))
  allocate(press(n2,js:je))
  allocate(rain(n2,js:je))
  allocate(snow(n2,js:je))
  allocate(snow_past(n2,js:je))
  allocate(snow_fut(n2,js:je))
  allocate(precip(n2,js:je))
  allocate(lai(n2,js:je))
  allocate(lai_past(n2,js:je))
  allocate(lai_fut(n2,js:je))
  allocate(press_past(n2,js:je))
  allocate(press_fut(n2,js:je))
  allocate(temp_past(n2,js:je))
  allocate(temp_fut(n2,js:je))
  allocate(qair_past(n2,js:je))
  allocate(qair_fut(n2,js:je))
  allocate(wind_past(n2,js:je))
  allocate(wind_fut(n2,js:je))
  allocate(icefactor(n2,js:je,nzg-14:nzg))
  allocate(topoera5(n2,js:je))

  allocate(wtdflux(n2,js:je,8))
  allocate(et_s_daily(n2,js:je,8))
  allocate(et_c_daily(n2,js:je,8))
  allocate(transptop(n2,js:je,8))
  allocate(infilk(n2,js:je,8))
  allocate(wtd_daily(n2,js:je,8))
  allocate(smoi_daily(n2,js:je,8))

if(writepar.eq.0)then
  allocate(varoutput(n2,js:je,nzg))
  allocate(varoutput2(n2,js:je,nzg))
  if(freedrain.eq.0)then
      allocate(varoutput3(n2,js:je,nvar_out))
  else
     allocate(varoutput3(n2,js:je,5))
  endif

  allocate(varhis(n2,js:je,nzg))
  allocate(varhis2(n2,js:je,0:nzg+1))
  allocate(varhis3(n2,js:je,3:7))
endif

allocate(reqforcing(numtasks-2))
allocate(varforcing(dimerax,dimeray,0:23,11))
!allocate(reqforcing2(numtasks-2))
!allocate(varforcing2(xdim,ydim,0:23,3))
!allocate(reqforcing3(numtasks-2))
!allocate(varforcing3(dimerax,dimeray))
allocate(reqforcing4(numtasks-2))
allocate(varforcing4(dimeralandx,dimeralandy,8))
!allocate(reqforcing5(numtasks-2))
!allocate(varforcing5(xdim,ydim,0:23,4))


!timestep for rivers
        niter_river = max(1,nint(deltat/deltatriver+0.4))
!        niter_river = 1
        dtlr = deltat / float(niter_river)

if(pid.eq.0)write(6,*)'niter_river',niter_river,dtlr


!counter to tell the output routine that it is first time

  istart=1
  istarthis=1
  istartforcing=1
  istartforcing2=1
  istartforcing3=1
  istartforcing4=1
  istartforcing5=1
  irec=0

!initialize some soil parameters

      call init_soil_param(fieldcp,nzg)

if(pid.eq.0)write(6,*)'slwilt',slwilt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!read in fixed fields

if(pid.eq.0)write(6,*)'reading soiltxt,topo and landmask'

       call READINITIAL(n2,n3,js,je,soiltxt,topo,fdepth,landmask,filesoil,filetopo,filef)

!soiltxt=13

!write(6,*)'mirar soiltxt',soiltxt(1,300,300),soiltxt(2,300,300)

       call READVEG(n2,n3,js,je,veg,fileveg)

       call READHVEG(n2,n3,js,je,hveg,filehveg)

       call READLATLON(n2,n3,js,je,lats,lons,area,dxy,swlat,swlon)

if(pid.eq.1)write(6,*)'SW corner of the domain',lats(1,1),lons(1,1)

       call READTOPOERA5(n2,n3,js,je,landmask,lats,lons,topoera5)

!initialize variables

   if(riverswitch.eq.1)then
         call READFLOWDIRECTION(n2,n3,js,je,fd,bfd,filerivers)
         call READRIVERPARAMETERS(n2,n3,js,je,riverlength,filerivers,2)
         call READRIVERPARAMETERS(n2,n3,js,je,riverwidth,filerivers,4)
         call READRIVERPARAMETERS(n2,n3,js,je,slope,filerivers,5)
         call READRIVERPARAMETERS(n2,n3,js,je,topoflood,filerivers,9)
         call READRIVERPARAMETERS(n2,n3,js,je,maxdepth,filerivers,11)
          riverarea = riverwidth*riverlength
          floodarea = max( area-riverarea , 0. )
          riverchannel = maxdepth*riverarea
   endif

!        open(56,file='testinput.dat' &
!            ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)

!                          write(56,rec=1) ((topoflood(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=2) ((float(soiltxt(1,i,j)),i=1,n2),j=1,n3)
!                          write(56,rec=3) ((float(soiltxt(2,i,j)),i=1,n2),j=1,n3)
!                          write(56,rec=4) ((fdepth(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=5) ((veg(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=6) ((hveg(i,j),i=1,n2),j=1,n3)
!           close(56)
!           write(6,*)'done with test'

!stop

if(pid.eq.0)write(6,*)'now to initialize soil depth and initial data'

       call INITIALIZESOILDEPTH(nzg,slz,dz)

       if(pid.eq.0)write(6,*)'soil layers',(slz(k),k=1,nzg+1),(dz(k),k=1,nzg)

       if(pid.eq.0)write(6,*)'soil nodes',(-0.5 * (slz(k) + slz(k+1)),k=1,nzg)

       pppendepthold=nzg+1
       wtdflux=0.
       et_s_daily=0.
       et_c_daily=0.
       transptop=0.
       infilk=nzg+1

       daypast = 1 !5day period number for output

       if(restart.eq.1)then
!                filename='/mnt/gluster/distributed/home/usc/fm/gmm/SAMERICA/outputSA/history_wt_01aug2013.nc'
!                filename='/mnt/netapp2/Store_uscfmgmm/GLUSTER/SAMERICA/outputSAinitial/history_wt_01jan2014.nc'
!               filename='/mnt/netapp2/Store_uscfmgmm/GLUSTER/SAMERICA/outputera5SAini/history_wt_01jan2014.nc'
                filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/outputera5SA/history_wt_01may2017.nc'
           call READHISTORYNC(n2,n3,js,je,nzg,smoi,smoiwtd,intercepstore,wtd,inactivedays,filename)
           if(freedrain.eq.0)then
!                 call READSMOIEQ(n2,n3,nzg,js,je,smoieq,filesmoieq)
                 call EQSOILMOISTUREtheor(n2,js,je,nzg,slz,dz,soiltxt,landmask,fdepth,smoieq)
           endif
           if(riverswitch.eq.1)then
                  call READHISTORYVARNC(n2,n3,js,je,riverflow,'RIVERFLOW',filename)
                  call READHISTORYVARNC(n2,n3,js,je,riverdepth,'RIVERDEPTH',filename)
                  call READHISTORYVARNC(n2,n3,js,je,floodheight,'FLOODHEIGHT',filename)
           endif

                  call READHISTORYBYTEVARNC(n2,n3,js,je,pppendepthold,'PPPENDEPTHOLD',filename)

                  call READHISTORYVARNC(n2,n3,js,je,wtdflux(1,js,1),'WTDFLUX',filename)
                  call READHISTORYVARNC(n2,n3,js,je,et_s_daily(1,js,1),'ET_S_DAILY',filename)
                  call READHISTORYVARNC(n2,n3,js,je,et_c_daily(1,js,1),'ET_C_DAILY',filename)
                  call READHISTORYVARNC(n2,n3,js,je,transptop(1,js,1),'TRANSPTOP',filename)
                  call READHISTORYBYTEVARNC(n2,n3,js,je,infilk(1,js,1),'INFILK',filename)

           where(veg.lt.0.5)landmask=0

!!!CHANGE THIS 
!           inactivedays=maxinactivedays+1
!           inactivedays(nzg+1,1:n2,js:je)=0
           
!           riverdepth=max(riverdepth,0.)
!           where(riverwidth.lt.1.00001)riverdepth=max(min(riverdepth,5.),floodheight)
!           where(riverwidth.lt.0.5)riverdepth=floodheight
!           where(riverwidth.lt.0.5)riverflow=0.

       else

           if(freedrain.eq.0)then
                 call READWTDNC(n2,n3,js,je,wtd,filewtd)
                 where(wtd.lt.-1.e5)wtd=0.
                 wtd=max(wtd,slz(1))
!                 call READSMOIEQ(n2,n3,nzg,js,je,smoieq,filesmoieq)
!                 call EQSOILMOISTURE(n2,js,je,nzg,slz,dz,deltat,soiltxt,landmask,smoieq)
                 call EQSOILMOISTUREtheor(n2,js,je,nzg,slz,dz,soiltxt,landmask,fdepth,smoieq)
           endif

           where(veg.lt.0.5)landmask=0
           where(landmask.eq.0)wtd=0.
           call INITIALIZE(n2,n3,js,je,nzg,freedrain,slz,dz,soiltxt,wtd,smoi,smoieq &
                          ,fdepth,topo,landmask,deltat,area)
           intercepstore = 0.

           inactivedays=maxinactivedays+1
           inactivedays(nzg+1,1:n2,js:je)=0

           floodheight = 0.

           if(riverswitch.eq.1)then
               call READRIVERPARAMETERS(n2,n3,js,je,riverflow,filerivers,6)
               call READRIVERPARAMETERS(n2,n3,js,je,riverdepth,filerivers,3)
           endif
       endif

if(pid.eq.1)write(6,*)'mirar smoi',smoi(5,5,5),landmask(5,5),smoieq(5,5,5)

       delsfcwat = 0.
       delsfcwatsum = 0.
       qrf = 0.
       qrfsum = 0.
       qsrun = 0.
       qsrunsum = 0.
       qsprings = 0.
       rech = 0.
       deeprech = 0.
       qslat = 0.
       qlat = 0.
       qlatsum = 0.
       et_s = 0.
       et_i = 0.
       et_c = 0.
       ppacum = 0.
       watext = 0.
       watextdeep = 0.
       waterdeficit = 0.
       pppendepth = 0.
       wind_past=0.
       wind_fut=0.
       temp_past=0.
       temp_fut=0.
       press_past=0.
       press_fut=0.
       qair_past=0.
       qair_fut=0.
       rain=0.
       snow=0.
       snow_past=0.
       snow_fut=0.
       rshort=0.
       netrad=0.

!initial time

hour=0
day=1
month=5
year=2017

       daysfromstart = daynumber(day,month,year)

daylai=1
monthlai=1



!read past temp,wind,qair,pres

!       icount = (day-1)*8 + hour/3+1
       icount = (day-1)*24 + hour +1

       write(filename,'(i4.4,a1,i2.2,a3)')year,'-',month,'.nc'

if(pid.eq.0)write(6,*)'reading now forcings first ',filename

       call READFORCINGS(n2,js,je,filename,icount,hour,topo,lats,lons,wind_past,temp_past,press_past,qair_past,landmask&
            ,topoera5,varforcing,request,istartforcing)

!       call READFORCINGSACC(n2,js,je,filename,icount,hour,lats,lons,rain,rshort,netrad,landmask &
!                        ,varforcing,reqforcing,istartforcing)

          daysn = day
          hoursn = hour
          yearsn = year
          monthsn = month

          icount = (daysn-1)*8 + hoursn/3 +1

          write(filename,'(i4.4,a1,i2.2,a3)')yearsn,'-',monthsn,'.nc'

          call READFORCINGSSNOW(n2,js,je,filename,icount,hoursn,lats,lons,snow_fut,landmask  &
                        ,varforcing4,requestsnow,istartforcing4)


!read lai. the first time of each file is always 1 jan

!       call READLAI(n2,js,je,lats,lons,year,monthlai,daylai,lai)

!read lai

        jday = julday(month,day,year)      

        jdaylaipast = jday - mod(jday-1,8)
        yearlai = year

            write(filename,'(i4.4,a13,i4.4,a1,i3.3,a3)')yearlai,'/SAhires_30s_',yearlai,'_',jdaylaipast,'.nc'
            filename='/mnt/netapp2/Store_uscfmgmm/LAI/'//filename(1:len_trim(filename))

            call READLAICHINA(n2,n3,js,je,filename,lai_past)

        jdaylaifut = jdaylaipast + 8 !this works because initial time is always the first day of a month
   

            write(filename,'(i4.4,a13,i4.4,a1,i3.3,a3)')yearlai,'/SAhires_30s_',yearlai,'_',jdaylaifut,'.nc'
            filename='/mnt/netapp2/Store_uscfmgmm/LAI/'//filename(1:len_trim(filename))

            call READLAICHINA(n2,n3,js,je,filename,lai_fut)



         tfact = float(jday -jdaylaipast) / float(jdaylaifut - jdaylaipast)
         lai = lai_past + (lai_fut - lai_past) * tfact

if(pid.eq.0)write(6,*)'read lai',jday,jdaylaipast,jdaylaifut

!write initial state
           if(restart.eq.0)then
                 if(freedrain.eq.0)then

                   write(date,'(i4,a1,i2.2,a1,i2.2,a9)')year,'-',month,'-',day,'_00:00:00'

                   if(writepar.eq.1)then

                   write(filename,'(a13,i2.2,a3,i4.4,a1,i3.3,a3)')'rootdaily_wt_',day,monthname(month),year,'_',pid,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/outputparera5/'//filename
                         call writeoutputnc_par(n2,n3,js,je,nzg,dz,smoi,waterdeficit,watext,watextdeep &
                                   ,wtd,smoiwtd,qsrunsum,rech,qsprings,qlatsum,et_s,et_i,et_c,ppacum,pppendepth &
                                   ,riverflow,qrfsum,delsfcwatsum &
                                   ,filename,irec,istart,req,req2,req3,date,nvar_out)
                   else

                   write(filename,'(a13,i2.2,a3,i4.4,a3)')'rootdaily_wt_',day,monthname(month),year,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/output/'//filename

                         call writeoutputnc(n2,n3,js,je,nzg,dz,smoi,waterdeficit,watext,watextdeep &
                                   ,wtd,smoiwtd,qsrunsum,rech,qsprings,qlatsum,et_s,et_i,et_c,ppacum,pppendepth &
                                   ,riverflow,qrfsum,delsfcwatsum &
                                   ,lai,filename,irec,istart,varoutput,varoutput2,varoutput3,req,req2,req3,date,nvar_out)
                   endif
                 else
!                         write(filename,'(a14,i2.2,a4)')'smoiwatext_fd_',month-1,'.dat'
!                         if(year.eq.1982)filename='smoiwatext_fd_12.dat'
                   write(filename,'(a13,i2.2,a3,i4.4,a6)')'rootdaily_fd_',day,monthname(month),year,'.grads'
                   filename='outputfd/'//filename

                   call writeoutputfd(n2,n3,js,je,nzg,smoi,waterdeficit,watext &
                                   ,qsrunsum,rech,et_c,ppacum &
                                   ,filename,irec,istart,req,req2,req3)
                 endif
           endif


!!!!!! do something if restarting on day day different than jan 1

DO WHILE(year.ne.2019)
!DO WHILE(month.ne.2)

!read icefactor before updating time
       icount = (day-1)*24 + hour +1

       write(filename,'(i4.4,a1,i2.2,a3)')year,'-',month,'.nc'

!       icount = (julday(month,day,year)-1)*8 + hour/3+1
!       write(filename,'(i4.4)')year
       call READFORCINGSSOILT(n2,js,je,nzg,filename,icount,hour,lats,lons,icefactor,landmask,topo,topoera5 &
                        ,varforcing,reqforcing,istartforcing)

!advance time to read forcing
      monthpast=month
      yearpast=year
      hour = hour + 1
      if(hour.ge.24)then
         hour=hour-24
         day=day+1
            lastday=mondays(month)
            if(month.eq.2.and.mod(year,4).eq.0)lastday=29
            if(day.eq.lastday+1)then
               day=1
               month=month+1
               if(month.eq.13)then
                    month=1
                    year=year+1
               endif
             endif
       daysfromstart = daynumber(day,month,year)
       endif


if(pid.eq.0)write(6,*)'Hour, day, month, year',hour,day,month,year,lastday,daysfromstart


       t1=mpi_wtime()

       icount = (day-1)*24 + hour +1
     
       write(filename,'(i4.4,a1,i2.2,a3)')year,'-',month,'.nc'

       call READFORCINGS(n2,js,je,filename,icount,hour,topo,lats,lons,wind_fut,temp_fut,press_fut,qair_fut,landmask &
                        ,topoera5,varforcing,request,istartforcing)

       press = 0.5 * ( press_past + press_fut )
       wind = 0.5 * ( wind_past + wind_fut )
       temp = 0.5 * ( temp_past + temp_fut )
       qair   = 0.5 * ( qair_past   + qair_fut   )

!move forcings to past
       press_past = press_fut
       wind_past = wind_fut
       temp_past = temp_fut 
       qair_past   = qair_fut


       call READFORCINGSACC(n2,js,je,filename,icount,hour,lats,lons,rain,rshort,netrad,landmask &
                        ,varforcing,reqforcing,istartforcing)

!read snowmelt for the next 3h

      hoursn = hour + 2 ! so, add 2 + 1 = 3h from present time

     if(mod(hoursn,3).eq.0)then

      daysn = day
      monthsn = month
      yearsn = year

      if(hoursn.ge.24)then
         hoursn=hoursn-24
         daysn=day+1
            lastday=mondays(month)
            if(month.eq.2.and.mod(year,4).eq.0)lastday=29
            if(daysn.eq.lastday+1)then
               daysn=1
               monthsn=month+1
               if(monthsn.eq.13)then
                    monthsn=1
                    yearsn=year+1
               endif
             endif
       endif

          icount = (daysn-1)*8 + hoursn/3 +1

          write(filename,'(i4.4,a1,i2.2,a3)')yearsn,'-',monthsn,'.nc'

          call READFORCINGSSNOW(n2,js,je,filename,icount,hoursn,lats,lons,snow_fut,landmask  &
                        ,varforcing4,requestsnow,istartforcing4)

          snow = ( snow_fut - snow_past ) * (deltat/(3.*3600.))

          if(hoursn.eq.0)then
              snow_past=0.
          else
              snow_past=snow_fut
          endif

       endif

       precip = rain + snow


!        open(56,file='testinput.dat' &
!           ,form='unformatted',convert='big_endian',access='direct',recl=4*n2*n3)


!                          write(56,rec=1) ((wind(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=2) ((temp(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=3) ((press(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=4) ((qair(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=5) ((rshort(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=6) ((netrad(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=7) ((rain(i,j),i=1,n2),j=1,n3)
!                          write(56,rec=8) ((topoera5(i,j),i=1,n2),j=1,n3)
!           close(56)
!stop
!           write(6,*)'done with test'

!write(6,*)press(100,100),wind(100,100),temp(100,100),qair(100,100),rshort(100,100),netrad(100,100),lai(100,100),veg(100,100)

       t2=mpi_wtime()

!now run the model
       if(freedrain.eq.0.and.hour.eq.nint(deltat/3600.)) then
              call LATERAL(n2,n3,js,je,soiltxt,wtd,qlat,fdepth,topo,landmask,deltatwtd,area,lats,dxy)
              qslat = 0.
              qlatsum = qlatsum + qlat*1.e3
       endif


       if(riverswitch.eq.1)then
              call GW2RIVER(n2,js,je,nzg,slz,deltat,soiltxt,landmask,wtd,maxdepth,riverdepth & 
                            ,riverwidth,riverlength,area,fdepth,qrf)
       endif

       t3=mpi_wtime()

       call ROOTDEPTH(freedrain,n2,js,je,nzg,slz,dz,deltat,landmask,veg,hveg,soiltxt,wind,temp,qair,press,netrad,rshort &
                      ,lai,precip,qsrun,smoi,smoieq,smoiwtd,wtd,waterdeficit,watext,watextdeep,rech,deeprech &
                      ,et_s,et_i,et_c,intercepstore,ppacum,pppendepth,pppendepthold &
                      ,qlat*deltat/deltatwtd,qslat,qsprings,inactivedays,maxinactivedays,fieldcp,fdepth,steps,floodheight &
                      ,qrf,delsfcwat,icefactor &
                      ,wtdflux(1,js,daypast),et_s_daily(1,js,daypast),et_c_daily(1,js,daypast),transptop(1,js,daypast),infilk(1,js,daypast))


!if(hour.eq.0)write(6,*)'mirar',daypast,wtdflux(50,50,daypast),et_s_daily(50,50,daypast),et_c_daily(50,50,daypast),transptop(50,50,daypast),infilk(50,50,daypast),wtd(50,50)

                       qrfsum = qrfsum + qrf*1.e3
                       qsrunsum=qsrunsum+qsrun*1.e3
                       delsfcwatsum=delsfcwatsum-delsfcwat*1.e3

                       if(mod(daysfromstart-1,5).eq.0.and.hour.eq.0)then
                                wtd_daily(:,:,daypast)=wtd(:,:)
                                smoi_daily(1:n2,js:je,daypast)=smoi(nzg,1:n2,js:je)
                                daysforoutput(daypast)=daysfromstart-1
                                daypast = daypast + 1 ! advance counter           
                       endif

       t4=mpi_wtime()

!this is not needed now, all layers are treaded the same       
!       if(freedrain.eq.0.and.hour.eq.3) call WTABLE(n2,n3,js,je,nzg,slz,dz,area,soiltxt,wtd,deeprech,rech,qslat,fdepth &
!                                     ,topo,landmask,deltatwtd,smoi,smoieq,smoiwtd,qsprings)

       if(riverswitch.eq.1)then

              call FLOODING(n2,js,je,deltat,fd,bfd,topoflood,area,riverwidth,riverlength,riverdepth,floodheight,delsfcwat)

          do iter_river=1,niter_river

              call RIVERS_KW_FLOOD(n2,js,je,deltat,dtlr,fd,bfd,riverflow,qsrun,qrf,delsfcwat &
                  ,slope,riverdepth,riverwidth,riverlength,maxdepth,area,riverarea,floodarea,riverchannel &
                  ,riverflowmean,floodheight,topoflood)

          enddo

              qsrun = 0.
              delsfcwat = 0.
        
       endif


       t5=mpi_wtime()

if(pid.eq.numtasks/2)write(6,'(a12,5f7.3)')'CPU(sec)',T5-T1,t2-t1,t3-t2,t4-t3,t5-t4


inpair(1)=T5-T1
inpair(2)=float(pid)

call MPI_REDUCE(inpair, outpair, 1, MPI_2REAL, MPI_MAXLOC, numtasks-1 , MPI_COMM_WORLD , ierr)

if(pid.eq.numtasks-1)write(6,'(a27,f7.3,f7.0)')'max CPU total step, pid',outpair(1),outpair(2)
!if(pid.eq.1)write(6,'(a12,5f7.3)')'pid 1 CPU(sec)',T5-T1,t2-t1,t3-t2,t4-t3,t5-t4




!output

!        if(year.ge.1981)then
!         if(day.eq.2)then

          if(day.eq.1.and.hour.eq.0)then

riverflowmean = riverflowmean/(float(lastday)*86400.)
delsfcwatsum = delsfcwatsum/float(lastday)
qrfsum = qrfsum/float(lastday)
qsrunsum=qsrunsum/float(lastday)
rech=rech/float(lastday)
qsprings=qsprings/float(lastday)
qlatsum=qlatsum/float(lastday)
et_s= et_s/float(lastday)
et_i=et_i/float(lastday)
et_c=et_c/float(lastday)
ppacum=ppacum/float(lastday)
waterdeficit=waterdeficit/float(lastday)
watext=watext/float(lastday)
watextdeep=watextdeep/float(lastday)

!              irec=0
!          endif
!              if(month.gt.1.or.year.eq.1983)then
            if(freedrain.eq.0)then

                write(date,'(i4,a1,i2.2,a1,i2.2,a9)')year,'-',month,'-',day,'_00:00:00'

                if(writepar.eq.1)then

                   write(filename,'(a13,i2.2,a3,i4.4,a1,i3.3,a3)')'rootdaily_wt_',day,monthname(month),year,'_',pid,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/outputparera5/'//filename

                         call writeoutputnc_par(n2,n3,js,je,nzg,dz,smoi,waterdeficit,watext,watextdeep &
                                   ,wtd,smoiwtd,qsrunsum,rech,qsprings,qlatsum,et_s,et_i,et_c,ppacum,pppendepth &
                                   ,riverflowmean,qrfsum,delsfcwatsum &
!                                   ,wind,temp,qair,press,netrad,rshort,pet,lai  &
                                   ,filename,irec,istart,req,req2,req3,date,nvar_out)

                   write(filename,'(a12,a3,i4.4,a1,i3.3,a3)')'dailyoutput_',monthname(monthpast),yearpast,'_',pid,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/outputparera5/'//filename

                         call writeoutputnc_daily_par(n2,n3,js,je,nzg,daypast-1 &
                                             ,wtdflux,et_s_daily,et_c_daily,transptop,infilk,smoi_daily,wtd_daily &
                                             ,daysforoutput,filename )


                   else

                   write(filename,'(a13,i2.2,a3,i4.4,a3)')'rootdaily_wt_',day,monthname(month),year,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/output/'//filename

                         call writeoutputnc(n2,n3,js,je,nzg,dz,smoi,waterdeficit,watext,watextdeep &
                                   ,wtd,smoiwtd,qsrunsum,rech,qsprings,qlatsum,et_s,et_i,et_c,ppacum,pppendepth &
                                   ,riverflowmean,qrfsum,delsfcwatsum &
                                   ,lai,filename,irec,istart,varoutput,varoutput2,varoutput3,req,req2,req3,date,nvar_out)


                   write(filename,'(a12,a3,i4.4,a3)')'dailyoutput_',monthname(monthpast),yearpast,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/output/'//filename

                         call writeoutputnc_daily_par(n2,n3,js,je,nzg,daypast-1 &
                                             ,wtdflux,et_s_daily,et_c_daily,transptop,infilk,smoi_daily,wtd_daily &
                                             ,daysforoutput,filename )


                 endif

              else
!                         write(filename,'(a14,i2.2,a4)')'smoiwatext_fd_',month-1,'.dat'
!                         if(year.eq.1982)filename='smoiwatext_fd_12.dat'
                   write(filename,'(a13,i2.2,a3,i4.4,a6)')'rootdaily_fd_',day,monthname(month),year,'.grads'
                   filename='outputfd/'//filename

                   call writeoutputfd(n2,n3,js,je,nzg,smoi,waterdeficit,watext &
                                   ,qsrunsum,rech,et_c,ppacum &
                                   ,filename,irec,istart,req,req2,req3)
                 endif
!               endif

               riverflowmean = 0.
               delsfcwatsum = 0.
               qrfsum = 0.
               qsrunsum = 0.
               rech = 0.
               qsprings = 0.
               qlatsum = 0.
               et_s = 0.
               et_i = 0.
               et_c = 0.
               ppacum = 0.
               waterdeficit = 0.
               watext = 0.
               watextdeep = 0.
               pppendepth = 0.

               wtdflux(:,:,1:daypast-1)=0.
               et_s_daily(:,:,1:daypast-1)=0.
               et_c_daily(:,:,1:daypast-1)=0.
               transptop(:,:,1:daypast-1)=0.
               infilk(:,:,1:daypast-1)=nzg+1



           endif
!        endif
!       endif


!history

          if(day.eq.1.and.hour.eq.0)then

                   write(date,'(i4,a1,i2.2,a1,i2.2,a9)')year,'-',month,'-',day,'_00:00:00'

             if(writepar.eq.1)then

                   write(filename,'(a11,i2.2,a3,i4.4,a1,i3.3,a3)')'history_wt_',day,monthname(month),year,'_',pid,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/outputparera5/'//filename
                   if(pid.eq.1)write(6,*)'writing history file',filename,date,daypast

                   call writehistorync_par(n2,n3,js,je,nzg,smoi,intercepstore,wtd,inactivedays &
                                          ,riverflow,riverdepth,floodheight &
                                          ,pppendepthold &
                                          ,wtdflux(1,js,daypast),et_s_daily(1,js,daypast),et_c_daily(1,js,daypast),transptop(1,js,daypast),infilk(1,js,daypast) &
!                                           ,riverflow,netrad,floodheight &
                                          ,filename,date)
              else

                   write(filename,'(a11,i2.2,a3,i4.4,a3)')'history_wt_',day,monthname(month),year,'.nc'
                   filename='/mnt/lustre/scratch/home/usc/fm/gmm/ROOTDEPTH/SAMERICA/output/'//filename
                   if(pid.eq.1)write(6,*)'writing history file',filename,date

                   call writehistorync(n2,n3,js,je,nzg,smoi,smoiwtd,intercepstore,wtd,inactivedays &
                                      ,riverflow,riverdepth,floodheight &
                                      ,filename,date,istarthis,varhis,varhis2,varhis3,reqhis,reqhis2,reqhis3)
          
               endif

               wtdflux(:,:,1)=wtdflux(:,:,daypast)
               et_s_daily(:,:,1)=et_s_daily(:,:,daypast)
               et_c_daily(:,:,1)=et_c_daily(:,:,daypast)
               transptop(:,:,1)=transptop(:,:,daypast)
               infilk(:,:,1)=infilk(:,:,daypast)

               wtdflux(:,:,daypast)=0.
               et_s_daily(:,:,daypast)=0.
               et_c_daily(:,:,daypast)=0.
               transptop(:,:,daypast)=0.
               infilk(:,:,daypast)=nzg+1

               daypast=1  !bring back counter to 1

          endif


!now LAI from MODIS
!       if(mod( julday(month,day,year) -1 , 4 ) .eq. 0 .and.hour.eq.0)then
!                 if(pid.eq.0)write(6,*)'time to read LAI',julday(month,day,year)
!                 call READLAI(n2,js,je,lats,lons,year,month,day,lai)
!       endif


  if(hour.eq.0)then

        jday = julday(month,day,year)
        if(jday.eq.1)jdaylaifut=1

        if(jday.eq.jdaylaifut)then

            jdaylaipast = jdaylaifut
            jdaylaifut = jday + 8
            yearlai = year

            lastday=365
            if(mod(year,4).eq.0)lastday=366

            if(jdaylaifut.gt.lastday)then
               jdaylaifut=lastday+1
               yearlai = year + 1
               write(filename,'(i4.4,a13,i4.4,a1,i3.3,a3)')yearlai,'/SAhires_30s_',yearlai,'_',1,'.nc'
            else
               write(filename,'(i4.4,a13,i4.4,a1,i3.3,a3)')yearlai,'/SAhires_30s_',yearlai,'_',jdaylaifut,'.nc'
            endif

            filename='/mnt/netapp2/Store_uscfmgmm/LAI/'//filename(1:len_trim(filename))

            lai_past = lai_fut
            call READLAICHINA(n2,n3,js,je,filename,lai_fut)

if(pid.eq.0)write(6,*)'jdaylaipast,jdaylaifut,yearlai',jdaylaipast,jdaylaifut,yearlai

          endif

         tfact = float(jday -jdaylaipast) / float(jdaylaifut - jdaylaipast)
         lai = lai_past + (lai_fut - lai_past) * tfact

   endif


ENDDO


   call MPI_FINALIZE(ierr)

end
