MODULE module_initial

implicit none

CONTAINS

!**********************************************************************************************
SUBROUTINE INITIALIZE(n2,n3,js,je,nzg,freedrain,slz,dz,soiltextures,wtd,smoi,smoieq &
                      ,fdepth,topo,landmask,deltat,area)
use module_parallel
use module_rootdepth
use module_wtable

integer :: n2,n3,js,je,nzg,nss,nee,i,j,iwtd,nsoil,nsoil1,k,freedrain
real , dimension(nzg+1) :: slz
real , dimension(nzg) :: dz,vctr4,vctr5,vctr6
integer, dimension(n2,js:je) :: landmask
real, dimension(n2,js:je) :: wtd,smoiwtd,fdepth,topo,area,qlat,klat
real, dimension(nzg,n2,js:je) :: smoi,smoieq
real :: deltat
real :: alpha,flux,wgpmid,wmid,k1,d0,d1,d2,z0,z1,z2,tol!,zbrent1,zbrent2,zbrent3
real :: smoisat,hydcon,psisat,smoicp,beta

integer, dimension(2,n2,js:je) :: soiltextures
integer, dimension(nzg) :: soiltxt



   nss=js+1
   if(pid.eq.1.or. numtasks.eq.1)nss=js

   nee=je-1
   if(pid.eq.numtasks-2.or.numtasks.eq.1)nee=je

         do k=1,nzg
            vctr4(k) = 0.5 * (slz(k) + slz(k+1))
         enddo
         do k = 2,nzg
            vctr5(k) = vctr4(k) - vctr4(k-1)
            vctr6(k) = 1. / vctr5(k)
         enddo


!Calculate lateral flow

IF(freedrain.eq.0.)then

qlat=0.

do j=js,je
   do i=1,n2
       nsoil=soiltextures(1,i,j)
       klat(i,j)=slcons(nsoil)*klatfactor(nsoil)
   enddo
enddo

call lateralflow(n2,n3,js,je,wtd,qlat,fdepth,topo,landmask,deltat,area,klat)

ENDIF

!now initialize soil moisture

smoi=0.
smoiwtd=0.


 DO j=nss,nee
   DO i=1,n2

where(slz.lt.-0.3)
     soiltxt=soiltextures(1,i,j)
elsewhere
     soiltxt=soiltextures(2,i,j)
endwhere

IF(landmask(i,j).gt.0)then

    if(freedrain.eq.1)then

                        do k=1,nzg
                               nsoil=soiltxt(k)
                               smoi(k,i,j)=0.5*(slmsts(nsoil)+soilcp(nsoil))
                        enddo
    else

      flux=min(qlat(i,j),0.)/deltat

!if(i.eq.51.and.j.eq.51)flux=0.

!flux=0.
      !check where the wtd is

      do k=1,nzg
        if(wtd(i,j)+1.e-6.lt.slz(k))exit
      enddo
      iwtd=k

!first below wtd 
              do k=1,iwtd-2
                nsoil = soiltxt(k)
                smoi(k,i,j)=slmsts(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)
!if(i.eq.51.and.j.eq.51)write(6,*)'mirar 1',k,smoi(k,i,j),fdepth(i,j),slz(k)
              enddo
!first wgp in the layer where the water table is
               k=iwtd-1
               nsoil = soiltxt(max(k,1))
               smoisat = slmsts(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)
               wgpmid = smoieq(k,i,j)
               smoi(k,i,j)=(wgpmid*(slz(k+1)-wtd(i,j))+smoisat*(wtd(i,j)-slz(k))) &
                                  / (slz(k+1)-slz(k))
!if(i.eq.51.and.j.eq.51)write(6,*)'mirar 2',k,smoi(k,i,j),smoisat,smoieq(k,i,j),wtd(i,j),flux,slz(k)
!then the rest 
            do k = iwtd,nzg
               nsoil = soiltxt(nzg)

               hydcon = slcons(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)
               smoisat = slmsts(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)
               psisat = slpots(nsoil)*min(max(exp(-(vctr4(k)+1.5)/fdepth(i,j)),1.),10.)
               smoicp = soilcp(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)

               z1 = vctr4(k) - vctr4(k-1)

               tol=0.0001
               wmid=zbrent(0.,1.,tol,slbs(nsoil),smoisat,hydcon,psisat,slz(k),vctr4(k),vctr6(k),z1,smoi(k-1,i,j),flux)

               smoi(k,i,j)=min(max(wmid,smoicp),smoisat)
!               smoi(k,i,j)=min(smoi(k,i,j),smoi(k-1,i,j))

!if(i.eq.51.and.j.eq.51)write(6,*)'mirar 3',k,smoi(k,i,j),smoisat,smoicp,wmid,flux

             enddo
!         endif


endif
ENDIF


   ENDDO
 ENDDO


end subroutine initialize
!     ******************************************************************

SUBROUTINE EQSOILMOISTUREtheor(n2,js,je,nzg,slz,dz,soiltextures,landmask,fdepth,smoieq)
use module_parallel
use module_rootdepth

integer :: n2,js,je,nzg,nss,nee,i,j,nsoil,nsoil1,k,iter
real , dimension(nzg+1) :: slz
real , dimension(nzg) :: dz,vctr2,vctr4,vctr5,vctr6
integer, dimension(n2,js:je) :: landmask
real, dimension(n2,js:je) :: fdepth
real, dimension(nzg,n2,js:je) :: smoieq
real :: wmid,flux,smoisat,hydcon,psisat,smoicp,z1,smoisatdw,tol

integer, dimension(2,n2,js:je) :: soiltextures
integer, dimension(0:nzg) :: soiltxt



   nss=js+1
   if(pid.eq.1.or. numtasks.eq.1)nss=js

   nee=je-1
   if(pid.eq.numtasks-2.or.numtasks.eq.1)nee=je

         do k=1,nzg
            vctr2(k) = 1. / dz(k)
            vctr4(k) = 0.5 * (slz(k) + slz(k+1))
         enddo
         do k = 2,nzg
            vctr5(k) = vctr4(k) - vctr4(k-1)
            vctr6(k) = 1. / vctr5(k)
         enddo
            vctr5(1) = dz(1)
            vctr6(1) = 1. / vctr5(1)

smoieq = 0.

 DO j=nss,nee
   DO i=1,n2

     if(landmask(i,j).eq.0)cycle

where(slz.lt.-0.3)
     soiltxt=soiltextures(1,i,j)
elsewhere
     soiltxt=soiltextures(2,i,j)
endwhere
     soiltxt(0)=soiltxt(1)

      flux=0.

    do k=1,nzg

               nsoil = soiltxt(max(k-1,1))
               smoisatdw=slmsts(nsoil)*max(min(exp((vctr4(max(k-1,1))+1.5)/fdepth(i,j)),1.),0.1)

               nsoil = soiltxt(k)

               hydcon = slcons(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)
               smoisat = slmsts(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)
               psisat = slpots(nsoil)*min(max(exp(-(vctr4(k)+1.5)/fdepth(i,j)),1.),10.)
               smoicp = soilcp(nsoil)*max(min(exp((vctr4(k)+1.5)/fdepth(i,j)),1.),0.1)


               tol=0.0001
               wmid=zbrent(0.,1.,tol,slbs(nsoil),smoisat,hydcon,psisat,slz(k),vctr4(k),vctr6(k),vctr5(k),smoisatdw,flux)

               smoieq(k,i,j)=min(max(wmid,smoicp),smoisat)

    enddo

   ENDDO
 ENDDO

end subroutine eqsoilmoisturetheor


!     ******************************************************************

SUBROUTINE EQSOILMOISTURE(n2,js,je,nzg,slz,dz,dtll,soiltextures,landmask,smoieq)
use module_parallel
use module_rootdepth

integer :: n2,js,je,nzg,nss,nee,i,j,nsoil,nsoil1,k,iter
real , dimension(nzg+1) :: slz
real , dimension(nzg) :: dz,vctr2,vctr4,vctr5,vctr6
integer, dimension(n2,js:je) :: landmask
real, dimension(nzg,n2,js:je) :: smoieq
real :: dtll,alpha,b,b1,ddw,dup,kfdw,kfup,CC,beta,DD,ff,dff,dx,smoimid,smoi,smoibotbc

integer, dimension(2,n2,js:je) :: soiltextures
integer, dimension(0:nzg) :: soiltxt



   nss=js+1
   if(pid.eq.1.or. numtasks.eq.1)nss=js

   nee=je-1
   if(pid.eq.numtasks-2.or.numtasks.eq.1)nee=je

         do k=1,nzg
            vctr2(k) = 1. / dz(k)
            vctr4(k) = 0.5 * (slz(k) + slz(k+1))
         enddo
         do k = 2,nzg
            vctr5(k) = vctr4(k) - vctr4(k-1)
            vctr6(k) = 1. / vctr5(k)
         enddo
            vctr5(1) = dz(1)
            vctr6(1) = 1. / vctr5(1)

smoieq = 0.

 DO j=nss,nee
   DO i=1,n2

     if(landmask(i,j).eq.0)cycle

where(slz.lt.-0.3)
     soiltxt=soiltextures(1,i,j)
elsewhere
     soiltxt=soiltextures(2,i,j)
endwhere
     soiltxt(0)=soiltxt(1)

    do k=1,nzg

     call SOILFLUXES_EQSMOI(i,j,nzg,k,dtll,slz,dz,soiltextures(1,i,j),smoibotbc)

     nsoil = soiltxt(k)
     b = slbs(nsoil)
     ddw = -slcons(nsoil)*slpots(nsoil)*b / slmsts(nsoil)**(b+3.)
     kfdw = slcons(nsoil) / slmsts(nsoil)**(2. * b + 3.)

     nsoil1 = soiltxt(k-1)
     b1 = slbs(nsoil1)
     dup = -slcons(nsoil1)*slpots(nsoil1)*b1 / slmsts(nsoil1)**(b1+3.)
     kfup = slcons(nsoil1) / slmsts(nsoil1)**(2. * b1 + 3.)


     CC = (slz(k)-vctr4(k))*vctr6(k)
     alpha = 1. + CC
     beta = -CC * slmsts(nsoil1)

    smoi =  slmsts(nsoil)

    do iter = 1, 100

     smoimid = alpha * smoi + beta
     DD =  0.5 * ( ddw * smoimid**(b+2.) + dup *smoimid**(b1+2.) )
!     ff =  DD * ( smoi-slmsts(nsoil1) ) * vctr6(k) &
     ff =  DD * ( smoi-smoibotbc ) * vctr6(k) &
          + 0.5 * ( kfdw * smoimid**(2. * b + 3.) +  kfup * smoimid**(2. * b1 + 3.) )

     dff = 0.5 * alpha * ( ddw * (b+2.) * smoimid**(b+1.) &
!                 + dup * (b1+2.) * smoimid**(b1+1.) ) * ( smoi-slmsts(nsoil1) ) * vctr6(k) & 
                 + dup * (b1+2.) * smoimid**(b1+1.) ) * ( smoi-smoibotbc ) * vctr6(k) &
           + DD * vctr6(k) &
           + 0.5 * alpha * ( kfdw * (2.*b+3.) * smoimid**(2.*b+2.) + kfup * (2.*b1+3.) * smoimid**(2.*b1+2.) )

     dx = ff/dff
     smoi = smoi - dx
     if (abs(dx) < 1.e-6) exit

    enddo
        
      smoieq(k,i,j) = min(max(smoi,soilcp(nsoil)),.99*slmsts(nsoil))

    enddo


   ENDDO
ENDDO

end subroutine eqsoilmoisture
!******************************************************************************
SUBROUTINE EQSOILMOISTUREiter(n2,js,je,nzg,slz,dz,dtll,soiltxt,landmask,fdepth,smoieq)
use module_parallel
use module_rootdepth

integer :: n2,js,je,nzg,nss,nee,i,j,nsoil,nsoil1,k,kk,iter,flag,freedrain,itermax
real , dimension(nzg+1) :: slz,vt3di
real , dimension(0:nzg+1) :: qlatflux
real , dimension(nzg) :: dz,vctr4,smoi,dsmoi
integer*1, dimension(nzg) :: icefactor
integer, dimension(n2,js:je) :: landmask
real, dimension(n2,js:je) :: fdepth
real, dimension(nzg,n2,js:je) :: smoieq
real :: dtll,smoiwtd,dsmoideep,wtd,rechstep,deeprech,ppdrip,petstep_s,etstep_s,runoff,pppendepth,smoisat,psisat,smoicp,qlat,qrf,qrfcorrect,flood

integer, dimension(2,n2,js:je) :: soiltxt

itermax=500
!if(pid.le.71.or.pid.ge.90)itermax=500000
if(pid.ne.50)itermax=500000


freedrain=0

         do k=1,nzg
            vctr4(k) = 0.5 * (slz(k) + slz(k+1))
         enddo


   nss=js+1
   if(pid.eq.1.or. numtasks.eq.1)nss=js

   nee=je-1
   if(pid.eq.numtasks-2.or.numtasks.eq.1)nee=je

if(pid.eq.numtasks-1)write(6,*)'mirar',nss,nee

DO j=nss,nee
if(pid.eq.80)write(6,*)'now doing j',j
  DO i=1,n2

          if(landmask(i,j).eq.0)cycle

DO k=2,nzg

nsoil=soiltxt(1,i,j)

if(k.gt.1)then
     wtd=slz(k)
else
     wtd=slz(k)-1.e-6
endif
if(k.gt.1)then
     do kk=k,nzg
         smoi(kk)=soilcp(nsoil)*max(min(exp((vctr4(kk)+1.5)/fdepth(i,j)),1.),0.1)
     enddo
     do kk=1,k-1
          smoi(kk)=slmsts(nsoil)*max(min(exp((vctr4(kk)+1.5)/fdepth(i,j)),1.),0.1)
     enddo
endif


flag=0
iter=0
do while(flag.eq.0)
!do iter =1, 2000

      iter=iter+1

dsmoi=0.
dsmoideep=0.
ppdrip=0.
petstep_s=0.
qlat=0.
qlatflux=0.
qrf=0.
qrfcorrect=0.
flood=0.
icefactor=0

      call soilfluxes(i,j,nzg,freedrain,dtll,slz,dz,soiltxt(1,i,j),smoiwtd,dsmoi,dsmoideep &
                      ,smoi(1),wtd,rechstep,deeprech,ppdrip,petstep_s,etstep_s,runoff,vt3di,fdepth(i,j) &
                      ,qlat,qlatflux,qrf,qrfcorrect,flood,icefactor)


!if(k.eq.1)then
!    if(abs(smoiwtd-slmsts(soiltxt(1))).lt.1.e-6)flag=1
!    smoiwtd=slmsts(soiltxt(1))
!else

!    if(1.-abs( smoi(k-1) /  (slmsts(nsoil) * max(min(exp((vctr4(k-1)+1.5)/fdepth(i,j)),1.),0.1))) &
!               .gt.0.999.and.abs(vt3di(k)).lt.1.e-5)flag=1
     if(smoi(k-1).eq.slmsts(nsoil) * max(min(exp((vctr4(k-1)+1.5)/fdepth(i,j)),1. ) ,0.1))flag=1

     do kk=1,k-1
          smoi(kk)=slmsts(nsoil)*max(min(exp((vctr4(kk)+1.5)/fdepth(i,j)),1.),0.1)
     enddo
!endif

if(iter.gt.itermax)exit

if(i.eq.n2/2.and.j.eq.4000.and.k.eq.nzg)write(6,*)'doing iter',iter

enddo


!    write(6,*)'done with layer',nsoil,k,iter,smoi(k),vt3di(k)
    smoieq(k,i,j)=smoi(k)


  ENDDO
           smoisat = slmsts(nsoil)*max(min(exp((vctr4(1)+1.5)/fdepth(i,j)),1.),0.1)
           psisat = slpots(nsoil)*min(max(exp(-(vctr4(1)+1.5)/fdepth(i,j)),1.),10.)
           smoicp = soilcp(nsoil)*max(min(exp((vctr4(1)+1.5)/fdepth(i,j)),1.),0.1)

                smoieq(1,i,j) = max( smoisat * ( psisat / &
                    (psisat - dz(1)) ) ** (1./slbs(nsoil)) , smoicp )

 ENDDO
ENDDO

end subroutine EQSOILMOISTUREiter

!******************************************************************************
subroutine SOILFLUXES_EQSMOI(i,j,nzg,ztop,dtll,slz,dz,soiltxt,smoibotbc)
use module_rootdepth

integer :: nzg,ztop,nsoil,nsoil1,k,iwtd,kwtd,i,j
real, dimension(nzg+1) :: slz
real, dimension(nzg) :: dz,vctr2,vctr4,vctr5,vctr6
real, dimension(nzg) :: smoi,kfmid,diffmid &
                        ,aa,bb,cc,rr
real, dimension(nzg+1) :: vt3di
integer, dimension(2) :: soiltxt
real :: wgpmid,kfup,kfdw,hydcon,smoiwtd,smoibotbc &
     ,fracliqwtd,wmid,wtdold,dzup,vt3dbdw,vt3dcdw,dtll,smoibot,dsmoi,icefac,ddw,dup


         do k=1,nzg
            vctr2(k) = 1. / dz(k)
            vctr4(k) = 0.5 * (slz(k) + slz(k+1))
         enddo
         do k = 2,nzg
            vctr5(k) = vctr4(k) - vctr4(k-1)
            vctr6(k) = 1. / vctr5(k)
         enddo


         kfmid = 0.
         diffmid = 0.

         vt3di = 0.

         do k = 1,ztop

            if(slz(k).lt.-0.30)then
                  nsoil=soiltxt(1)
            else
                  nsoil=soiltxt(2)
            endif

            smoi(k)=slmsts(nsoil)

         enddo

         smoiwtd = smoi(1)


         do k = 2,ztop

!gmmdiffusivity and conductivity at the interlayers

            wgpmid=smoi(k)+(smoi(k)-smoi(k-1))*(slz(k)-vctr4(k))*vctr6(k)

            if(slz(k).lt.-0.30)then
                  nsoil=soiltxt(1)
            else
                  nsoil=soiltxt(2)
            endif

            hydcon=slcons(nsoil)
!            icefac=fracliq(k)** (2. * slbs(nsoil) + 3.)
            icefac=1.
            kfdw =   icefac * hydcon  &
               * (wgpmid  / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)
            ddw =-icefac * (hydcon*slpots(nsoil)*slbs(nsoil)/slmsts(nsoil))  &
                       * (wgpmid/slmsts(nsoil)) **(slbs(nsoil)+2.)

            if(slz(k-1).lt.-0.30)then
                  nsoil=soiltxt(1)
            else
                  nsoil=soiltxt(2)
            endif

            hydcon=slcons(nsoil)
!            icefac=fracliq(k-1)** (2. * slbs(nsoil) + 3.)
            icefac=1.
            kfup =   icefac * hydcon  &
               * (wgpmid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)
            dup =-icefac * (hydcon*slpots(nsoil)*slbs(nsoil)/slmsts(nsoil))  &
                       * (wgpmid/slmsts(nsoil)) **(slbs(nsoil)+2.)

            if(kfdw.gt.0.)then
                 kfmid(k)=0.5*(kfdw+kfup)
            else
                 kfmid(k)=0.
            endif

            if(smoi(k).le.smoi(k-1).and.dup.eq.0.)then
                    diffmid(k)=0.
            elseif(smoi(k).ge.smoi(k-1).and.ddw.eq.0)then
                    diffmid(k)=0.
            else
                    diffmid(k)=0.5*(ddw+dup)
            endif


         enddo


!calculate tridiagonal matrix elements

       do k=2,ztop-1

            aa(k) = diffmid(k)*vctr6(k)
            cc(k) = diffmid(k+1)*vctr6(k+1)
            bb(k) = -( aa(k) + cc(k) + dz(k)/dtll )
            rr(k) = -smoi(k)*dz(k)/dtll -kfmid(k+1) +kfmid(k) 

        enddo

!boundary conditions

!top boundary
        if(ztop.gt.1)then
            aa(ztop) = diffmid(ztop)*vctr6(ztop)
            bb(ztop) = -aa(ztop) -dz(ztop)/dtll
            rr(ztop) = vt3di(ztop+1)/dtll -smoi(ztop)*dz(ztop)/dtll +kfmid(ztop)
        endif

!now bottom boundary condition


!interpolate smoiwtd up only to calculate k and diff

            wgpmid=smoi(1)-(smoi(1)-smoiwtd)*&
                            (slz(1)-vctr4(1))/(slz(1)-0.5*dz(1)-vctr4(1))

            smoibot=smoiwtd


            nsoil=soiltxt(1)
!            if(deeptemp.le.273.16)then
!                  fracliqwtd=0.
!            else
                  fracliqwtd=1.
!            endif

            hydcon=slcons(nsoil)
!            icefac=fracliq(1)** (2. * slbs(nsoil) + 3.)
            icefac=1.
            kfdw =   icefac * hydcon  &
               * (wgpmid  / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)
            ddw = -icefac * (hydcon*slpots(nsoil)*slbs(nsoil)/slmsts(nsoil))  &
                       * (wgpmid/slmsts(nsoil)) **(slbs(nsoil)+2.)

            icefac=fracliqwtd** (2. * slbs(nsoil) + 3.)
            kfup =   icefac * hydcon  &
               * (wgpmid / slmsts(nsoil)) ** (2. * slbs(nsoil) + 3.)
            dup=-icefac * (hydcon*slpots(nsoil)*slbs(nsoil)/slmsts(nsoil))  &
                       * (wgpmid/slmsts(nsoil)) **(slbs(nsoil)+2.)

            if(kfdw.gt.0)then
                 kfmid(1)=0.5*(kfdw+kfup)
            else
                 kfmid(1)=0.
            endif
            if(smoi(1).le.smoibot.and.dup.eq.0.)then
                    diffmid(1)=0.
            elseif(smoi(k).ge.smoi(k-1).and.ddw.eq.0)then
                    diffmid(1)=0.
            else
                    diffmid(1)=0.5*(ddw+dup)
            endif

       if(ztop.gt.1)then
            aa(1) = diffmid(1)*vctr6(2)
            cc(1) = diffmid(2)*vctr6(2)
            bb(1) = -( aa(1) + cc(1) + dz(1)/dtll )
            rr(1) = -smoi(1)*dz(1)/dtll -kfmid(2) +kfmid(1) - aa(1)*smoibot 


!solve tridiagonal system and update smoi

            call tridag(aa,bb,cc,rr,smoi,ztop)

        else
            aa(1) = diffmid(1)*vctr6(2)
            bb(1) = -( aa(1) + dz(1)/dtll )
            rr(1) = -smoi(1)*dz(1)/dtll + kfmid(1) - aa(1)*smoibot

            smoi(1) = rr(1) / bb(1)
        endif

            smoibotbc=smoi(ztop)

end subroutine soilfluxes_eqsmoi

!     ******************************************************************
FUNCTION func(x,bexp,smoisat,hydcon,psisat,slz,vctr4,vctr6,dz,smoi1,flux)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x,bexp,dz,smoi1,flux,smoisat,hydcon,psisat,slz,vctr4,vctr6
REAL(SP) :: func,wgpmid,d1,k1

            wgpmid=x+(x-smoi1)*(slz-vctr4)*vctr6
            wgpmid=min(wgpmid,smoisat)
            k1=hydcon*(wgpmid/smoisat)**(2.*bexp+3.)
            d1=-(hydcon*psisat*bexp/smoisat)  &
                       * (wgpmid/smoisat) **(bexp+2.)

       func = d1 * (x - smoi1) / dz + k1 + flux


END FUNCTION func

!     ******************************************************************
FUNCTION zbrent(x1,x2,tol,bexp,smoisat,hydcon,psisat,slz,vctr4,vctr6,dz,smoi1,flux)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x1,x2,tol,bexp,smoisat,hydcon,psisat,slz,vctr4,vctr6,dz,smoi1,flux
REAL(SP) :: zbrent
!INTERFACE
!FUNCTION func(x,pot,wgp1,bexp,alpha,dz,wsat)
!USE nrtype
!IMPLICIT NONE
!REAL(SP), INTENT(IN) :: x,pot,wgp1,bexp,alpha,dz,wsat
!REAL(SP) :: func
!END FUNCTION func
!END INTERFACE
INTEGER(I4B), PARAMETER :: ITMAX=100
REAL(SP), PARAMETER :: EPS=epsilon(x1)
!Using Brent's method, find the root of a function func known to lie
!between x1
!and x2.
!The root, returned as zbrent, will be refined until its accuracy is
!tol.
!Parameters: Maximum allowed number of iterations, and machine
!floating-point
!precision.
INTEGER(I4B) :: iter
REAL(SP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
   a=x1
   b=x2
   fa=func(a,bexp,smoisat,hydcon,psisat,slz,vctr4,vctr6,dz,smoi1,flux)
   fb=func(b,bexp,smoisat,hydcon,psisat,slz,vctr4,vctr6,dz,smoi1,flux)
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
!write(6,*)'mirar',fa,fb,a,b,bexp,smoisat,hydcon,psisat,slz,vctr4,vctr6,dz,smoi1,flux
!      call nrerror('root must be bracketed for zbrent')
  zbrent=-1.
  RETURN
  endif
  c=b
  fc=fb
      do iter=1,ITMAX
          if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
             c=a !Rename a, b, c and adjust bounding interval d.
             fc=fa
             d=b-a
             e=d
          end if
          if (abs(fc) < abs(fb)) then
             a=b
             b=c
             c=a
             fa=fb
             fb=fc
             fc=fa
          end if
     tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol !Convergence check.
     xm=0.5_sp*(c-b)
          if (abs(xm) <= tol1 .or. fb == 0.0) then
             zbrent=b
             RETURN
          end if
          if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
              s=fb/fa !Attempt inverse quadratic interpolation.
              if (a == c) then
                 p=2.0_sp*xm*s
                 q=1.0_sp-s
              else
                 q=fa/fc
                 r=fb/fc
                 p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
                 q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
              end if
          if (p > 0.0) q=-q !Check whether in bounds.
              p=abs(p)
              if (2.0_sp*p < min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
                 e=d !Accept interpolation.
                 d=p/q
              else
                 d=xm !Interpolation failed; use bisection.
                 e=d
              end if
          else !Bounds decreasing too slowly; use bisection.
              d=xm
              e=d
          end if
     a=b !Move last best guess to a.
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 ) !Evaluate new trial root.
     fb=func(b,bexp,smoisat,hydcon,psisat,slz,vctr4,vctr6,dz,smoi1,flux)
end do
call nrerror('zbrent: exceeded maximum iterations')
zbrent=b
END FUNCTION zbrent
SUBROUTINE nrerror(string)
!Report a message, then die.
CHARACTER(LEN=*), INTENT(IN) :: string
write (*,*) 'nrerror: ',string
STOP 'program terminated by nrerror'
END SUBROUTINE nrerror
!     ******************************************************************
FUNCTION func1(x,ks,bexp,alpha,wsat,qz)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x,ks,bexp,alpha,wsat,qz
REAL(SP) :: func1
func1 = -ks * ((x+wsat)/(2.*wsat))**(2.*bexp+3.) * &
        (alpha*((wsat/x)**bexp - 1.) +1.) -qz
END FUNCTION func1

!     ******************************************************************
FUNCTION zbrent1(x1,x2,tol,ks,bexp,alpha,wsat,qz)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x1,x2,tol,ks,bexp,alpha,wsat,qz
REAL(SP) :: zbrent1
!INTERFACE
!FUNCTION func1(x,ks,bexp,alpha,wsat,qz)
!USE nrtype
!IMPLICIT NONE
!REAL(SP), INTENT(IN) :: x,ks,bexp,alpha,wsat,qz
!REAL(SP) :: func1
!END FUNCTION func1
!END INTERFACE
INTEGER(I4B), PARAMETER :: ITMAX=100
REAL(SP), PARAMETER :: EPS=epsilon(x1)
!Using Brent's method, find the root of a function func known to lie
!between x1
!and x2.
!The root, returned as zbrent, will be refined until its accuracy is
!tol.
!Parameters: Maximum allowed number of iterations, and machine
!floating-point
!precision.
INTEGER(I4B) :: iter
REAL(SP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
   a=x1
   b=x2
   fa=func1(a,ks,bexp,alpha,wsat,qz)
   fb=func1(b,ks,bexp,alpha,wsat,qz)
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
    write(6,*)'mirar',a,b,x1,x2,fa,fb,ks,bexp,alpha,wsat,qz
      call nrerror('root must be bracketed for zbrent1')
  endif
  c=b
  fc=fb
      do iter=1,ITMAX
          if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
             c=a !Rename a, b, c and adjust bounding interval d.
             fc=fa
             d=b-a
             e=d
          end if
          if (abs(fc) < abs(fb)) then
             a=b
             b=c
             c=a
             fa=fb
             fb=fc
             fc=fa
          end if
     tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol !Convergence check.
     xm=0.5_sp*(c-b)
          if (abs(xm) <= tol1 .or. fb == 0.0) then
             zbrent1=b
             RETURN
          end if
          if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
              s=fb/fa !Attempt inverse quadratic interpolation.
              if (a == c) then
                 p=2.0_sp*xm*s
                 q=1.0_sp-s
              else
                 q=fa/fc
                 r=fb/fc
                 p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
                 q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
              end if
          if (p > 0.0) q=-q !Check whether in bounds.
              p=abs(p)
              if (2.0_sp*p < min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
                 e=d !Accept interpolation.
                 d=p/q
              else
                 d=xm !Interpolation failed; use bisection.
                 e=d
              end if
          else !Bounds decreasing too slowly; use bisection.
              d=xm
              e=d
          end if
     a=b !Move last best guess to a.
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 ) !Evaluate new trial root.
     fb=func1(b,ks,bexp,alpha,wsat,qz)
end do
call nrerror('zbrent: exceeded maximum iterations')
zbrent1=b
END FUNCTION zbrent1
!     ******************************************************************
!     ******************************************************************
FUNCTION func2(x,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz
REAL(SP) :: func2,a,b,c
REAL(DP) :: d,e,f
a=wsat/x
b=wsat/(wgp1+(x-wgp1)*slp2)
c=(wgp1+(x-wgp1)*slp1)/wsat
d=dble(ksubs)*dble(c)**(2.*dble(bexp)+3.)
e= dble(alpha)*( dble(b)**dble(bexp) -dble(a)**dble(bexp) ) - 1.
f=d*e-dble(qz)
!func2 = ksubs * ( (wgp1+(x-wgp1)*slp1)/wsat )**(2.*bexp+3.) &
!        * (alpha*( b**bexp -a**bexp) +1.) -qz
func2=sngl(f)

!write(6,*)'mirar',x,a,b,c,d,e,f
END FUNCTION func2

!     ******************************************************************
FUNCTION zbrent2(x1,x2,tol,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x1,x2,tol,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz
REAL(SP) :: zbrent2
!INTERFACE
!FUNCTION func2(x,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz)
!USE nrtype
!IMPLICIT NONE
!REAL(SP), INTENT(IN) :: x,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz
!REAL(SP) :: func2
!END FUNCTION func2
!END INTERFACE
INTEGER(I4B), PARAMETER :: ITMAX=100
REAL(SP), PARAMETER :: EPS=epsilon(x1)
!Using Brent's method, find the root of a function func known to lie
!between x1
!and x2.
!The root, returned as zbrent, will be refined until its accuracy is
!tol.
!Parameters: Maximum allowed number of iterations, and machine
!floating-point
!precision.
INTEGER(I4B) :: iter
REAL(SP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
   a=x1
   b=x2
   fa=func2(a,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz)
   fb=func2(b,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz)
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
!     write(6,*)'puta
!     2',fa,fb,a,b,wgp1,qz,wgp1+(a-wgp1)*slp1,wgp1+(a-wgp1)*slp2,&
!      alpha,alpha*((wsat/a)**bexp -(wsat/(wgp1+(a-wgp1)*slp2))**bexp),&
!       -ksubs * ( (wgp1+(a-wgp1)*slp1)/wsat )**(2.*bexp+3.),&
!       -ksubs * ( (wgp1+(a-wgp1)*slp1)/wsat )**(2.*bexp+3.) * &
!    (alpha*((wsat/a)**bexp -(wsat/(wgp1+(a-wgp1)*slp2))**bexp) + 1.)
!    -qz
      if(qz.gt.0..and.fa.lt.0.)then
         zbrent2=a
         return
      else
      write(6,*)'problema en zbrent2',fa,fb,wgp1,qz,a,b
      call nrerror('root must be bracketed for zbrent2')
      endif
  endif
  c=b
  fc=fb
      do iter=1,ITMAX
          if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
             c=a !Rename a, b, c and adjust bounding interval d.
             fc=fa
             d=b-a
             e=d
          end if
          if (abs(fc) < abs(fb)) then
             a=b
             b=c
             c=a
             fa=fb
             fb=fc
             fc=fa
          end if
     tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol !Convergence check.
     xm=0.5_sp*(c-b)
          if (abs(xm) <= tol1 .or. fb == 0.0) then
             zbrent2=b
             RETURN
          end if
          if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
              s=fb/fa !Attempt inverse quadratic interpolation.
              if (a == c) then
                 p=2.0_sp*xm*s
                 q=1.0_sp-s
              else
                 q=fa/fc
                 r=fb/fc
                 p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
                 q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
              end if
          if (p > 0.0) q=-q !Check whether in bounds.
              p=abs(p)
              if (2.0_sp*p < min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
                 e=d !Accept interpolation.
                 d=p/q
              else
                 d=xm !Interpolation failed; use bisection.
                 e=d
              end if
          else !Bounds decreasing too slowly; use bisection.
              d=xm
              e=d
          end if
     a=b !Move last best guess to a.
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 ) !Evaluate new trial root.
     fb=func2(b,wgp1,ksubs,bexp,wsat,alpha,slp1,slp2,qz)
end do
call nrerror('zbrent: exceeded maximum iterations')
zbrent2=b
END FUNCTION zbrent2
!     ******************************************************************
FUNCTION func3(x,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz
REAL(SP) :: func3
func3 = -0.5 * ( ks1 * (wgp1+alpha*(x-wgp1))**(2.*bexp1+3.) + &
                 ks2 * (wgp1+alpha*(x-wgp1))**(2.*bexp2+3.) ) * &
         ( slp2 * x**(-bexp2) - slp1 + 1. )  -qz
END FUNCTION func3

!     ******************************************************************
FUNCTION zbrent3(x1,x2,tol,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz)
USE nrtype
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x1,x2,tol,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz
REAL(SP) :: zbrent3
!INTERFACE
!FUNCTION func3(x,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz)
!USE nrtype
!IMPLICIT NONE
!REAL(SP), INTENT(IN) :: x,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz
!REAL(SP) :: func3
!END FUNCTION func3
!END INTERFACE
INTEGER(I4B), PARAMETER :: ITMAX=100
REAL(SP), PARAMETER :: EPS=epsilon(x1)
!Using Brent's method, find the root of a function func known to lie
!between x1
!and x2.
!The root, returned as zbrent, will be refined until its accuracy is
!tol.
!Parameters: Maximum allowed number of iterations, and machine
!floating-point
!precision.
INTEGER(I4B) :: iter
REAL(SP) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
   a=x1
   b=x2
   fa=func3(a,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz)
   fb=func3(b,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz)
  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0))then
      if(qz.gt.0..and.fa.lt.0.)then
         zbrent3=a
         return
      else
      write(6,*)'puta',fa,fb,wgp1,qz,a,b
      call nrerror('root must be bracketed for zbrent3')
      endif
  endif
  c=b
  fc=fb
      do iter=1,ITMAX
          if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
             c=a !Rename a, b, c and adjust bounding interval d.
             fc=fa
             d=b-a
             e=d
          end if
          if (abs(fc) < abs(fb)) then
             a=b
             b=c
             c=a
             fa=fb
             fb=fc
             fc=fa
          end if
     tol1=2.0_sp*EPS*abs(b)+0.5_sp*tol !Convergence check.
     xm=0.5_sp*(c-b)
          if (abs(xm) <= tol1 .or. fb == 0.0) then
             zbrent3=b
             RETURN
          end if
          if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
              s=fb/fa !Attempt inverse quadratic interpolation.
              if (a == c) then
                 p=2.0_sp*xm*s
                 q=1.0_sp-s
              else
                 q=fa/fc
                 r=fb/fc
                 p=s*(2.0_sp*xm*q*(q-r)-(b-a)*(r-1.0_sp))
                 q=(q-1.0_sp)*(r-1.0_sp)*(s-1.0_sp)
              end if
          if (p > 0.0) q=-q !Check whether in bounds.
              p=abs(p)
              if (2.0_sp*p < min(3.0_sp*xm*q-abs(tol1*q),abs(e*q))) then
                 e=d !Accept interpolation.
                 d=p/q
              else
                 d=xm !Interpolation failed; use bisection.
                 e=d
              end if
          else !Bounds decreasing too slowly; use bisection.
              d=xm
              e=d
          end if
     a=b !Move last best guess to a.
     fa=fb
     b=b+merge(d,sign(tol1,xm), abs(d) > tol1 ) !Evaluate new trial root.
     fb=func3(b,wgp1,ks1,ks2,bexp1,bexp2,alpha,slp1,slp2,qz)
end do
call nrerror('zbrent3: exceeded maximum iterations')
zbrent3=b
END FUNCTION zbrent3


END MODULE MODULE_INITIAL
