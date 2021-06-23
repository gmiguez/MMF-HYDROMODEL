MODULE module_parallel
implicit none
include 'mpif.h'

!maximumn number of processors
integer, parameter :: npmax=1000

integer, save, dimension (0:npmax) :: domblock,domblock2dint,domblock3d,domblock3dint &
                        ,domblocksmall,domblocksmallint,domblock3dsmall &
                        ,nini,nend,rcountblock,rcountblocksmall,disp,domblockbyte
integer, save :: columntype,columntype2
integer, save :: pid,numtasks
integer :: status(MPI_STATUS_SIZE),ierr

!integer, parameter :: n2big=7320,n3big=8520,nw=4401,ne=4500,ns=4481,nn=4580
!manaosinteger, parameter :: n2big=7320,n3big=8520,nw=4001,ne=4500,ns=6391,nn=6890
!pantanalinteger, parameter :: n2big=7320,n3big=8520,nw=4201,ne=4500,ns=4501,nn=4800
!testinteger, parameter :: n2big=7320,n3big=8520,nw=4601,ne=4700,ns=6301,nn=6400
!integer, parameter :: n2big=7320,n3big=8520,nw=2851,ne=3250,ns=3251,nn=3650
!argentinainteger, parameter :: n2big=7320,n3big=8520,nw=2651,ne=2750,ns=2801,nn=2900
!integer, parameter :: n2big=7320,n3big=8520,nw=2950,ne=3051,ns=2801,nn=2900
!integer, parameter :: n2big=7320,n3big=8520,nw=3901,ne=4800,ns=4101,nn=5000
integer, parameter :: n2big=7320,n3big=8520,nw=1,ne=n2big,ns=1,nn=n3big


CONTAINS

SUBROUTINE INITIALIZEDOMAIN(n2,n3,nzg,filetopo)
integer :: n2,n3,nzg
integer :: n,nmax
integer :: tasktype
integer :: request
integer, allocatable :: req(:),stats(:,:)
character (len = 300) :: filetopo

   call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

if(numtasks.eq.1)then
   nini(0)=1
   nend(0)=n3
   return
endif

rcountblocksmall(0)=n2
rcountblock(0)=0
disp(0)=0

rcountblocksmall(numtasks-1)=n2
rcountblock(numtasks-1)=0
disp(numtasks-1)=0

call MPI_TYPE_CONTIGUOUS(numtasks,MPI_INTEGER,tasktype,ierr)
call MPI_Type_commit(tasktype,ierr)

if(pid.eq.0)then

   nini(0)=0
   nini(numtasks-1)=n3+1

   call dividedomain(n2,n3,nini,filetopo)

   allocate(req(numtasks-1))
   allocate(stats(MPI_STATUS_SIZE,numtasks-1))

   do n=1,numtasks-1
     call MPI_isend(nini(0),1,tasktype,n,1,MPI_COMM_WORLD,req(n),ierr)
   enddo
    if(numtasks.gt.1)call MPI_waitall(numtasks-1,req,stats,ierr)

   deallocate(req,stats)

else
     call MPI_irecv(nini(0),1,tasktype,0,1,MPI_COMM_WORLD,request,ierr)
     call MPI_wait(request,status,ierr)
endif

call MPI_TYPE_FREE (tasktype,ierr)

nend(0)=-1
nend(numtasks-1)=n3
nend(numtasks-2)=n3

do n=2,numtasks-2
  nend(n-1)=nini(n)+1
enddo

!gmmdeclare pieces to be send and received

do n=1,numtasks-2

   nmax=nend(n)-nini(n)+1

if(pid.eq.0)write(6,*)nini(n),nend(n),nmax,n,pid

  rcountblocksmall(n)=n2*(nmax-2)
  if(n.eq.1.or.n.eq.numtasks-2)rcountblocksmall(n)=n2*(nmax-1)
  rcountblock(n)=n2*nmax

  disp(n)=rcountblocksmall(n-1)+disp(n-1)

call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_REAL,domblock(n),ierr)
call MPI_Type_commit(domblock(n),ierr)
call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_INTEGER,domblock2dint(n),ierr)
call MPI_Type_commit(domblock2dint(n),ierr)
call MPI_TYPE_CONTIGUOUS(2*n2*nmax,MPI_INTEGER,domblock3dint(n),ierr)
call MPI_Type_commit(domblock3dint(n),ierr)
call MPI_TYPE_CONTIGUOUS(rcountblocksmall(n),MPI_REAL,domblocksmall(n),ierr)
call MPI_Type_commit(domblocksmall(n),ierr)
call MPI_TYPE_CONTIGUOUS(rcountblocksmall(n),MPI_REAL,domblocksmallint(n),ierr)
call MPI_Type_commit(domblocksmallint(n),ierr)
call MPI_TYPE_CONTIGUOUS(n2*nmax,MPI_BYTE,domblockbyte(n),ierr)
call MPI_Type_commit(domblockbyte(n),ierr)

enddo

   call MPI_Type_CONTIGUOUS(n2,MPI_REAL,columntype,ierr)
   call MPI_Type_commit(columntype,ierr)

   call MPI_Type_CONTIGUOUS(2*n2,MPI_REAL,columntype2,ierr)
   call MPI_Type_commit(columntype2,ierr)


end subroutine initializedomain
!**********************************************************************

subroutine dividedomain(n1,n2,ini,filetopo)
integer :: n1,n2
integer :: ini(0:numtasks-1)
real, dimension(n2big,n3big) :: varreadbig
real, dimension(n1,n2) :: varread
integer,  dimension(n2) :: ncells
integer :: ntotal,ncount,i,j,n
real :: nperpid
character (len = *) :: filetopo

write(6,*)'reading soil data'
!read in topo data to use as mask

      open(21,file=filetopo &
      ,access='direct',convert='big_endian',recl=4*n2big*n3big)


     read(21,rec=9)((varreadbig(i,j),i=1,n2big),j=1,n3big)

     varread(1:n1,1:n2)=varreadbig(nw:ne,ns:nn)

     close(21)

!varread=0.
write(6,*)'test soil data',varread(2,2),varread(n1/2,n2/2)

      ntotal=count(varread>=-1.E+05)

      nperpid=float(ntotal)/float(numtasks-2)

      write(6,*)'total number of land cells',ntotal

      ncells=count(varread>=-1.E+05,1)

      ncount=0
      ini(1)=1
      n=2
      do j=1,n2
         ncount=ncount+ncells(j)
         if(ncount.ge.nint(float(n-1)*nperpid))then
!              write(6,*)ncount,ncount-ncells(j-1),j,n
              ini(n)=j-1
!              ncount=ncells(j)
              n=n+1
         endif
         if(n.eq.numtasks-1)exit
      enddo

!      do n=0,numtasks-1
!         write(6,*)nini(n),n
!      enddo

!test

!      ncount=0
!      n=0
!      do j=nini(10),nini(11)-1
!        n=ncells(j)+n
!        do i=1,n1
!         if(varread(i,j).gt.0.5)ncount=ncount+1
!        enddo
!      enddo

!     write(6,*)ncount,ntotal/numtasks,n


end subroutine dividedomain

!     ******************************************************************
subroutine SENDBORDERS(n2,js,je,wtd,reqsu,reqsd,reqru,reqrd)
integer :: n2,js,je,reqsu,reqsd,reqru,reqrd
real, dimension(n2,js:je):: wtd

if(pid.eq.1)then
    call MPI_isend(wtd(1,je-1),1,columntype,2,200,MPI_COMM_WORLD,reqsu,ierr)
    call MPI_irecv(wtd(1,je),1,columntype,2,201,MPI_COMM_WORLD,reqru,ierr)

elseif(pid.eq.numtasks-2)then
    call MPI_isend(wtd(1,js+1),1,columntype,pid-1,201,MPI_COMM_WORLD,reqsd,ierr)
    call MPI_irecv(wtd(1,js),1,columntype,pid-1,200,MPI_COMM_WORLD,reqrd,ierr)

elseif(pid.gt.1.and.pid.lt.numtasks-2)then
    call MPI_isend(wtd(1,je-1),1,columntype,pid+1,200,MPI_COMM_WORLD,reqsu,ierr)
    call MPI_isend(wtd(1,js+1),1,columntype,pid-1,201,MPI_COMM_WORLD,reqsd,ierr)
    call MPI_irecv(wtd(1,js),1,columntype,pid-1,200,MPI_COMM_WORLD,reqrd,ierr)
    call MPI_irecv(wtd(1,je),1,columntype,pid+1,201,MPI_COMM_WORLD,reqru,ierr)

endif

end subroutine sendborders
!     ******************************************************************
subroutine SENDBORDERSFLOOD(n2,js,je,wtd,reqsu,reqsd,reqru,reqrd)
integer :: n2,js,je,reqsu,reqsd,reqru,reqrd
real, dimension(n2,js:je):: wtd

if(pid.eq.1)then
    call MPI_isend(wtd(1,je),1,columntype,2,200,MPI_COMM_WORLD,reqsu,ierr)
    call MPI_irecv(wtd(1,je-1),1,columntype,2,201,MPI_COMM_WORLD,reqru,ierr)

elseif(pid.eq.numtasks-2)then
    call MPI_isend(wtd(1,js),1,columntype,pid-1,201,MPI_COMM_WORLD,reqsd,ierr)
    call MPI_irecv(wtd(1,js+1),1,columntype,pid-1,200,MPI_COMM_WORLD,reqrd,ierr)

elseif(pid.gt.1.and.pid.lt.numtasks-2)then
    call MPI_isend(wtd(1,je),1,columntype,pid+1,200,MPI_COMM_WORLD,reqsu,ierr)
    call MPI_isend(wtd(1,js),1,columntype,pid-1,201,MPI_COMM_WORLD,reqsd,ierr)
    call MPI_irecv(wtd(1,js+1),1,columntype,pid-1,200,MPI_COMM_WORLD,reqrd,ierr)
    call MPI_irecv(wtd(1,je-1),1,columntype,pid+1,201,MPI_COMM_WORLD,reqru,ierr)

endif

end subroutine sendbordersflood

!     ******************************************************************


END MODULE MODULE_PARALLEL
