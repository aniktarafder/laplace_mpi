!*************************************************
! Laplace MPI Fortran Version
!
! Temperature is initially 0.0
! Boundaries are as follows:
!                                  ________________
!                                 |       3       |
!      0 _________100             |_______________| 
!       |        | 100            |       2       |
!       |        |                |_______________|
!  T=0  |        |                |       1       |
!       |        |                |_______________|
!       |________|                |       0       |
!       0  T=0   0                |_______________|
!                                             
! Each Processor works on a sub grid and then sends its
! boundaries to neighbours
!
!  John Urbanic, PSC 2014
!  Modified by Anik  2021 (Same problem as John, just rotate 90 degree clockwise)
!*************************************************
program mpi
      implicit none

      include    'mpif.h'

      !Size of plate   ! Be sure to use nx>=128, ny>=128 otherwise solution doesn't converge.
      integer, parameter             :: nx=128
      integer, parameter             :: ny=128

      !these are the new parameters for parallel purposes
      integer, parameter             :: total_pes=4
      integer, parameter             :: nydim=ny/total_pes
      integer, parameter             :: bottom=100, top=101

      !usual mpi variables
      integer                        :: mype, npes, ierr
      integer                        :: status(MPI_STATUS_SIZE)

      double precision, parameter    :: max_temp_error=0.01

      integer                        :: i, j, max_iterations, iteration=1
      double precision               :: temp_err, temp_err_global=100.0
      real                           :: start_time, stop_time

      double precision               :: temperature(nx+2, nydim+2), temperature_last(nx+2, nydim+2)
      real                           :: Lx, Ly, hxi, hyi, time
      double precision, parameter    :: dt=0.0001, last_time=10.0
      double precision               :: grid_topo(4)	  
      !usual mpi startup routines
      call MPI_Init(ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, npes, ierr) ! This will give amount of processes requested for the job and will be stored in npes. You will give this by mpirun -n 8. 8 is npes.
      call MPI_Comm_rank(MPI_COMM_WORLD, mype, ierr) ! This will give rank or idemtification number of each processor involved in the work. i.e 0,1,2 etc

      !It is nice to verify that proper number of PEs are running
      if ( npes /= total_pes ) then
         if( mype == 0 ) then
            print *,'This example is hardwired to run only on ', total_pes, ' PEs'
         endif
         call MPI_Finalize(ierr)
         stop
      endif

      Lx = 10.0
      Ly = 10.0
	  Ly = Ly/npes
      hxi= (nx)/Lx           
      hyi= (nydim)/Ly	  
      time=0.0

      !Only one PE should prompt user
      if( mype == 0 ) then
         print*, 'Temperature array', size(Temperature)
         print*, 'Maximum iterations [100-4000]?'
         read*,   max_iterations
      endif
      ! Now, as you have only taken this input user information to PE 0. Other PEs need to know that information. So, you need to send this to everyone. Bcast is the efficient way to do so.
      ! Other PEs need to recieve this information
      call MPI_Bcast(max_iterations, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) ! Basically this tells send max_iterations,1, integer value that is in rank or id 0 processor to all the processors that are in the communication world.

      call cpu_time(start_time)

      call initialize(temperature_last, npes, mype, nx, ny, total_pes, nydim )
!      call initialize(temperature, npes, mype, nx, ny, total_pes, nydim )

      call print_field(temperature_last, iteration, nx+2, nydim+2, hxi, hyi, time, npes, mype, Lx, Ly)

      !do until global error is minimal or until maximum steps
      do while ( temp_err_global > max_temp_error .and. iteration <= max_iterations)
! Always make sure that your first index (here i is the 1st index and j is the 2nd) is the inner most loop. It will increase performance in fortran. In C it's the total opposite the 2nd index should be the outer most loop.
         do j=2,nydim+1
            do i=2,nx+1
               temperature(i,j)=0.25*(temperature_last(i+1,j)+temperature_last(i-1,j)+ &   
                                      temperature_last(i,j+1)+temperature_last(i,j-1) )
            enddo
         enddo

         ! COMMUNICATION PHASE
         ! send messages to right
		 ! This whole thing is simple. Think as a parallel boundary. Coloumn nxdim is sending to coloumn 1. nxdim+1 is receiving from 2.
         ! PE 0,1,2 will send message to the right.
		 ! In Fortran 1st index is contiguous. So, Becareful while sending an array of information. You will give the address of the 1st value of the array and number of information that needs to be send.
		 ! Fortran will move along 1st index 2 --> nx contiguously.
         ! Send message to Top		 
         if (mype < npes-1) then
            call MPI_Send(temperature(2,nydim+1), nx, MPI_DOUBLE_PRECISION, &
                          mype+1, TOP, MPI_COMM_WORLD, ierr)
         endif
         if (mype /= 0) then
            call MPI_Recv(temperature_last(2,1),  nx, MPI_DOUBLE_PRECISION, &     
                          MPI_ANY_SOURCE, TOP, MPI_COMM_WORLD, status, ierr)
         endif

         ! send messages to bottom
         if (mype /= 0) then
            call MPI_Send(temperature(2,2), nx, MPI_DOUBLE_PRECISION, &
                          mype-1, BOTTOM, MPI_COMM_WORLD, ierr)
         endif
         if (mype /= npes-1) then
            call MPI_Recv(temperature_last(2,nydim+2), nx, MPI_DOUBLE_PRECISION, &
                          MPI_ANY_SOURCE, BOTTOM, MPI_COMM_WORLD, status, ierr)
         endif

         temp_err=0.0

         do j=2,nydim+1
            do i=2,nx+1
               temp_err = max( abs(temperature(i,j) - temperature_last(i,j)), temp_err )
               temperature_last(i,j) = temperature(i,j)
            enddo
         enddo
		 
		 
!         do j=1,nydim+2
!            do i=1,nx+2
!           
!               temperature_last(i,j) = temperature(i,j)
!            enddo
!         enddo		 
		 
		 
		 
		 

         !Need to determine and communicate maximum error
         call MPI_Reduce(temp_err, temp_err_global, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
         call MPI_Bcast(temp_err_global, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

         !periodically print test values - only for PE in lower corner. That means last pe.
		 time = dt*iteration
         if( mod(iteration,100).eq.0 ) then
            if( mype == npes-1 ) then
               call track_progress(temperature_last, iteration, npes, mype, nx, ny, total_pes, nydim)
            endif
			
!			if( mype == 2 ) then
!               print*, 'magic point=', temperature(901,251)
!            endif

			call print_field(temperature_last, iteration, nx+2, nydim+2, hxi, hyi, time, npes, mype, Lx, Ly)
         endif

         iteration = iteration+1

      enddo

      !Slightly more accurate timing and cleaner output
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      call cpu_time(stop_time)

      if( mype == 0 ) then
         print*, 'Max error at iteration ', iteration-1, ' was ',temp_err_global
         print*, 'Total time was ',stop_time-start_time, ' seconds.'
      endif

      call MPI_Finalize(ierr)

end program mpi

!Parallel version requires more attention to global coordinates
subroutine initialize(temperature_last, npes, mype, nx, ny, total_pes, nydim )
      implicit none


      integer                        :: nx, ny
	  integer                        :: total_pes
	  integer                        :: nydim
      integer                        :: i,j
      integer                        :: npes,mype

      double precision               :: temperature_last(nx+2, nydim+2)
      double precision               :: tmin, tmax

      temperature_last = 0

      !Top and Bottom Boundaries      
      if( mype == 0 ) then
         do i=1,nx+2
            temperature_last(i,1) = 0.0
         enddo
      endif
      if( mype == npes-1 ) then
         do i=1,nx+2
            temperature_last(i,nydim+2) = (100.0/nx) * (i-1)
         enddo
      endif
      !Left and Right Boundaries      
      tmin =  mype    * 100.0/npes
      tmax = (mype+1) * 100.0/npes
      do j=1,nydim+2
         temperature_last(1,j)       = 0.0
         temperature_last(nx+2,j)    = tmin + ((tmax-tmin)/(nydim)) * (j-1)
      enddo

end subroutine initialize

subroutine track_progress(temperature, iteration, npes, mype, nx, ny, total_pes, nydim)
      implicit none

      integer             :: nx, ny
      integer             :: nydim
      integer             :: mype, npes, total_pes
      integer             :: i,iteration

      double precision    :: temperature(nx+2, nydim+2)

!Parallel version uses global coordinate output so users don't need
!to understand decomposition
      print *, '---------- Iteration number: ', iteration, ' ---------------'
      do i=5,1,-1
         write (*,'("("i4,",",i4,"):",f6.2,"  ")',advance='no') &
                   nx-i,ny-i,temperature(nx+1-i,nydim+1-i)
      enddo
      print *
end subroutine track_progress


subroutine print_field(temperature, iteration, nxp2, nyp2, hxi, hyi, time, npes, mype, Lx, Ly)
      implicit none

      integer             :: nxp2
      integer             :: nyp2
      integer             :: i,j,iteration
	  logical             :: fext
	  real                :: hxi, hyi, time, Lx, Ly
	  character*24        :: temperature_field
	  character*2         ::  num(4)
      data (num(i),i=1,4)/'01','02','03','04'/
      real                :: offsety	  
	  
	  
	  integer             :: npes, mype

      double precision    :: temperature(nxp2, nyp2)
      
	  offsety = mype*Ly

	  
	  temperature_field = ' '
	  temperature_field(1:18) = 'temperature_field '
	  temperature_field(18:20)=  num(mype+1)
	  temperature_field(20:24)= '.dat'

      INQUIRE(FILE=Temperature_field, EXIST=fext)
      if(.not.fext) then

         open(30,file=Temperature_field,form='formatted')

         write(30,*)'VARIABLES = "X","Y","Temperature"'      

         write(30,*)'ZONE ','T= "',time,'" I=',nxp2,' J=',nyp2,&
        ' F=BLOCK',', SOLUTIONTIME=',time
     
		  
         write(30,'(8e16.8)')(((DBLE(i-1))/hxi,i=1,nxp2),j=1,nyp2)
         write(30,'(8e16.8)')(((DBLE(j-1))/hyi+offsety,i=1,nxp2),j=1,nyp2)   
         write(30,'(8e16.8)')temperature


       close(30)
      else
       open(30,file=Temperature_field,form='formatted',status='old', &
               position='append')
      write(30,*) 'ZONE ','T="',time,'" I=',nxp2,' J=',nyp2, &
       ' F=BLOCK, ','VARSHARELIST=([1-2]=1)',', SOLUTIONTIME=',time
     
         write(30,'(8e16.8)')temperature 
	  
       close(30)
      endif    	  
	  
	  

end subroutine print_field
