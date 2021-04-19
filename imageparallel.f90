program imageparallel
        use mpi
        use pgmio

        implicit none 
        
        character(len=64)                                    :: filename = 'edge192x128.pgm'
        integer                                              :: nx, ny
        double precision, dimension(:,:), allocatable        :: masterbuf, buf
        double precision, dimension(:,:), allocatable        :: old, new, edge
        integer                                              :: i, j, it
        integer, parameter                                   :: maxiter=1500, printfreq=150
        double precision                                     :: criterion, criterion_loc, tol=1.0d-1 

        !MPI variables
        integer                                              :: rank, size, comm, comm1D, ierror
        integer                                              :: tag=1, master = 0, nbrequest = 4
        integer, dimension(MPI_STATUS_SIZE)                  :: statut
        integer, dimension(MPI_STATUS_SIZE, 4)               :: status
        integer, dimension(4)                                :: requests
        integer, dimension(2)                                :: sizes, subsizes, starts 
        integer                                              :: type_overlap
        !Cartesian topology variable
        integer, parameter                                   :: ndims=1, up=1, down=2
        logical                                              :: reorder
        integer                                              :: direction, step
        integer, dimension(2)                                :: neighbor
        integer, dimension(ndims)                            :: dims, coords
        logical, dimension(ndims)                            :: periods


        !local variables
        integer                                              :: ny_loc
        
        !time variables
        double precision                                     :: tStart, tStop

        !intialisation MPI
        call MPI_INIT(ierror)
        
        comm = MPI_COMM_WORLD

        !call MPI_COMM_RANK(comm, rank, ierror)
        call MPI_COMM_SIZE(comm, size, ierror)

        !initialization cartesian topology
        periods(1)=.false. ; reorder=.false.
        dims(1)=0
        
        call MPI_DIMS_CREATE(size, ndims, dims, ierror)

        !create cartesian topology
        call MPI_CART_CREATE(comm, ndims, dims, periods, reorder, comm1D, ierror)

        !processor rank in comm2D
        call MPI_COMM_RANK(comm1D, rank, ierror)

        !initialize neighbor on MPI_PROC_NULL
        neighbor(:) = MPI_PROC_NULL

        !research neighbors West and East
        direction = 0; step = 1
        call MPI_CART_SHIFT(comm1D, direction, step, neighbor(down), neighbor(up), ierror)
       ! print*, "processus", rank, "mon voisin en haut", neighbor(up)
       ! print*, "processus", rank, "mon voisin en bas", neighbor(down)

       
        !read data size
        if (rank == master) then
                call pgmsize(filename, nx, ny)

                write(*,*) 'Processing ', nx, ' x ' , ny, ' image'
                write(*,*) 'Number of iterations = ', maxiter
        endif

        

        !scatter data size amongst the processes
        call MPI_BCAST(nx, 1, MPI_INTEGER, master, comm1D, ierror)
        call MPI_BCAST(ny, 1, MPI_INTEGER, master, comm1D, ierror)
        
        !compute ny_local
        ny_loc = ny/size
     
        !master process allocate masterbuffer, old and new
        allocate(masterbuf(1:nx, 1:ny))
     
        allocate(buf( 1:nx,       1:ny_loc    ))
        allocate(old( 0:nx+1,     0:ny_loc+1  ))
        allocate(new( 0:nx+1,     0:ny_loc+1  ))
        allocate(edge(0:nx+1,     0:ny_loc+1  ))

        !create a new type 
        sizes(1) = nx + 2
        sizes(2) = ny + 2

        starts(:) = 0

        subsizes(1) = nx
        subsizes(2) = 1
       ! call MPI_type_create_subarray(2, sizes, subsizes, starts, &
       !                               MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
       !                               type_overlap, ierror)

       ! call MPI_type_commit(type_overlap, ierror)
        
        !tStart = MPI_WTIME()

        !read  the edges data file into the buf array
        if (rank == master) then
                write(*,*)
                call pgmread(filename, masterbuf)
                write(*,*)
        endif
        
        call MPI_BARRIER(comm1D, ierror)
        tStart = MPI_WTIME()

        !scatter masterbuf amongst process
       ! call MPI_SCATTER(masterbuf, nx*ny_loc, MPI_DOUBLE_PRECISION, & 
       !                  buf,       nx*ny_loc, MPI_DOUBLE_PRECISION, &
       !                  master   , comm1D, ierror)

        call MPI_SCATTER(masterbuf, nx*ny_loc, MPI_DOUBLE_PRECISION, &
                         buf,       nx*ny_loc, MPI_DOUBLE_PRECISION, &
                         master   , comm1D, ierror)

        !put buff to edge
        do j = 1, ny_loc 
           do i = 1, nx
              edge(i,j) = buf(i,j)
           enddo 
        enddo
        
        !set the entire old array to 255  include the halos 
        do j = 0, ny_loc+1
           do i = 0, nx+1
              old(i,j) = 255.0
           enddo
        enddo
        
        
        do it = 1, maxiter
              
              
              !send to the top, receive from the bottom
              call MPI_SENDRECV(old(1, ny_loc),   nx, MPI_DOUBLE_PRECISION, neighbor(up),   tag, & ! send row before ghost
                                old(1, 0),        nx, MPI_DOUBLE_PRECISION, neighbor(down), tag, & ! receive 1st row (ghost)
                                comm1D, statut, ierror)
              !send to the bottom, receive from the  up
              call MPI_SENDRECV(old(1, 1),        nx, MPI_DOUBLE_PRECISION, neighbor(down), tag, & ! send 1st row after ghost 
                                old(1, ny_loc+1), nx, MPI_DOUBLE_PRECISION, neighbor(up),   tag, & ! receive last row  ghost
                                comm1D, statut, ierror)

             ! call MPI_IRECV(old(1, 0)      , 1, type_overlap, neighbor(down), tag, comm1D, requests(1), ierror)
             ! call MPI_ISEND(old(1, ny_loc) , 1, type_overlap, neighbor(up),   tag, comm1D, requests(2), ierror) 

             ! call MPI_IRECV(old(1, ny_loc+1), 1, type_overlap, neighbor(up),  tag, comm1D, requests(3), ierror)
             ! call MPI_ISEND(old(1, 1)       , 1, type_overlap, neighbor(down),tag, comm1D, requests(4), ierror)
             
             ! call MPI_WAITALL(nbrequest, requests, status, ierror)

              do j = 1, ny_loc
                   do i = 1, nx
                        new(i,j) = 0.25*(old(i-1,j) + old(i+1,j) + old(i,j-1) + old(i,j+1) - edge(i,j))
                   enddo
              enddo
              
              !stop criterion compute
              criterion_loc = maxval(new - old)
              
              !max criterion in global 
              call MPI_ALLREDUCE(criterion_loc, criterion, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm1D, ierror) 
              
              if(criterion < tol) exit
         
              if (rank == master) then
                if(mod(it, printfreq) == 0) then
                        write(*,*) "iteration =", it, "====> stop criterion =", criterion
                endif
              endif

              !set the old array equal to new
              do j = 1, ny_loc
                   do i = 1, nx
                       old(i,j) = new(i,j)
                   enddo
              enddo 
        enddo
        
        if (rank == master) then
                write(*,*)
                write(*,*) "Finishing", it-1, "iterations"
                write(*,*)
        endif
        
        !copy the old array back to buf excluding the halos
        do j = 1, ny_loc
             do i = 1, nx
                 buf(i,j) = old(i,j)
              enddo
        enddo
        
        !gather the data from all the buf arrays back to masterbuf
        call MPI_GATHER(buf,       nx*ny_loc, MPI_DOUBLE_PRECISION, &
                        masterbuf, nx*ny_loc, MPI_DOUBLE_PRECISION, &
                        master,  comm1D, ierror)   
        
        call MPI_BARRIER(comm1D, ierror)
        tStop = MPI_WTIME()  
        !write out the final image by passing buf to pgmwrite
        if (rank == master) then 

            filename='image192×128.pgm'
            call pgmwrite(filename, masterbuf)

        endif

       ! tStop =  MPI_WTIME()

        if (rank == master) write(*,*) "Runtime:" , tStop - tStart


        !deollocate data
        deallocate(masterbuf)
        deallocate(old)
        deallocate(new)
        deallocate(buf)
        deallocate(edge)
        
       ! call MPI_Type_free(type_overlap, ierror)
        call MPI_Comm_free(comm1D, ierror)

        call MPI_FINALIZE(ierror)

end program imageparallel
