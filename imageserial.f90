program imageserial
        use pgmio

        implicit none 
        
        character(len=64)                                    :: filename = 'edge192x128.pgm'
        integer                                              :: nx, ny
        double precision ,dimension(:,:), allocatable        :: buf
        double precision, dimension(:,:), allocatable        :: old, new, edge
        integer                                              :: i, j, it
        integer, parameter                                   :: maxiter=1500, printfreq=150
        real                                                 :: start, finish  
        double precision                                     :: tol=1.0d-4
  

        !read data size
        call pgmsize(filename, nx, ny)

        write(*,*) 'Processing ', nx, ' x ' , ny, ' image'
        write(*,*) 'Number of iterations = ', maxiter
        
     
        !allocate buffer, old and new
        allocate(buf(1:nx, 1:ny     ))
        allocate(old(0:nx+1, 0:ny+1 ))
        allocate(new(0:nx+1, 0:ny+1 ))
        allocate(edge(0:nx+1, 0:ny+1))

        !read  the edges data file into the buf array
        write(*,*)
        call pgmread(filename, buf)
        write(*,*)
   
        call cpu_time(start)

        !put buff to edge
        do j = 1, ny 
           do i = 1, nx
              edge(i,j) = buf(i,j)
           end do 
        end do
        
        !set the entire old array to 255  include the halos 
        do j = 0, ny+1
           do i = 0, nx+1
              old(i,j) = 255.0
           end do
        end do

        do it = 1, maxiter
              
              !if(mod(it, printfreq) == 0) then
              !       write(*,*) "iteration =", it
              !end if

              do j = 1, ny
                   do i = 1, nx
                        new(i,j) = 0.25*(old(i-1,j) + old(i+1,j) + old(i,j-1) + old(i,j+1) - edge(i,j))
                   end do
              end do
              
              if(mod(it, printfreq) == 0) then
                     write(*,*) "iteration =", it, "---> stop criterion =",  maxval(new-old)
              end if
              
              !compute stop criterion
              if (maxval(new-old) <= tol) exit
              
              !set the old array equal to new
              do j = 1, ny
                   do i = 1, nx
                       old(i,j) = new(i,j)
                   end do
              end do 
        end do
        
        write(*,*)
        write(*,*) "Finishing", it-1, "iterations"
        write(*,*)
        !copy the old array back to buf excluding the halos
        do j = 1, ny
             do i = 1, nx
                 buf(i,j) = old(i,j)
              end do
        end do
         
        call cpu_time(finish)
        !write out the final image by passing buf to pgmwrite
        filename='image192×128.pgm'
        call pgmwrite(filename, buf)

        print '("Time = ",f6.3," seconds.")', finish-start

        !deollocate data
        deallocate(old)
        deallocate(new)
        deallocate(buf)
        deallocate(edge)
end program imageserial
