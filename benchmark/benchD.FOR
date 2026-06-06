      program ThreePoint
      implicit none
      doubleprecision h
      doubleprecision E_trial
      doubleprecision mass
      doubleprecision xmin, xmax
      integer count_nodes, total_nodes
      integer   n_of_nodes, niter, nitermax, nsteps
      external function count_nodes
      
      mass = 2.0d0
      xmin = -5.5d0
      xmax = -xmin
      
c      write(*,'(1x,a)') 'Enter step size h = '
c      read(*,*) h
      h=0.001d0
      nsteps = int((xmax-xmin)/h)
c      write(*,'(1X,A,I10)') "Number of steps = ", 
       
      nitermax = 60000
      total_nodes = 0.d0
      do 20 niter = 0, nitermax
         E_trial = (1.d0/3.d0) + dble(20*niter)/dble(nitermax)
         n_of_nodes = count_nodes(mass, xmin, xmax, h, E_trial)
         total_nodes = total_nodes + n_of_nodes
20    continue
      
      write(*,*) 'Finished, total_nodes = ', total_nodes
      write(*,'(a, I10, I10,E14.4)') 'nsteps, nitermax, product = ', 
     +            nsteps, nitermax, 
     +            dble(nsteps)*dble(nitermax)
      
      end
      
      function count_nodes(mass, xmin, xmax, h, E)
      implicit none
      integer count_nodes
      doubleprecision xmin, xmax, h, E
      doubleprecision f0, f1, f2
      doubleprecision x, V, mass
      external V 
      
      x = xmin
      f0 = 0.d0

      x = x + h
      f1 = 0.0000001d0
      
      count_nodes = 0
40    continue    
        x = x + h
        f2 = 2.0*f1*(1.0-h*h*mass*(E-V(x))) - f0
        if(sign(1.0d0,f2).ne.sign(1.0d0,f1)) count_nodes=count_nodes+1
        f0 = f1
        f1 = f2
        if( x .ge. xmax) goto 41
      goto 40
41    continue
      
      end

      function V(x)
      implicit none
      doubleprecision x
      doubleprecision V
      
      V = x*x
      
      end
      
