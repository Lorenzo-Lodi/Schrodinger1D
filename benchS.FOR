      program ThreePoint
      implicit none
      real h
      real E_trial
      real mass
      real xmin, xmax
      integer count_nodes, total_nodes
      integer   n_of_nodes, niter, nitermax, nsteps
      external function count_nodes
      
      mass = 2.0
      xmin = -5.5
      xmax = -xmin
      
c      write(*,'(1x,a)') 'Enter step size h = '
c      read(*,*) h
      h=0.001
      nsteps = int((xmax-xmin)/h)
c      write(*,'(1X,A,I10)') "Number of steps = ", 
       
      nitermax = 60000
      total_nodes = 0.
      do 20 niter = 0, nitermax
         E_trial = (1./3.) + real(20*niter)/real(nitermax)
         n_of_nodes = count_nodes(mass, xmin, xmax, h, E_trial)
         total_nodes = total_nodes + n_of_nodes
20    continue
      
      write(*,*) 'Finished, total_nodes = ', total_nodes
      write(*,'(a, I10, I10,E14.4)') 'nsteps, nitermax, product = ', 
     +            nsteps, nitermax, 
     +            real(nsteps)*real(nitermax)
      
      end
      
      function count_nodes(mass, xmin, xmax, h, E)
      implicit none
      integer count_nodes
      real xmin, xmax, h, E
      real f0, f1, f2
      real x, V, mass
      external V 
      
      x = xmin
      f0 = 0.

      x = x + h
      f1 = 0.0000001
      
      count_nodes = 0
40    continue    
        x = x + h
        f2 = 2.0*f1*(1.0-h*h*mass*(E-V(x))) - f0
        if(sign(1.0,f2).ne.sign(1.0,f1)) count_nodes=count_nodes+1
        f0 = f1
        f1 = f2
        if( x .ge. xmax) goto 41
      goto 40
41    continue
      
      end

      function V(x)
      implicit none
      real x
      real V
      
      V = x*x
      
      end
      
