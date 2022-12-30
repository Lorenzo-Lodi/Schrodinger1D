      program ThreePoint
      implicit none
      real h
      real E_trial, E_k, E_exact, E_low, E_high
      real mass
      real xmin, xmax, x
      integer   n_of_desired_nodes
      real find_eig, V
      external function find_eig
      
      mass = 2.0
      xmin = -5.5
      xmax = -xmin
      
      write(*,'(1x,a)') 'Enter step size h = '
      read(*,*) h
c      h=0.01
      write(*,'(1X,A,I10)') "Number of steps = ",int((xmax-xmin)/h)
      
      E_high = -1e10
      E_low  =  1e10
      x = xmin
      
10    continue
        if( x .ge. xmax) goto 11
        if( V(x) .gt. E_high) E_high = V(x)
        if( V(x) .lt. E_low)  E_low = V(x)
        x = x + h
      goto 10
11    continue    
      
      E_low  = E_low  - 0.05*(E_high-E_low)
      E_high = E_high + 0.05*(E_high-E_low)
      
      E_trial = (E_high+E_low)/2.
      write(*,*) E_low, E_high, E_trial
      
      write(*,'(A)') " h, n_desired_nodes, approx, exact, rel.err."
      
      do 20 n_of_desired_nodes = 0, 20
        E_exact = 0.5 + real(n_of_desired_nodes)
        E_k = find_eig(mass, xmin, xmax, h, n_of_desired_nodes, 
     +        E_low, E_high,E_trial)
        write(*,'(F12.6, I6, 2F23.14, E15.2)') 
     +  h, n_of_desired_nodes, E_k, E_exact, (E_k - E_exact)/E_exact
20    continue
      
      E_k = V(-1.e10)
      
      end
      
      real function find_eig(mass, xmin, xmax, h, 
     + n_of_desired_nodes, E_low_input, E_high_input, E_trial)
      integer   n_of_desired_nodes
      real E_trial, E_low_input, E_high_input
      real E, E_high, E_low
      integer   i, imax
      real target_relative_error
      integer count_nodes
      real xmin, xmax, h, mass
      parameter   (target_relative_error = 2.*1.19209290E-07)      
      external count_nodes
      
      E = E_trial
      E_low = E_low_input
      E_high = E_high_input
      
      imax = 100
      do 30 i=1, imax
          E = (E_low + E_high)*0.5
          if( (E_high - E_low)/abs(E) .lt. target_relative_error) exit
          if( count_nodes(mass, xmin, xmax, h, E) .gt. 
     +     n_of_desired_nodes) then
            E_high = E
          else
            E_low = E
          endif
30    continue
      
      find_eig = E
      
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
      
