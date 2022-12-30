      program ThreePoint
      implicit none
      real    h
      real    E_trial, E_k, E_exact, E_low, E_high
      real    mass = 2.0
      real    xmin, xmax, x
      integer   n_of_desired_nodes

      xmin = -5.5
      xmax = -xmin
      h=0.001
      write(*,'(A, I10)') 'Number of steps = ', int( (xmax-xmin)/h )
      
      E_high = -1e10
      E_low  =  1e10
      x = xmin
      do
        if( x >= xmax) exit
        if( V(x) > E_high) E_high = V(x)
        if( V(x) < E_low) E_low = V(x)
        x = x + h
      enddo
      
      E_low  = E_low  - 0.05*(E_high-E_low)
      E_high = E_high + 0.05*(E_high-E_low)
      
      E_trial = (E_high+E_low)/2.
      
      do n_of_desired_nodes = 0, 10
        E_exact = 0.5 + n_of_desired_nodes
        E_k = find_eigenvalue(n_of_desired_nodes, E_low, E_high, E_trial)
        write(*,'(A, ES12.4,I6, 2F22.14, ES14.2)') 'h, n_of_desired_nodes, approx, exact, rel. err. = ', &
           h, n_of_desired_nodes, E_k, E_exact, (E_k - E_exact)/E_exact
      enddo
      
      E_k = V(-1.e10)
      
      end program

      real  function find_eigenvalue(n_of_desired_nodes, E_low_input, E_high_input, E_trial)
      integer   n_of_desired_nodes
      real    E_trial, E_low_input, E_high_input
      real    E, E_high, E_low
      integer   i, imax
      real , parameter   target_relative_error = 2.*1.19209290E-07      
      
      E = E_trial
      E_low = E_low_input
      E_high = E_high_input
      
      imax = 100
      do i=1, imax
          E = (E_low + E_high)*0.5
          if( (E_high - E_low)/abs(E)  < target_relative_error) exit
          if( count_nodes(xmin, xmax, h, E) > n_of_desired_nodes) then
            E_high = E
          else
            E_low = E
          endif
          
      enddo
      
      find_eigenvalue = E
      
      end

      integer function count_nodes(xmin, xmax, h, E)
      real , intent(in)   xmin, xmax, h, E
      real    f0, f1, f2
      real    x
      
          x = xmin
          f0 = 0.

          x = x + h
          f1 = 0.0000001
      
          count_nodes = 0
          do
            x = x + h
            f2 = 2.0*f1*(1.0-h*h*mass*(E-V(x))) - f0
            if( sign(1.0,f2) /= sign(1.0,f1) )  count_nodes = count_nodes +1
            f0 = f1
            f1 = f2
            if( x >= xmax) exit
          enddo
      
      end

      function V(x)
      real , intent(in)   x
      real    V
      
      V = x*x
      
      end    
      
      
