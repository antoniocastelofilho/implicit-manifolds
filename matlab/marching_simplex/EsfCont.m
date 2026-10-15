function EsfCont
  
    First      = [-0.6 -0.6 -0.6 -1.1 -1.1];
    Last       = [ 0.6  0.6  0.6  1.1  1.1];
    Division   = [ 4 4 4 4 4]; 
    FirstPoint = [ 0.5  0.0  0.0  1.0 0.0];
    
    ContinuationSimplex(5, 2, First, Last, Division, FirstPoint, @S2xS1, 'S2xS1_CS.pol');
    
    return 

    function [f] = reta2(n,v) 
      f(1) = v(1) - 0;
      return
    end 

    function [f] = esfera2(n,v) 
      f(1) = v(1)*v(1) + v(2)*v(2) - 1;
      return
    end 

    function [f] = esfera3(n,v) 
      f(1) = v(1)*v(1) + v(2)*v(2) + v(3)*v(3) - 1;
      return
    end 

    function [f] = esfera4(n,v) 
      f(1) = v(1)*v(1) + v(2)*v(2) + v(3)*v(3) + v(4)*v(4) - 1;
      return
    end 

    function [f] = esfera5(n,v) 
      f(1) = v(1)*v(1) + v(2)*v(2) + v(3)*v(3) + v(4)*v(4) + v(5)*v(5) - 1;
      return
    end 

    function [f] = S2xS1(n,v) 
      f(1) = v(1)*v(1) + v(2)*v(2) + v(3)*v(3) - 0.25;
      f(2) = v(4)*v(4) + v(5)*v(5) - 1.0;
      return
    end
    
    function [f] = kleinBottle(n,v) 
      R = 1;
      P = 0.5;
      
      theta = v(1);
      ve = v(2);
      weird_e = v(3);
      
      f(1) = R * (cos(theta/2) * cos(ve) - sin(theta/2) * sin(2*ve)) - v(4);
      f(2) = R * (sin(theta/2) * cos(ve) + cos(theta/2) * sin(2*ve)) - v(5);
      f(3) = P * cos(theta) * (1 + weird_e * sin(ve)) - v(6);
      f(4) = P * sin(theta) * (1 + weird_e * sin(ve)) - v(7);
    end

 end 