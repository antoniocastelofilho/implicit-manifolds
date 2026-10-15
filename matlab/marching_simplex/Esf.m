function Esf
  
    First    = [-1.1 -1.1 -0.6 -0.6];
    Last     = [1.1 1.1 0.6 0.6];
    Division = [10 10 6 6]; 
  
    MarchingSimplex(4, 2, First, Last, Division, @torus4, 'torus4_MS.pol');
    
    return 

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

    function [f] = torus4(n,v) 
      f(1) = v(1)*v(1) + v(2)*v(2) - 1;
      f(2) = v(3)*v(3) + v(4)*v(4) - 1/4;
      return
    end 

 end 