for a = {'MS','GMHC'}
  runjob(a{1},3,1,[-1.9 -1.9 -1.1],[1.9 1.9 1.1],[32 32 22],[0 0 0],@knotF,['knot_' a{1} '.pol']);
end
