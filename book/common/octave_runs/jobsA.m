sph3 = @(x) x(1)^2+x(2)^2+x(3)^2-1;
for a = {'MS','MHC','GMHC','CS','GCHC'}
  runjob(a{1},3,1,[-1.1 -1.1 -1.1],[1.1 1.1 1.1],[12 12 12],[1 0 0],sph3,['sph3_' a{1} '.pol']);
end
sph4 = @(x) x(1)^2+x(2)^2+x(3)^2+x(4)^2-1;
for a = {'GMHC','GCHC','CS','MS'}
  runjob(a{1},4,1,-1.1*ones(1,4),1.1*ones(1,4),[8 8 8 8],[1 0 0 0],sph4,['sph4_' a{1} '.pol']);
end
