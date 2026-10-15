function f = knotF(x)
  p = 30; m = 150; r = 0.24;
  t = 4*pi*(0:m)/m;
  x0 = (1+0.4*cos(1.5*t)).*cos(t); y0 = (1+0.4*cos(1.5*t)).*sin(t); z0 = 0.5*sin(1.5*t);
  e = (x(1)-x0).^2 + (x(2)-y0).^2 + (x(3)-z0).^2 - r*r + 1;
  f = sum(e.^(-p))^(-1/p) - 1;
end
