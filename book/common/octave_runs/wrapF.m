function h = wrapF(f)
  % accept both calling conventions Func(n,x) and Func(x)
  h = @(varargin) f(varargin{end});
end
