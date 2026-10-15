function r = random(dist, mu, sigma, m, n)
  % shim for the Statistics toolbox: only 'Normal' is used by the codes
  r = mu + sigma*randn(m, n);
end
