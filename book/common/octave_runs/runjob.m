function runjob(alg, n, k, First, Last, D, P0, f, fname)
  % uniform wrapper: alg in {MS, CS, MHC, CHC, GMHC, GCHC}; writes timing line
  F = wrapF(f);
  t0 = tic;
  switch alg
    case 'MS',   MarchingSimplex(n,k,First,Last,D,F,fname);
    case 'CS',   ContinuationSimplex(n,k,First,Last,D,P0,F,fname);
    case 'MHC',  MarchingHyperCube(n,First,Last,D,F,fname);
    case 'CHC',  ContinuationHypercube(n,k,First,Last,D,P0,F,fname);
    case 'GMHC', GeneralizedMarchingHyperCube(n,k,First,Last,D,fname,F);
    case 'GCHC', GeneralizedContinuationHyperCube(n,k,First,Last,D,P0,F,fname);
  end
  t = toc(t0);
  fid = fopen('timings.txt','a'); fprintf(fid,'%s %s n=%d k=%d D=%s %.2f\n', fname, alg, n, k, mat2str(D), t); fclose(fid);
end
