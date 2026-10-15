F = wrapF(@(x) x(1)^2+x(2)^2+x(3)^2-1);
First=[-1.1 -1.1 -1.1]; Last=-First; D=[8 8 8]; P0=[1 0 0];
tic; MarchingSimplex(3,1,First,Last,D,F,'s3_MS.pol'); t(1)=toc;
tic; ContinuationSimplex(3,1,First,Last,D,P0,F,'s3_CS.pol'); t(2)=toc;
tic; MarchingHyperCube(3,First,Last,D,F,'s3_MHC.pol'); t(3)=toc;
tic; ContinuationHypercube(3,1,First,Last,D,P0,F,'s3_CHC.pol'); t(4)=toc;
tic; GeneralizedMarchingHyperCube(3,1,First,Last,D,'s3_GMHC.pol',F); t(5)=toc;
tic; GeneralizedContinuationHyperCube(3,1,First,Last,D,P0,F,'s3_GCHC.pol'); t(6)=toc;
fprintf(2,'TIMES %g %g %g %g %g %g\n',t);
