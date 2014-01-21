function [Useq] = plotPoissonSeq(n)
  fu = fopen('u-seq.m', 'r');
  Useq = fread(fu, [n,n], "double");
  imagesc(Useq);
  print('poissonSeq.eps', '-deps');
end
