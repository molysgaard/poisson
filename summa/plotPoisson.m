function [U] = plotPoisson(n)
  fu = fopen('u.m', 'r');
  U = fread(fu, [n,n], 'double');
  imagesc(U);
  print('poissonPar.eps', '-deps');
end
