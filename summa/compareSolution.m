function [err] = compareSolution(n)
  fa = fopen("a.m", "r");
  fb = fopen("b.m", "r");
  fc = fopen("c.m", "r");

  A = fread(fa, [n,n], "double");
  B = fread(fb, [n,n], "double");
  C = fread(fc, [n,n], "double");

  err = max(max(abs(A'*B' - C')));
end
