n = 10;
h = 1.0;

G = zeros(n,n);
% initial condition
for i=5:6
  for j=5:6
    G(i,j) = 1;
  end
end
% scaling to grid scale
G = h^2*G;

% calculate the diagonalization of our discrete laplace opperator
o = ones(n,1);
T = spdiags([o,-2*o,o],[-1,0,1],n,n);

[Q,D] = eig(T);

Gm = Q'*G*Q;

Um = zeros(n,n);
for i=1:n
  for j=1:n
    Um(i,j) = Gm(i,j)/(D(i,i)+D(j,j));
  end
end

U = Q*Um*Q';
save out.mat U
