function x = grid_system( b, k, rho)
% this function solves (L^k + rho*I ) x = b for L being the laplacian of a
% grid using DCT.

[n1,n2]=size(b);
% get eigenvalues
lambs1 = 4*sin(pi*((1:n1)-1)/(2*n1)).^2;
lambs2 = 4*sin(pi*((1:n2)-1)/(2*n2)).^2;
lambs = lambs1'*ones(1,n2) + ones(n1,1)*lambs2;
sigma = rho*lambs.^k + ones(n1,n2);

tmp=dct2(b);
tmp=tmp.*(1./sigma);
x=idct2(tmp);

end
