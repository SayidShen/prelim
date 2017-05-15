function [ x ] = gtf_admm( y,D,k,lambda, rho, varargin )
%GTF_ADMM Summary of this function goes here
%   ADMM version of the graph trend filtering
% lambda is tuning parameter, rho is the numerical paramter
% D is the incidence matrix with dimension edge times node
% k is the order.
if nargin == 6
    idx=varargin{1};
else
    i=(1-mod(k,2))+1;
    idx=false(size(D,i),1);
end
n=length(y);
m=size(D,1);

if size(D,2)==2% it is indices
    edges1=D(:,1);
    edges2=D(:,2);
    D=sparse(m,n);
    %Make sparse D matrix
    D = sparse((1:m)',edges1,1,m,n,10*m);
    D = D+ sparse((1:m)',edges2,-1,m,n);
else % D is already the incidence matrix
    edges1=zeros(m,1);
    edges2=zeros(m,1);
    for i=1:m
    Iy=find(D(i,:)~=0);
    edges1(i)=Iy(1);
    edges2(i)=Iy(2);
    %[Ix,Iy]=find(D==1);
    %edges2=Iy;
    end
    D=sparse(D);
end

O=construct_O(D,k);

tol_abs=1e-5;
tol_rel=1e-4;

O=O(~idx,:);

% pre-compute L
L=chol(rho*(O'*O)+speye(n),'lower');

% initialization
conv=0;
x=y;
m=size(O,1);
z=zeros(size(O,1),1);
u=z;
iter = 1;
maxiter=1000;
while ~conv
    x=L\(L'\(O'*(rho*z-u)+y));
    
    Ox=O*x;
    z_new=soft_thresh(Ox+u/rho,lambda/rho);
    s= rho*norm(O'*(z_new-z));
    z=z_new;
    
    u=u+rho*(Ox-z);
    r=norm(Ox-z);
    
    eps_pri= sqrt(n)*tol_abs + tol_rel*max(norm(Ox),norm(z));    
    eps_dual= sqrt(m)*tol_abs  + tol_rel*norm(O'*u);
   
    fprintf('[%d] [r, eps_pri]=%f, %f, [s, eps_dual]=%f, %f. \n',...
        iter,r,eps_pri,s, eps_dual);
    if r < eps_pri && s< eps_dual
        conv=1;
        fprintf('converged.\n')
    elseif iter >= maxiter
        conv=1;
        fprintf('Reached maxiter.\n')
    end
    
    
    
    iter=iter+1;
end


end

