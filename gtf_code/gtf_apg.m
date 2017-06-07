function [ x, history] = gtf_apg( y,D,k,lambda, varargin )
addpath('TVexact')
%GTF_APG Summary of this function goes here
%   Accelerated proximal gradient that solves
% 1. project y into the nullspace of O. Solve the least square for that
% part, x'.
% 2. For the remaining part of y: y' solve the following
% min_w 0.5\|y'- O^+ w\|^2 + \lambda |w|_1
%   using optimal accelerated proximal gradient algorithms, or other first order methods.
%   Construct gradient using Laplacian solvers
% 3. Get the solution x= O^+ w +  x';

%GTF_APG Summary of this function goes here
%   APG version of the graph trend filtering
% lambda is tuning parameter, rho is the numerical paramter
% D is the incidence matrix with dimension edge times node
% k is the order.

%the k notation is different from signal processing literature by 1
%k=0:  O=D
%k=1:  O=L
%k=2:  O=DL





flag=~mod(k,2);% if k is even, then we have a D in front

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

% Choose the right algorithm
if k==0 % direct solution
    fprintf('k=0: Direct Solution.\n');
    x=graphtv(y,edges1,edges2,lambda);
    history=0;
    return;
end


history=[];

L=D'*D;
pfun=cmg_sdd(L);

Lk=speye(n);
for i=1:ceil(k/2)
    Lk=L*Lk;
end


[S, C] = graphconncomp(L);
% project y
yy=y;
for i=1:S
    idx=C==i;
    y(idx)=y(idx)-mean(y(idx));
end
%component in the null space
yy=yy-y;
xx=yy;

% initialization
conv=0;
x=y;%zeros(n,1);
if ~flag
    kk=k;
else
    kk=k-1;
end
w=Times_S(D,L,kk,x);


%stepsize
eta=1/1.7266e4;
tol=1e-4;

% x=y;
% m=n;
% z=zeros(size(A,1),1);
% u=z; % u=Lx
iter = 1;
maxiter=1000;
while ~conv
    %solve O^T grad = (O^+ w -y) 
    grad =times_invST(D,L,kk,x-y,pfun); 

    %eta=0.01/iter;
    %proximal operator
    if flag %even
        w=graphtv(w-eta*grad,edges1,edges2,lambda);
    else
        w=soft_thresh(w-eta*grad,lambda); 
    end
%    w=w-mean(w);
    %dual certificate
    u=grad*min(1,lambda/max(abs(grad)));
    % solve for it v= O'\u
    v=times_invST(D,L,kk,-u,pfun);
    x=times_invS(D,L,kk,w,pfun);
    
    if ~flag%odd
        Pobj=0.5*norm(y-x)^2 + lambda*norm(w,1);
    else
        Pobj=0.5*norm(y-x)^2 + lambda*norm(D*w,1);
    end
    Dobj= -0.5*(v'*v) + v'*y; 
    dualgap=Pobj - Dobj;
    history=[history, dualgap];
    fprintf('[%d] Primal=%f, Dual=%f, gap = %f.\n',iter,Pobj,Dobj, dualgap);
    if dualgap < tol
        conv=1;
        fprintf('converged.\n')
    elseif iter >maxiter
        conv=1;
        fprintf('Reached maxiter.\n')
    end
    iter=iter+1;

end
%add the
x=x+xx;
end

function x = times_invS(D,L,k,u,pfun)
%pfun is the CMG preconditioner for L
%solve S x = u (assume u is in the range of S so it has a solution)
tol_pcg=1e-10;
    if ~mod(k,2)%even
        x=D'*u;
        for i=1:(k/2)
            [x,pcg_flg]=pcg(L,x,tol_pcg,100,pfun);
            x=x-mean(x);
        end
    else %odd
        x=u;
        for i=1:(k+1)/2
            [x,pcg_flg]=pcg(L,x,tol_pcg,100,pfun);
            x=x-mean(x);
        end
    end  
end

function x=times_invST(D,L,k,u,pfun)
%pfun is the CMG preconditioner for L
%solve S^T x = u (assume u is in the range of S so it has a solution)
tol_pcg=1e-10;
    if ~mod(k,2)%even
        a=u;
        for i=1:(k/2)+1
            [a,pcg_flg]=pcg(L,a,tol_pcg,100,pfun);
            a=a-mean(a);
        end
        x=D*a;
    else %odd (the same as before)
        x=u;
        for i=1:(k+1)/2
            [x,pcg_flg]=pcg(L,x,tol_pcg,100,pfun);
            x=x-mean(x);
        end
    end 
end

function [STu]=Times_ST(D,L,k,u)
%multiply S'
    if ~mod(k,2)%even
        STu=D'*u;
        for i=1:(k/2)
            STu=L*STu;
        end
    else %odd
        STu=u;
        for i=1:(k+1)/2
            STu=L*STu;
        end
    end  
end

function [Sx]=Times_S(D,L,k,x)
%multiply S
    Sx=x;
    if ~mod(k,2)%even
        for i=1:(k/2)
            Sx=L*Sx;
        end
        Sx=D*Sx;
    else %odd
        Sx=x;
        for i=1:(k+1)/2
            Sx=L*Sx;
        end
    end  
end