function [ x, history] = sgtf_admm( y,D,k,lambda1, lambda2,rho1,rho2, varargin )
addpath('TVexact')
%SGTF_ADMM Sparse Graph Trend Filtering
%   Detailed explanation goes here
% min 0.5||x-y||^2  + lambda1*|O*x|_1 + lambda2*|x|_!

%GTF_ADMM Summary of this function goes here
%   ADMM version of the graph trend filtering
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


if k==0 % direct solution
    fprintf('Direct Solution.\n');
    x=graphtv(y,edges1,edges2,lambda1);
    x=soft_thresh(x,lambda2);
    history=0;
    return;
end

history=[];

L=D'*D;
Lk=speye(n);
for i=1:ceil(k/2)
    Lk=L*Lk;
end

%preprocessing for SDD solver
A=(rho1*(Lk*Lk)+(1+rho2)*speye(n));
tol_pcg=1e-6;
pfun=cmg_sdd(A);


tol_abs=1e-5;
tol_rel=1e-4;


% initialization
conv=0;
x=y;
m=n;
z=zeros(size(A,1),1);
u=z; % u=Lx
w=zeros(size(x,1),1);
v=w;

iter = 1;
maxiter=1000;
while ~conv
    %[x,pcg_flg]=pcg(A,(Lk*(rho1*z-u)+(rho2*w-v)+y),tol_pcg,100,pfun);
    [x,pcg_flg]=pcg(A,(Lk*(rho1*z-u)+(rho2*w-v)+y),tol_pcg,100);
%     if pcg_flg~=0
%         fprintf('Warning: PCG did not converge.\n');
%     end
    Lkx=Lk*x;
    
    if flag
        z_new=graphtv(Lkx+u/rho1,edges1,edges2,lambda1/rho1);
    else
        z_new=soft_thresh(Lkx+u/rho1,lambda1/rho1);
    end

    w_new=soft_thresh(x+v/rho2,lambda2/rho2);
    
    s= norm([rho1*Lk*(z_new-z);rho2*(w_new-w)]);
    z=z_new;
    w=w_new;
    
    u=u+rho1*(Lkx-z);
    v=v+rho2*(x-w);
    r=norm([Lkx-z;x-w]);
    
    eps_pri= sqrt(n)*tol_abs + tol_rel*max(norm(Lkx),norm(z));    
    eps_dual= sqrt(m)*tol_abs  + tol_rel*norm(L*u);
   
    %fprintf('[%d] [r, eps_pri]=%f, %f, [s, eps_dual]=%f, %f. \n',...
    %    iter,r,eps_pri,s, eps_dual);
    if ~mod(iter,50)
        fprintf('[%d] [r,s]=%f, %f, [eps_pri, eps_dual]=%f, %f. \n',...
            iter,r,s, eps_pri,eps_dual);
    end
    if r < eps_pri && s< eps_dual
        conv=1;
        fprintf('converged.\n')
    elseif iter >= maxiter
        conv=1;
        fprintf('Reached maxiter.\n')
    end
    
    history=[history,[s;r]];
    
    iter=iter+1;
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
