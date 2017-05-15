function [ x, history] = gtf_admm_grid( y,sz,D,k,lambda, rho, varargin )
addpath('TVexact')
%GTF_ADMM_grid Summary of this function goes here

%   ADMM version of the graph trend filtering
% for specialized graphs (2D grid)
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

n1=sz(1);
n2=sz(2);


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

% Else just use ADMM
history=[];

L=D'*D;
Lk=speye(n);
for i=1:ceil(k/2)
    Lk=L*Lk;
end



tol_abs=1e-5;
tol_rel=1e-4;


% initialization
conv=0;
x=y;
m=n;
z=zeros(n,1);
u=z; % u=Lx
iter = 1;
maxiter=1000;
while ~conv

   %[x,pcg_flg]=pcg(A,(Lk*(rho*z-u)+y),tol_pcg,1000,pfun);
   %[x,pcg_flg]=pcg(A,(Lk*(rho*z-u)+y),tol_pcg,k*100);

      x = grid_system(reshape((Lk*(rho*z-u)+y),[n1,n2]),ceil(k/2)*2,rho);
   x=x(:);
   
    Lkx=Lk*x;
    
    if flag
        z_new=graphtv(Lkx+u/rho,edges1,edges2,lambda/rho);
    else
        z_new=soft_thresh(Lkx+u/rho,lambda/rho);
    end
    s= rho*norm(Lk*(z_new-z));
    z=z_new;
    
    u=u+rho*(Lkx-z);
    r=norm(Lkx-z);
    
    eps_pri= sqrt(n)*tol_abs + tol_rel*max(norm(Lkx),norm(z));    
    eps_dual= sqrt(m)*tol_abs  + tol_rel*norm(Lk'*u);
   
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
    
    
    tau=2;
    if r>10*s
        rho=rho*tau;

    elseif s>10*r
        rho=rho/2;

    end
    
    
    
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
