function [ x, history] = gtc_admm_v2( y,mask,D,w,k,lambda, rho, varargin )
addpath('TVexact')
%GTC_ADMM_V2 Summary of this function goes here
%   Detailed explanation goes here
%   Graph trend completion. with observed data in mask.

%GTC_ADMM_V2 Summary of this function goes here
%   ADMM version of the semi-supervised graph trend filtering
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
D=spdiags(sqrt(w),0,m,m)*D;


% Choose the right algorithm
% if k==0 % direct solution
%     fprintf('k=0: Direct Solution.\n');
%     x=graphtv(y,edges1,edges2,lambda,w);
%     history=0;
%     return;
% end

% Else just use ADMM
history=[];

L=D'*D;
Lk=speye(n);
for i=1:ceil(k/2)
    Lk=L*Lk;
end

%preprocessing for SDD solver
A=(rho*(Lk*Lk)+spdiags(mask,0,n,n));
%LF = chol(A,'lower');
tol_pcg=1e-6;
%[pfun, H]=cmg_sdd(A);
CHOLESKY=0;



tol_abs=1e-5;
tol_rel=1e-4;


% initialization
conv=0;
x=y;
m=n;
z=zeros(size(A,1),1);
u=z; % u=Lx
iter = 1;
maxiter=1000;
while ~conv

   %[x,pcg_flg]=pcg(A,(Lk*(rho*z-u)+y),tol_pcg,1000,pfun);
   %[x,pcg_flg]=pcg(A,(Lk*(rho*z-u)+mask.*y),tol_pcg,k*200+100);
    %if pcg_flg~=0
    %    fprintf('Warning: PCG did not converge.\n');
    %end
   if ~CHOLESKY
   [x,pcg_flg]=pcg(A,(Lk*(rho*z-u)+mask.*y),tol_pcg,k*100,spdiags(diag(A),0,size(A,1),size(A,2)));
    if pcg_flg~=0
        fprintf('Warning: PCG did not converge. Start using Cholesky\n');
        CHOLESKY=1;
        [LF,p,S] = chol(A,'lower','matrix');
    end
   end
   if CHOLESKY   
    xx=Lk*(rho*z-u)+mask.*y;
    x=S*(LF'\(LF\(S'*xx)));
   end
    
    

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
%         tau=2;
%     if r>100*s
%         rho=rho*tau;
%         A=A*tau-(tau-1)*speye(n);
%         if CHOLESKY
%             [LF,p,S] = chol(A,'lower','matrix');
%         end
%     elseif s>100*r
%         rho=rho/2;
%         A=A/tau+(1-1/tau)*speye(n);%(rho*(Lk*Lk)+speye(n));
%         if CHOLESKY
%             [LF,p,S] = chol(A,'lower','matrix');
%         end
%     end
%     
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
