function [ x,history ] = gtf_proj_newton1( y,D,k,lambda)
%GTF_ADMM Summary of this function goes here
%   ADMM version of the graph trend filtering
% lambda is tuning parameter, rho is the numerical paramter
% D is the incidence matrix with dimension edge times node
% k is the order.
% if nargin == 6
%     idx=varargin{1};
% else
%     i=(1-mod(k,2))+1;
%     idx=false(size(D,i),1);
% end
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



L=D'*D;

% construct operator O 
D=sparse(D);
O=D;
for i=1:k
    if mod(i,2)
        O= D'*O;
    else
        O=D*O;
    end
end

% form it explicitly now... for test of concept
H=O*O';
H=H;

diagH=spdiags(H,0);

%pre-compute the combinatorial preconditioner
%v=\delta-u
tol_pcg=1e-5;
%pfun=cmg_sdd(L);


kmax=n-1;%ceil(n/2);

%line-search parameter
alpha=0.1; %intepolation parameter
beta=0.8; %backtracking parameter


%tolerance and maxiter
tol=max(1e-4,tol_pcg);
maxiter=30;
%tol_pcg=1e-5;
history=[];

% initialization
conv=0;
%z=zeros(size(O,1),1);
u=zeros(size(O,1),1);
fx=0;% when u=0
iter = 1;
deltaIXold=randn(size(u));

while ~conv   
    %
    grad=compute_grad(D,L,k,y,u);
    
    [IX]=free_set(u,grad,lambda);
    idx=find(IX>0);
    sz=sum(IX>0);
    %idx=randperm(sz);
    %if sz>=kmax
    %    IX(idx(kmax+1:end))=0;
    %end
    
    
    delta=zeros(size(u));
    
    
    %delta(IX) = lsqr(H(IX,IX),grad(IX));
    %opts.display=0;
    %[pfun, H_pcd]=cmg_sdd(H(IX,IX),opts);
    
    % use "@(x)Times_H_IX(D,L,k,IX,x)" to replace H(IX,IX) when k is large.
    
    %if k==1
%        [pfun, H_pcd]=cmg_sdd(H(IX,IX));
 %       [delta(IX), pcgflag]= pcg(H(IX,IX),grad(IX),tol_pcg,200,pfun);
        [delta(IX), pcgflag, relres] = pcg(H(IX,IX),grad(IX),tol_pcg,200);
%         if pcgflag~=0 && relres>0.1
%             fprintf('Warning: PCG did not converge. Start using Cholesky\n');
%             %CHOLESKY=1;
%             [LF,p,S] = chol(H(IX,IX),'lower','matrix');
%             if p>0
%                 fprintf('Warning! Not PSD');
%             end
%             delta(IX)=S*(LF'\(LF\(S'*grad(IX))));
%             if norm(delta-deltaIXold)<10*eps
%                 fprintf('Stagnant iterations. quiting.\n')
%                 conv=1;
%             end
%         end
        deltaIXold=delta;
%        [delta(IX), pcgflag] = pcg(H(IX,IX),grad(IX),tol_pcg,200,spdiags(1./diagH,0,sz,sz));
    %else
    %    [delta(IX), pcgflag] = pcg(H(IX,IX),grad(IX),tol_pcg,200);
    %end
    %delta(IX) = pinv(full(H(IX,IX)))*grad(IX);
    
    %delta=v+u-mean(u);
    %Projected Newton step with backtracking line search
    gamma=1;%stepsize
    %u_new=proj_lambda(u-gamma*(v+u),lambda);
    u_new=proj_lambda(u-gamma*delta,lambda);
    fx_new=obj_eval(D,L,k,y,u_new); 
     while  fx_new > fx - alpha*gamma*grad'*(delta)
         gamma=beta*gamma;
         u_new=proj_lambda(u-gamma*delta,lambda);
         fx_new=obj_eval(D,L,k,y,u_new);
     end
     
     %update parameters
     fx=fx_new;
     u=u_new;
    iter=iter+1;
    STu=Times_ST(D,L,k,u);
    x=-STu+y;
    dgap=lambda*norm(Times_S(D,L,k,x),1)-STu'*x;
    fprintf('iter = %d. duality_gap= %f\n',iter,dgap);    
    history=[history, dgap];
    if dgap<tol
        conv=1;
        fprintf('converged.\n')
    elseif iter>maxiter
        conv=1;
        fprintf('reach maxiter.\n')
    end
end


end

function [IX]=free_set(u,grad,lambda)
    IX=~((grad>10*eps & u==-lambda) | (grad<-10*eps & u==lambda));
end


function [u_proj]=proj_lambda(u,lambda)
u_proj=min(u,lambda);
u_proj=max(u_proj,-lambda);
end

function [fx]=obj_eval(D,L,k,y,u)  
    STu=Times_ST(D,L,k,u); 
    fx=0.5*(STu'*STu)-STu'*y;
end

function [grad]=compute_grad(D,L,k,y,u)
    grad=Times_S(D,L,k,Times_ST(D,L,k,u)-y);
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

function [H_IX_x]=Times_H_IX(D,L,k,IX,x)
    if ~mod(k,2)  %even
        y=D(IX,:)'*x;   
        for i=1:k
            y=L*y;
        end
        H_IX_x=D(IX,:)*y;
    else % odd
        y=L(:,IX)*x;
        for i=1:k-1
            y=L*y;
        end
        H_IX_x=L(IX,:)*y;
    end
    
end
