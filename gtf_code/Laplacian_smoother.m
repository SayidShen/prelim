function [x]=Laplacian_smoother(y,D,lambda,k,normalized)
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

if normalized
    Deg=sum(abs(D))/2;
    D=D*spdiags((Deg.^(-1/2))',0,n,n);
end

L=D'*D;
O=L;
for i=1:k
    O=O*L;
end

       [LF,p,S] = chol(speye(n)+lambda*O,'lower','matrix');
            if p>0
                fprintf('Warning! Not PSD');
            end
            x=S*(LF'\(LF\(S'*y)));
 
    %x=(speye(n)+lambda*O)\y;

    
end