function [ x ] = gtc_CVX( y,mask,D,w,k,lambda, rho, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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

L=D'*D;
Lk=speye(n);
for i=1:ceil(k/2)
    Lk=L*Lk;
end
if ~mod(k,2)%even
    O=D*Lk;
else
    O=Lk;
end


cvx_begin
    variable x(n)
    minimize( 0.5*((x-y).*sqrt(mask))'*((x-y).*sqrt(mask)) +  lambda*norm(O*x))
cvx_end





end

