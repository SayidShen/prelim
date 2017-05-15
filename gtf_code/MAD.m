function [ Yhat ] = MAD( Y, D, w, lambda, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin>4
    WEIGHTED=1;
    varTrain=varargin{1};
else
    WEIGHTED=0;
end

[n,C]=size(Y);
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

mask=(Y~=-1);
L=D'*diag(w)*D;

Yhat=zeros(size(Y));
tol_pcg=1e-8;
for i=1:C
    if ~WEIGHTED
        I=mask(:,i);
    else
        I=mask(:,i)./varTrain;
    end
    %xm=(lambda*L+speye(sum(I)))\(Y(mask,i));
    Yhat(:,i) = pcg(lambda*L+diag(I),Y(:,i).*I,tol_pcg,200);
end

end

