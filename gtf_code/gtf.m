function [ x, history] = gtf( y,D,k,lambda, rho, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
    fprintf('k=0: Direct Solution.\n');
    x=graphtv(y,edges1,edges2,lambda);
    history=0;
    return;
elseif mod(k,2) % odd (assume sparse graph)
    fprintf('k=1: Run projected newton.\n')
    [ x,history ] = gtf_proj_newton1( y,[edges1,edges2],k,lambda);
    return;
else
    fprintf('others: Run ADMM.\n')
    [ x ,history] = gtf_admm_v2( y,[edges1,edges2],k,lambda,rho);
    return;
end



