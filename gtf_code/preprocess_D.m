function [ D, edges ] = preprocess_D(D,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%   Input:   D: can either be edgelist  or incidence matrix
%            n: number of nodes
%   Output:  D: Sparse incidence matrix
%            edges: Edgelist
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
edges=[edges1,edges2];

end

