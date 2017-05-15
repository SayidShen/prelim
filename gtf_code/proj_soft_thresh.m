function [ x ] = proj_soft_thresh( y,lambda )
%UNTITLED Summary of this function goes here
%   Algorithmically solve the "mean-0" constrained l1 prox operator.
%   assume input y obeys sum(y)=0;
n=length(y);
[ysort, idx]=sort(y);%ascending order
l=0;r=n;
lbound=-inf; rbound=ysort(1);
s=0;

fmax=0;
vmax=0;


v=(l-r-sumS)/(n-s);
if lbound<v && v<rbound
    f=fobj(v,lambda,l,r,s,n,sumS,sumR,sumL,sumS2);
else
    fl=fobj(lbound,lambda,l,r,s,n,sumS,sumR,sumL,sumS2);
    fr=fobj(rbound,lambda,l,r,s,n,sumS,sumR,sumL,sumS2);
    if fl>fr
        f=fl; v=fl;
    else
        f=fr;v=fr;
    end
end

if f>fmax
    % update the current solution
    fmax=f;
    vmax=v;
    x=soft_thresh(y-v,lambda);
end






end

function f=fobj(v,lambda,l,r,s,n,sumS,sumR,sumL,sumS2)
    f = -(l+r)*lambda^2/2 + (l-r-sumS)*v...
        + 0.5*(s/2-n/2)*v^2 + 0.5*sumS2+lambda*(sumR-sumL);
end