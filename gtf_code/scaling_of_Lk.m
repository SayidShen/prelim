x=zeros(10,1);
x2=zeros(10,1);
Lk=speye(n);
for i=1:10
    Lk=Lk*L;
    x(i)=norm(Lk,'fro');
    s=svds(Lk);
    x2(i) =s(1); 
end

%%
y=zeros(n,1);
y(101:105)=randn(5,1);
y=y-mean(y);
yy=pcg(L,y);
yy=yy/norm(yy);

y1=zeros(n,1);
y1(101:200)=randn(100,1);
y1=y1-mean(y1);
yy1=pcg(L,y1);
yy1=yy1/norm(yy1);

y2=zeros(n,1);
y2=randn(n,1);
y2=y2-mean(y2);
yy2=pcg(L,y2);
yy2=yy2/norm(yy2);

%%
hold off
semilogy(1:10,x,'-bo',...
    1:10,x2,'-rx',...
    1:10,fs,'--g',...
    1:10,fs1,'-g',...
1:10,fs2,'--k')
legend('Frobenius norm','spectral norm','5-sparse-l1','100-sparse-l1','dnese-l1');






fs=zeros(10,1);
fs1=zeros(10,1);
fs2=zeros(10,1);

Lk=speye(n);
for i=1:10
    Lk=Lk*L;
    fs(i)=norm(Lk*yy,1);
     fs1(i)=norm(Lk*yy1,1);
     fs2(i)=norm(Lk*yy2,1);
end