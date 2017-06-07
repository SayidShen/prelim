%% Generating Speed Comparisons
% 1:  For ProjNewton vs. ADMM on k=1
% 2:  For ADMM2  vs. ADMM 1 on k=2


%FLAG=1;% for exp 1
FLAG=0; % uncomment this for exp2

%% Test the admm on an realistic image
addpath('CMG')
%% Test script for general graph



img=imresize(imread('TVexact/lena.bmp'),0.25);
figure(1)
imshow(img);

img1=double(img)/255+0.1*randn(size(img));
figure(2)
imshow(img1);

%%

[ht,wt]=size(img1(:,:,1));
n=ht*wt;
m=(n*4-2*(ht+wt))/2;
edges1=zeros(m,1);
edges2=zeros(m,1);


%%
count=0;
idx = sub2ind([ht wt], kron(ones(wt-1,1),(1:ht)'), kron((1:(wt-1))',ones(ht,1)));
sz=length(idx);
count=count+sz;
edges1(1:sz)=idx;
%
idx = sub2ind([ht wt], kron(ones(wt-1,1),(1:ht)'), kron((2:(wt))',ones(ht,1)));
edges2(1:sz)=idx;

%
idx = sub2ind([ht wt], kron(ones(wt,1),(1:(ht-1))'), kron((1:wt)',ones(ht-1,1)));
sz=length(idx);
edges1(count+1:count+sz)=idx;
idx = sub2ind([ht wt], kron(ones(wt,1),(2:ht)'), kron((1:wt)',ones(ht-1,1)));
edges2(count+1:count+sz)=idx;

%%

D = sparse((1:m)',edges1,1,m,n,10*m);
D = D+ sparse((1:m)',edges2,-1,m,n);



%% 

if FLAG % run the proj newton exp
img2=img1;

for i=1%:3
lambda=0.2;
y=img1(:,:,i);
y=y(:);
%y=mu;
k=2;
%
if ~mod(k-1,2)
    rho=10*lambda;
else
    rho=5*lambda;
end
rho=lambda;
w=ones(m,1);
mask=rand(n,1)>0.9;
tic;
%[ x ,history1] = gtc_admm_v2( y,mask,[edges1,edges2],w,k-1,lambda,rho);
[ x ,history1] = gtf_admm_v2( y,[edges1,edges2],k-1,lambda,rho);
time_elapsed1=toc
tic;
%[ x ,history1] = gtf_apg( y,[edges1,edges2],k-1,lambda);
[ x2,history2 ] = gtf_proj_newton1( y,[edges1,edges2],k-1,lambda);
time_elapsed2=toc
img2(:,:,i)=reshape(x,[ht,wt]);
end



%%
h=figure(4)
hold off
step=time_elapsed1/length(history1);
step2=time_elapsed2/length(history2);
semilogy(step:step:time_elapsed1,0.5*(history1(1,:)+history1(2,:))','b-',...
    step2:step2:time_elapsed2,history2','r-','linewidth',1.5,'markersize',15)
l=legend('ADMM residual', 'ProjNewton duality gap');
set(l,'fontsize',14,'location','best')
%title('k=5 L^3 prox')
xlabel('Clock time in second','fontsize',14)
grid on

saveTightFigure(h,'PN_vs_ADMM_clocktime.pdf')

%%
h=figure(5)
hold off
step=time_elapsed1/length(history1);
step2=time_elapsed2/length(history2);
semilogy(0.5*(history1(1,:)+history1(2,:))','b-','linewidth',1.5)
hold on
semilogy(history2','r-','linewidth',1.5,'markersize',15)
l=legend('ADMM residual', 'ProjNewton duality gap');
set(l,'fontsize',14,'location','best')
%title('k=5 L^3 prox')
xlabel('Iterations','fontsize',14)
grid on

saveTightFigure(h,'PN_vs_ADMM_iteration.pdf')

else % run the admm vs admm exp
    
    
for i=1%:3
lambda=0.2;
y=img1(:,:,i);
y=y(:);
%y=mu;
k=3;
%
if ~mod(k-1,2)
    rho=10*lambda;
else
    rho=5*lambda;
end
rho=lambda;
w=ones(m,1);
mask=rand(n,1)>0.9;
tic;
%[ x ,history1] = gtc_admm_v2( y,mask,[edges1,edges2],w,k-1,lambda,rho);
[ x2 ,history2] = gtf_admm_v1( y,[edges1,edges2],k-1,lambda,rho);
time_elapsed2=toc
tic;
%[ x ,history1] = gtf_apg( y,[edges1,edges2],k-1,lambda);
[ x,history1 ] = gtf_admm_v2( y,[edges1,edges2],k-1,lambda,rho);
time_elapsed1=toc
img2(:,:,i)=reshape(x,[ht,wt]);
end

    
    h=figure(4)
hold off
step=time_elapsed1/length(history1);
step2=time_elapsed2/length(history2);
semilogy(step:step:time_elapsed1,0.5*(history1(1,:)+history1(2,:))','b-',...
    step2:step2:time_elapsed2,0.5*(history2(1,:)+history2(2,:))','c-','linewidth',1.5,'markersize',15)
l=legend('Special ADMM residual', 'Naive ADMM residual');
set(l,'fontsize',14,'location','best')
%title('k=5 L^3 prox')
xlabel('Clock time in second','fontsize',14)
grid on

saveTightFigure(h,'ADMM_vs_ADMM_clocktime.pdf')

h=figure(5)
hold off
step=time_elapsed1/length(history1);
step2=time_elapsed2/length(history2);
semilogy(0.5*(history1(1,:)+history1(2,:))','b-','linewidth',1.5)
hold on
semilogy(0.5*(history2(1,:)+history2(2,:))','c-','linewidth',1.5,'markersize',15)
l=legend('Special ADMM residual', 'Naive ADMM residual');
set(l,'fontsize',14,'location','best')
%title('k=5 L^3 prox')
xlabel('Iterations','fontsize',14)
grid on

saveTightFigure(h,'ADMM_vs_ADMM_iteration.pdf')
    
end


