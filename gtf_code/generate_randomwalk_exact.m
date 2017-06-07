%% Generate facebook data by random walk

load('../facebook.txt')
addpath('RandomWalks')

facebook=facebook+1;%change to 1 based indexing
n=max(facebook(:));
[ D, edges ] = preprocess_D(facebook,n);

%% generate underlying signals
rng(15213)

headnodes=[0,107,1684,1912,3437,348,3980,414,686,698];
headnodes=headnodes+1;
%s=length(headnodes);
%idx=headnodes;

s=10;
idx=randperm(n);
%idx=[idx(1:s), headnodes(1),headnodes(3)];
s=length(idx);


s=20;
idx=[idx(1:s), headnodes];
s=length(idx);


L=D'*D; A=-L;
A = spdiags(zeros(n,1),0,A); 


Deg_sqrtinv=spdiags(1./sqrt(sum(A))',0,n,n);
[Q,Lambda]=eig(full(Deg_sqrtinv*A*Deg_sqrtinv));
tau=0.25;%decay rate
lambdas=diag(Lambda);
ee=zeros(n,1);ee(idx)=1;
y=Deg_sqrtinv*( Q*(diag(1./(1-lambdas.*tau))*(Q'*(spdiags(sqrt(sum(A))',0,n,n)*ee))));


%% (uncomment to) Generate inhomogeneous Random Walk data

tau_list=rand(10,1);
num_list=rand(10,1)*1000;

y=zeros(n,1);
for i=1:10
ee = sparse(idx(i),1,num_list(i),n,1);
y=y+ Deg_sqrtinv*( Q*(diag(1./(1-lambdas.*tau_list(i)))*(Q'*(spdiags(sqrt(sum(A))',0,n,n)*ee))));
end

%save('Facebook_data_inhomogeneous2.mat','y')

%% (uncomment to) Generating homogeneous Random Walk data

tau=0.95;
tau_list=ones(10,1)*0.7;%95;
num_list=ones(10,1)*500;

y=zeros(n,1);
for i=1:10
ee = sparse(idx(i),1,num_list(i),n,1);
y=y+ Deg_sqrtinv*( Q*(diag(1./(1-lambdas.*tau))*(Q'*(spdiags(sqrt(sum(A))',0,n,n)*ee))));
end
y=Deg_sqrtinv*( Q*(diag(1./(1-lambdas.*tau_list(i)))*(Q'*(spdiags(sqrt(sum(A))',0,n,n)*ee))));


%save('Facebook_data_inhomogeneous2.mat','y')
save('Facebook_data_homogeneous.mat','y')
