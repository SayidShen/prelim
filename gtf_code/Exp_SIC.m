addpath('SIC2004')
addpath('Spectral Clustering')
addpath('gridfitdir')

X=csvread('SIC2004_input.csv');
X=csvread('SIC2004_joker.csv');
Y=csvread('SIC2004_out.csv');
Z=csvread('1st_file_true_values.csv');
Z=csvread('2nd_file_true_values.csv');




%%

y0=Z(:,4);

xy1=X(:,2:3);
xy2=Y(:,2:3);
figure(4);
plot(xy1(:,1),xy1(:,2),'ob',...
    xy2(:,1),xy2(:,2),'xr')

%% preprocessing. 
xy=[xy1;xy2];

n=length(xy1)+length(xy2);
mask=false(n,1);mask(1:length(xy1))=1;
y=zeros(n,1);
y(mask)=X(:,4);
y(~mask)=-1;

%also construct the ground truth
y00=y;
y00(~mask)=y0;



%% Variance estimate
TrainData=zeros(n,1);
for i=1:9
    T=csvread(['sic2004_0',int2str(i),'.csv']);
    TrainData(mask,i)=T(:,4);
end
T=csvread('sic2004_10.csv');
TrainData(mask,10)=T(:,4);

varTrain=var(TrainData,1,2);
varTrain(~mask)=-1;

%% getting graphs.
affinity = CalculateAffinity(xy/10000,1);
A=affinity;
A(logical(eye(n))) = 0;
K=5;
for i=1:n
    [a,idx]=sort(A(:,i),'descend');
    A(idx(K+1:end),i)=0;
end

[edges1,edges2, w]=find(A);
A=A'+A;



%%
%varTrain=ones(n,1);

lambda_list=exp(-10:0.5:10); k=1;
sz=length(lambda_list);
mse=zeros(sz,2);

for i=1:length(lambda_list)
    lambda=lambda_list(i);
    [ yhat ] = MAD( y, [edges1,edges2], w , lambda, varTrain );
    [ yhat1 ] = MAD_gtf( y, [edges1,edges2], w, lambda,k, 1,varTrain);

    mse(i,:)=[norm(yhat(~mask)-y0),norm(yhat1(~mask)-y0)];
    mae(i,:)=[norm(yhat(~mask)-y0,1),norm(yhat1(~mask)-y0,1)];
end

%%
%
[minmse, ii]=min(mse(:,2));
%ii=13
lambda=lambda_list(ii);
[ yhat1 ] = MAD_gtf( y, [edges1,edges2], w, lambda,k, 1,varTrain);

[minmse, ii]=min(mse(:,1));
%ii=11
lambda=lambda_list(ii);
[ yhat ] = MAD( y, [edges1,edges2], w, lambda,varTrain);

%%
xnodes=min(xy(:),1):(max(xy(:,1))-min(xy(:,1)))/50:max(xy(:,1));
ynodes=min(xy(:),1):(max(xy(:,2))-min(xy(:,1)))/50:max(xy(:,2));
[xi,yi] = meshgrid(xnodes, ynodes);

F = scatteredInterpolant(xy,yhat1);
F.Method = 'nearest';
zi=F(xi,yi);



%xnodes=min(xy(:),1):(max(xy(:,1))-min(xy(:,1)))/50:max(xy(:,1));
%ynodes=min(xy(:),1):(max(xy(:,2))-min(xy(:,1)))/50:max(xy(:,2));
%[xi,yi] = meshgrid(xnodes, ynodes);
figure(3)
hold off

%zi = griddata(xy(:,1),xy(:,2),yhat1,xi,yi);
surf(xi,yi,zi);

%[c,h] = contour(xi,yi,zi);
%clabel(c,h);
%plot3(xy(mask,1),xy(mask,2),yhat1(mask),'b.',...,
%    xy(~mask,1),xy(~mask,2),yhat1(~mask),'bo')
xlabel('Longitude')
ylabel('Latitude')
zlabel('Radiation Level')
%
hold on
%plot3(xy(:,1),xy(:,2),yhat1,'.')
zlim([0,1600])


figure(2)
hold off

F = scatteredInterpolant(xy,yhat);
F.Method = 'nearest';
zi=F(xi,yi);

surf(xi,yi,zi);
%[c,h] = contour(xi,yi,zi);
%clabel(c,h);
xlabel('Longitude')
ylabel('Latitude')
zlabel('Radiation Level')
hold on
%plot3(xy(:,1),xy(:,2),yhat,'.')

zlim([0,1600])

figure(1)
hold off

F = scatteredInterpolant(xy,y00);
F.Method = 'nearest';
zi=F(xi,yi);
surf(xi,yi,zi);
%[c,h] = contour(xi,yi,zi);
%clabel(c,h);
xlabel('Longitude')
ylabel('Latitude')
zlabel('Radiation Level')
hold on
%plot3(xy(:,1),xy(:,2),y00,'.')

zlim([0,1600])

%%
figure(4)
semilogx(lambda_list,mse)
legend('Laplacian Interpolation','Graph Trend Interpolation')
