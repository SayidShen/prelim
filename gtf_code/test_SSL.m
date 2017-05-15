%% data generation
addpath('Spectral Clustering')
data=load('jain.txt');%spiral.txt, jain.txt,  R15.txt
true_labels=data(:,end);
datapoints=data(:,1:end-1);
[n,p]=size(datapoints);

affinity = CalculateAffinity(datapoints);
C=max(true_labels);

cc=colormap;
sz=length(colormap);
figure(1)
hold off
for i=1:C
    plot(datapoints(true_labels==i,1),...
        datapoints(true_labels==i,2),'.','Color',cc(floor(i*sz/C),:));
    hold on
end
hold off;

%%
y=zeros(n,C);
for i=1:C
    y(:,i)=true_labels==i;
end
mask=rand(n,1)>0.95;
Y=y;
Y(~mask,:)=-1;

K=25;

A=affinity;
A(logical(eye(n))) = 0;
for i=1:n
    [a,idx]=sort(A(:,i),'descend');
    A(idx(K+1:end),i)=0;
end

[edges1,edges2, w]=find(A);
A=A'+A;
%%

%w=ones(length(edges1),1);% we can either use weight or not
lambda=0.01;
[ Yhat ] = MAD( Y, [edges1,edges2], w , lambda );
k=0;
[ Yhat1 ] = MAD_gtf( Y, [edges1,edges2], w, lambda,k );
[tmp, est_labels]= max(Yhat,[],2);
[tmp, est_labels1]= max(Yhat1,[],2);


%(Yhat(:,1)>Yhat(:,2))+1;
%%
figure(1)
hold off
for i=1:C
    plot(datapoints(true_labels==i,1),...
        datapoints(true_labels==i,2),'.','Color',[0.5,0.5,0.5]);
    hold on
end

for i=1:C
    plot(datapoints(true_labels==i&mask,1),...
        datapoints(true_labels==i&mask,2),'o','Color',cc(floor(i*sz/C),:),'linewidth',2);
    hold on
end

hold off;

figure(2)
hold off
for i=1:C
    plot(datapoints(est_labels==i,1),...
        datapoints(est_labels==i,2),'.','Color',cc(floor(i*sz/C),:));
    hold on
end

for i=1:C
    plot(datapoints(est_labels==i&mask,1),...
        datapoints(est_labels==i&mask,2),'o','Color',cc(floor(i*sz/C),:),'linewidth',2);
    hold on
end

hold off;

%%
Yhat=bsxfun(@rdivide,Yhat,sum(Yhat,2));
Yhat1=bsxfun(@rdivide,Yhat1,sum(Yhat1,2));
figure(3)
imagesc(Yhat);
figure(4)
imagesc(Yhat1);