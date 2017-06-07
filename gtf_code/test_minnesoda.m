%% Simulating graph signal on real geographic graph


% Load a road network to use for statistical computations
load_gaimc_graph('minnesota');
figure(1)
gplot(A,xy);


%%

% implement A to edge1,edge2 conversion

idx=xy(:,1)<-92.5 & xy(:,1)>-94 & xy(:,2)<45.3 & xy(:,2)>44.4;
hold on;
plot(xy(idx,1),xy(idx,2),'rx','markersize',10);

n=size(A,1);
% for i=1:n
%     gplot(A,xy);
% 
%     hold on;
%     
%     fprintf('%d th node\n',i);
%     hold off;
%     pause;
% end
    
%% solve the problem
[edges1,edges2,val]=find(triu(A));
y=zeros(size(A,1),1);

y(idx)=2; y(~idx)=1;
y1=y+randn(n,1);
yrange=[0,3];
gplot_value(A,xy,y1, yrange);

k=3;lambda=1;rho=1;
[ x ,history1] = gtf_admm_v2( y,[edges1,edges2],k-1,lambda,rho);

%%
figure(2)
gplot_value(A,xy,x, yrange);


%%
figure(3)
mesh(xy(:,1),xy(:,2),x);

figure(4)
mesh(xy(:,1),xy(:,2),y);
