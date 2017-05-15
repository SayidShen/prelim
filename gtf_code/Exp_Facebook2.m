%% Run the Chicago crime data set

load('../facebook.txt')

facebook=facebook+1;%change to 1 based indexing
n=max(facebook(:));
[ D, edges ] = preprocess_D(facebook,n);

%% generate underlying signals
rng(15213)

headnodes=[0,107,1684,1912,3437,348,3980,414,686,698];
headnodes=headnodes+1;
%s=length(headnodes);
%idx=headnodes;

s=20;
idx=randperm(n);
idx=[idx(1:s), headnodes];
s=length(idx);

%0 based k
k=1;%O should be L
O=construct_O( D,k );

ss=sparse(n,1);
ss(idx)=randn(s,1)*100;
ss(idx)=ss(idx)-mean(ss(idx));

tol_pcg=1e-6;
pfun=cmg_sdd(O);
y=pcg(O,ss,tol_pcg,100,pfun);

% %dense vector
% dd=randn(n,1);
% dd=dd-mean(dd);
% dd=dd/norm(dd)*norm(ss(idx));
% 
% y=pcg(O,dd,tol_pcg,100,pfun);


%% add noise 
rng(43651243)
sigma_list=[0.01,0.1,0.3,0.5,0.7,0.9];
%sigma_list=[0.01, 0.05];
szz=length(sigma_list);
mse_list=zeros(szz,1);

base=norm(D,1);

%% k=0
for i=1:szz
    sigma=sigma_list(i);
    yhat=y+sigma*randn(size(y));
    
    %%
    lambdalist=-6:0.2:6;
    %lambdalist=2:0.5:3;%-6:0.2:6;
    lambdalist=10.^lambdalist;
    %lambdalist=0.02;
    sz=length(lambdalist);
    mse_list_inner=zeros(sz,1);
    reg_weight=zeros(sz,1);
    for j=1:sz
        m=size(D,1);
        scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
        lambda=scalefactor*lambdalist(j); rho=10*lambda;
        
        [ x ,history1] = gtf_admm_v2( yhat,edges,k-1,lambda,rho);
        %[ x,history ] = gtf_proj_newton1( yhat,edges,k,lambda);

        mse_list_inner(j)=norm(x-y)^2/n;
        reg_weight(j)=norm(O*x);
        

 
        
    end
            loglog(lambdalist,mse_list_inner);
        hold on;
    mse_list(i)=min(mse_list_inner);
           pause(1);
    %save('temp.mat');
end

%% To compare against the Laplacian smoother 
mse_listL=zeros(szz,1);

rng(43651243)
for i=1:szz
        sigma=sigma_list(i);
    yhat=y+sigma*randn(size(y));
        lambdalist=-6:0.2:6;
    lambdalist=10.^lambdalist;
    sz=length(lambdalist);
    mse_list_inner=zeros(sz,1);
    reg_weight=zeros(sz,1);
    
    for j=1:sz
               
        scalefactor=base/norm(O,1);%base/norm(O,'fro')^2;
        lambda=lambdalist(j)*scalefactor;
        %lambda=exp(j/2*log(lambdalist(j)));
        x=Laplacian_smoother(yhat,edges,lambda,0,0);
        mse_list_inner(j)=norm(x-y)^2/n;
    end
                loglog(lambdalist,mse_list_inner);
        hold on;
    mse_listL(i)=min(mse_list_inner);
               pause(1);
end


% To compare against James's Wavelets
mse_listW=zeros(szz,1);
load FacebookWavletBasis.mat %got matrix B


rng(43651243)
for i=1:szz
        sigma=sigma_list(i);
    yhat=y+sigma*randn(size(y));
        lambdalist=-6:0.2:6;
    lambdalist=10.^lambdalist;
    sz=length(lambdalist);
    mse_list_inner=zeros(sz,1);
    reg_weight=zeros(sz,1);
    
    for j=1:sz
               
        scalefactor=base/norm(O,1);%base/norm(O,'fro')^2;
        lambda=lambdalist(j)*scalefactor;
        
             yrot=B*yhat;
     yrot=max(0,yrot-lambda*sqrt(2*log(n)))...%soft thresholding
         + min(0,yrot+lambda*sqrt(2*log(n))); 
     %)(abs(yrot)<lambda*sqrt(2*log(n))); 
     x=B'*yrot;
        
        mse_list_inner(j)=norm(x-y)^2/n;
    end
                loglog(lambdalist,mse_list_inner);
        hold on;
    mse_listW(i)=min(mse_list_inner);
               pause(1);
end


%% k=1
mse_list1=zeros(szz,1); k=1;
for i=1:szz
    sigma=sigma_list(i);
    yhat=y+sigma*randn(size(y));
    
    %%
    lambdalist=-6:0.5:3;
    %lambdalist=2:0.5:3;%-6:0.2:6;
    lambdalist=10.^lambdalist;
    %lambdalist=0.02;
    sz=length(lambdalist);
    mse_list_inner=zeros(sz,1);
    reg_weight=zeros(sz,1);
    for j=1:sz
        m=size(D,1);
        scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
        lambda=scalefactor*lambdalist(j); rho=10*lambda;
        
        [ x ,history1] = gtf( yhat,edges,k,lambda,rho);
        %[ x,history ] = gtf_proj_newton1( yhat,edges,k,lambda);

        mse_list_inner(j)=norm(x-y)^2/n;
        reg_weight(j)=norm(O*x);
        

 
        
    end
            loglog(lambdalist,mse_list_inner);
        hold on;
    mse_list1(i)=min(mse_list_inner);
           pause(1);
    %save('temp.mat');
end

%% k=2
mse_list2=zeros(szz,1); k=2;
for i=1:szz
    sigma=sigma_list(i);
    yhat=y+sigma*randn(size(y));
    
    %%
    lambdalist=-6:0.5:6;
    %lambdalist=2:0.5:3;%-6:0.2:6;
    lambdalist=10.^lambdalist;
    %lambdalist=0.02;
    sz=length(lambdalist);
    mse_list_inner=100*ones(sz,1);
    reg_weight=zeros(sz,1);
    for j=1:sz
        m=size(D,1);
        scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
        lambda=scalefactor*lambdalist(j); rho=0.5*lambda;
        
        [ x ,history1] = gtf( yhat,edges,k,lambda,rho);
        %[ x,history ] = gtf_proj_newton1( yhat,edges,k,lambda);

        mse_list_inner(j)=norm(x-y)^2/n;
        reg_weight(j)=norm(O*x);
        

 
        
    end
            loglog(lambdalist,mse_list_inner);
        hold on;
    mse_list2(i)=min(mse_list_inner);
           pause(1);
    %save('temp.mat');
end

% k=3
% mse_list3=zeros(szz,1); k=3;
% for i=1:szz
%     sigma=sigma_list(i);
%     yhat=y+sigma*randn(size(y));
%     
%     %%
%     lambdalist=-6:0.5:6;
%     %lambdalist=2:0.5:3;%-6:0.2:6;
%     lambdalist=10.^lambdalist;
%     %lambdalist=0.02;
%     sz=length(lambdalist);
%     mse_list_inner=100*ones(sz,1);
%     reg_weight=zeros(sz,1);
%     for j=1:sz
%         m=size(D,1);
%         scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
%         lambda=scalefactor*lambdalist(j); rho=0.5*lambda;
%         
%         [ x ,history1] = gtf( yhat,edges,k,lambda,rho);
%         %[ x,history ] = gtf_proj_newton1( yhat,edges,k,lambda);
% 
%         mse_list_inner(j)=norm(x-y)^2/n;
%         reg_weight(j)=norm(O*x);
%         
% 
%  
%         
%     end
%             loglog(lambdalist,mse_list_inner);
%         hold on;
%     mse_list2(i)=min(mse_list_inner);
%            pause(1);
%     %save('temp.mat');
% end



%%

%save('Exp_Facebook_results.mat','y','sigma_list','mse_list','mse_list2','mse_list3')
%%
h=figure;
plot(sigma_list,mse_list,'ro-',...
    sigma_list,mse_list1,'r*-',...
    sigma_list,mse_list2,'r^-',...
    sigma_list,mse_listL,'kx--',...
    sigma_list,mse_listW,'bs-.','linewidth',2,'markersize',12);
xlabel('Noise level \sigma','fontsize',14)
ylabel('MSE','fontsize',14)
lg=legend('Trend filtering k=0','Trend filtering k=1','Trend filtering k=2','Laplacian smoothing','Wavelet smoothing');
set(lg,'fontsize',14,'location','best')
grid on;

%saveTightFigure(h,'facebook_poisson_sparse.pdf')