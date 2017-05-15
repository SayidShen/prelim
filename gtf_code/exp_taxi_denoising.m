% Taxi denoising experiments

datapath='/Users/yxiangw/Documents/bitbucket/TaxiGTF/data'

%% Get the graph

n=3874;
m=7070;
t=24;

xy=load([datapath,'/junction_xy.txt']); 
tmp=xy(:,1);
xy(:,1)=xy(:,2);
xy(:,2)=tmp;


edges=load([datapath,'/road_graph.txt']);

w=edges(:,3);%edge weight
edges=edges(:,1:2)+1;%change to 1 based indexing in Matlab

% spatial temporal graphs
nn=n*t;
mm=m*t+n*(t-1);
EDGES=kron(ones(t,1),edges)+kron((0:t-1)'*n,ones(m,2));
EDGES2=[kron(ones(t-1,1),(1:n)')+kron((0:t-2)'*n,ones(n,1)),...
        kron(ones(t-1,1),(1:n)')+kron((1:t-1)'*n,ones(n,1))];
EDGES=[EDGES;EDGES2];


%% get all the days
%represent date in dta number format


DateStart=datenum(2011,1,1);
DateEnd=datenum(2012,12,31);
windowsz=8;

N=DateEnd-DateStart+1;
for iter=408:N %Day iteration
        
    % get current day
    dnum=datevec(iter+DateStart-1);
    Year=dnum(1);
    Mon=dnum(2);
    Day=dnum(3);
    yhat=log(load_a_day(Year,Mon,Day,n,t)+1);

    %construct seasonal averages for the day of the week
    [year,month,day]=find_shiftingwindow(DateStart,DateEnd,iter,windowsz);
    sz=length(year);
    y=zeros(size(yhat));
    for i=1:sz
        %note the data preprocessing
        y=y+log(load_a_day(year(i),month(i),day(i),n,t)+1);
    end
    y=y/sz;
    


%% Hold out data


% Cross validation


% Run GTF