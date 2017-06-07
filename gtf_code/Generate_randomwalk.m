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

s=20;
idx=randperm(n);
idx=[idx(1:s), headnodes];
s=length(idx);



L=D'*D; A=-L;
A = spdiags(zeros(n,1),0,A); 

% signal from head nodes and 20 other nodes.
decay_rate=0.5; %0.5 decay every 3 stesp.

y=zeros(n,1);

for j=1:10
    for i=1:length(idx)
        Ind = idx(i);
        F = @(t,X) 0;
        G = @(t,X) 1;
        S = sde(F,G,'startState',Ind);
        X = S.simByEuler(10*j,'ntrials',ceil(10000*0.5^j),'Z',@(t,X) RandomGraphMove(X,A)-X);
        X = squeeze(X);
        
        Visited = unique(X(end,:));
        for ii = 1:numel(Visited)
            y(ii)=y(ii)+sum(X(end,:) == Visited(ii));
        end
    end
end