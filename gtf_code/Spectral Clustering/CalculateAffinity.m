function [affinity] = CalculateAffinity(data,varargin)

% set the parameters
if nargin>1
    sigma=varargin{1};
else
    sigma = 1;
end

for i=1:size(data,1)    
    for j=1:size(data,1)
        dist = sqrt((data(i,1) - data(j,1))^2 + (data(i,2) - data(j,2))^2); 
        affinity(i,j) = 1/dist;
        %exp(-dist^2/(2*sigma^2));
    end
end


