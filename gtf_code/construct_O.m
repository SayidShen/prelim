function [ O ] = construct_O( D,k )
%UNTITLED4 Summary of this function goes here
%   k is 0 based.
O=D;
for i=1:k
    if mod(i,2)
        O= D'*O;
    else
        O=D*O;
    end
end

end

