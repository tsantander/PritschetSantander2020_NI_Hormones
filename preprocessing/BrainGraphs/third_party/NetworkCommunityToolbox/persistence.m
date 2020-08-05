function [ pers ] = persistence( S )
% this function computes the persistence for a given multilayer partition S
% defined in  "Community detection in temporal multilayer networks, and ...
% its application to correlation networks"
% (http://arxiv.org/abs/1501.00040)
% input:
%    S: pxn matrix, n is the numebr of partitions, p is the number of nodes
% output:
%    pers: the persistence value  
pers = sum(sum(S(:,2:end) == S(:,1:end-1)));
end

