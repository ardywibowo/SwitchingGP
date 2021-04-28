function [m,tie]=majority(x,varargin)
%MAJORITY Return majority values in each column on a matrix
% [m,tie]=majority(x,<checktie>)
% return the majority value in each column of x. If there is a tie the lowest value is returned
% if checktie is true, then tie identifies which columns have ties.
checktie=false;
if nargin==2
    checktie=varargin{1};
end
for i=1:size(x,2);
    u=unique(x(:,i)); s=zeros(1,length(u));
    for j=1:length(u)
        s(j)=length(find(x(:,i)==u(j)));
    end
    [val ind]=max(s);
    m(i)=u(ind);
    if checktie
        tie(i)=sum(s==val)>1;
    end
end
if ~checktie
    tie=[];
end