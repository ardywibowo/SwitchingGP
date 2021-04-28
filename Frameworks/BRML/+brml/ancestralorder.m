function [ord noparents_vars]=ancestralorder(A)
%ANCESTRALORDER Return the ancestral order for the DAG A (oldest first)
% ord=ancestralorder(A)
% If A is not a DAG, return the empty set for the order
N=size(A,1); AA=A; ord=[]; noparents_vars=[];
done=zeros(1,N);
while length(ord)<N
    oldlength=length(ord);
	for i=1:N
		nochildren = isempty(find(AA(i,:)));
		noparents = isempty(find(AA(:,i)));
		noparentsA = isempty(find(A(:,i)));
		if noparentsA; noparents_vars = unique([noparents_vars i]); end
		if ~(noparents & nochildren)
			if nochildren % i has no children
				AA(:,i)=0;
				ord=[ord i];
				done(i)=1;
			end
		end
	end
	if all(done); break; end
	rest=setdiff(1:N,find(done));
	if isempty(find(AA(:,rest)))
		break;
    end
    if length(ord)==oldlength
        ord=[]; % not a DAG
        disp('not a DAG')
        return
    end
  
end
ord=fliplr([ord rest]);