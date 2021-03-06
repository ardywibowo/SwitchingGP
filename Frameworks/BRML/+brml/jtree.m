function [jtpot jtsep infostruct]=jtree(pot)
%JTREE Setup a Junction Tree based on a set of potentials
% [jtpot jtsep infostruct]=jtree(pot)
% For the potential pot return a set of initialised Junction Tree potentials.
% Input : pot  -- a cell vector of potentials.
% Output:
% jtpot : the clique potentials.
% jtsep : the separator potentials.
% infostruct.cliquetree : the Junction Tree as an adjacency matrix
% infostruct.sepind : sepind(i,j)=separator index between cliques i and j
% infostruct.EliminationSchedule : Each row contains in the first position
% the clique to be eliminated, and the second position the clique that accepts this eliminated clique.
% Note that the variables in the potentials must be numbered 1,2,3,... with
% no missing numbers. Use squeezepots to achieve this if necessary
% See demoJTree.m
import brml.*
pot=str2cell(pot);
variables=potvariables(pot);
if variables(end)~=length(variables); warning('potential variables used are not numbered 1:end. Use squeezepots.'); end

% triangulation:
[dum1 cl] = triangulate(markov(pot)); % just need the cliques
C=length(cl);
% get the edge weights (number of variables in each separator):
w=sparse(C,C);
for c1=1:C-1
    for c2=c1+1:C
        w(c1,c2)=length(intersect(cl{c1}.variables,cl{c2}.variables));
        w(c2,c1)=w(c1,c2);
    end
end
Atree=spantree(w);

% Setup the potentials:
Aup=triu(Atree); sepind=Aup;
sepind(find(sepind))=1:sum(sepind(:)); % number of separators must be C-1 for a tree
sepind=sepind+sepind'; % separator index : sepind(jtclique_i,jtclique_j)=k, k=1:C-1
infostruct.separator=[];
for s=1:max(sepind(:))
    [c1,c2]=find(triu(sepind)==s);
    sepvariables = intersect(cl{c1}.variables,cl{c2}.variables);
    infostruct.separator{s}=sepvariables;
end
for c=1:C
    infostruct.cliquevariables{c}=cl{c}.variables;
end
infostruct.cliquetree=Atree; % connectivity structure of JT cliques
infostruct.sepind=sepind; % the indices of separators sepind(i,j)=separator index
[tree , elimseq, sched]=istree(Atree); % choose an elimination order
infostruct.EliminationSchedule=sched;
[jtpot jtsep]=jtassignpot(pot,infostruct); % initialise potentials