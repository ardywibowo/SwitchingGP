function newpot= setpot(pot,evvariables,evidstates)
%SETPOT sets potential variables to specified states
% newpot = setpot(pot,variables,evidstates)
%
% set variables in potential to evidential states in evidstates
% Note that the new potential does not contain the evidential variables
import brml.*
[vars nstates] = potvariables(pot);
if isempty(intersect(vars,evvariables))
    newpot=pot;
else
    nonevidentialvariables=setminus(vars,evvariables);
    thispotevidentialvariables=setminus(vars,nonevidentialvariables);   
    [a bind]=ismember(thispotevidentialvariables,vars);
    tmppot=brml.array;
    tmppot.variables=thispotevidentialvariables;
    tmppot.table=myzeros(nstates(bind));
    tmppot2=setstate(tmppot,evvariables,evidstates,1);
    newpot = sum(prod([tmppot2 pot]),thispotevidentialvariables);
end