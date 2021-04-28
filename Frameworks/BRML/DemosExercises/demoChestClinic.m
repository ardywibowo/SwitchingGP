function demoChestClinic
%DEMOCHESTCLINIC demo of inference for the chest clinic belief net
import brml.*

[asia smoker tub lcancer bronch xray dys tlc]=assign(1:8); % Variable order is arbitary
yes=1; no=2; % define states, starting from 1
varnames={'asia' 'smoker' 'tub' 'lcancer' 'bronch' 'xray' 'dys' 'tlc'};
for v=1:length(varnames);
    variable(v).name=varnames(v); variable(v).domain={'yes','no'};  % all variables have same domain
end

% For a BN associate p(i|pa(i)) with potential i: 
pot(asia)=array(asia); pot(asia).table(yes)=0.01; pot(asia).table(no)=1-pot(asia).table(yes);

pot(smoker)=array(smoker); pot(smoker).table(yes)=0.5; pot(smoker).table(no)=1-pot(smoker).table(yes);

pot(tub)=array([tub asia]); % define array below using this variable order
tmptable(yes, yes)=0.05; tmptable(yes, no)=0.01;  
tmptable(no,:)=1-tmptable(yes,:); % due to normalisation
pot(tub).table=tmptable; 

pot(lcancer)=array([lcancer smoker]); % define array below using this variable order
tmptable(yes, yes)=0.1; tmptable(yes, no)=0.01;  
tmptable(no,:)=1-tmptable(yes,:); % due to normalisation
pot(lcancer).table=tmptable; 

pot(bronch)=array([bronch smoker]); % define array below using this variable order
tmptable(yes, yes)=0.6; tmptable(yes, no)=0.3;  
tmptable(no,:)=1-tmptable(yes,:); % due to normalisation
pot(bronch).table=tmptable; 

pot(xray)=array([xray tlc]); % define array below using this variable order
tmptable(yes, yes)=0.98; tmptable(yes, no)=0.05;  
tmptable(no,:)=1-tmptable(yes,:); % due to normalisation
pot(xray).table=tmptable;  

pot(dys)=array([dys tlc bronch]); % define array below using this variable order
tmptable(yes, yes, yes)=0.9; tmptable(yes, yes, no)=0.7;
tmptable(yes, no, yes)=0.8; tmptable(yes, no, no)=0.1;
tmptable(no,:,:)=1-tmptable(yes,:,:); % due to normalisation
pot(dys).table=tmptable; 

pot(tlc)=array([tlc tub lcancer]); % define array below using this variable order
tmptable=ones([2 2 2]);
tmptable(yes, no, no)=0;  tmptable(no,:,:)=1-tmptable(yes,:,:); % due to normalisation
pot(tlc).table=tmptable; 

drawNet(dag(pot),variable);

disp('Sum the full joint probability table (inefficient):');
jointpot = multpots(pot(1:8)); % joint distribution
% p(dys=yes)
%margpot = sumpot(jointpot,setdiff(1:8,dys));
margpot = sumpot(jointpot,dys,0); % simpler syntax
disp(['p(dys=yes) ' num2str(margpot.table(yes)./sum(margpot.table))]);

% p(dys=yes|smoker=yes)
margpot = sumpot(setpot(jointpot,smoker,yes),dys,0);
disp(['p(dys=yes|smoker=yes) ',num2str(margpot.table(yes)./sum(margpot.table))]);

% p(dys=yes|smoker=no) % alternative way using condpot
cpot = condpot(setpot(jointpot,smoker,no),dys);
disp(['p(dys=yes|smoker=no) ' num2str(cpot.table(yes))]);
disp('cond pot:'); disptable(condpot(jointpot,dys,smoker),variable);