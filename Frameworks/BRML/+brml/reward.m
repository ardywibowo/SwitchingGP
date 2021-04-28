function r=reward(s,par)
r=0;
if any(s==par.rewardstate)
    r=1;
elseif any(s==par.penaltystate)
    r=-1;
end