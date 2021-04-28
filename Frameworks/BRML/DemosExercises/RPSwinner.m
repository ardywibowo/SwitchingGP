function w=RPSwinner(o,c)
% see if o beats c
import brml.*
[R P S]=assign(1:3);
if  (o==R)&(c==S) | (o==P)&(c==R) | (o==S)&(c==P)
    w=1;
elseif  (o==R)&(c==R) | (o==P)&(c==P) | (o==S)&(c==S)
    w=0; % draw
else
    w=-1; % loss
end