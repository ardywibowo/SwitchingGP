function c=compat(v,F,h,Gx,Gy)
%COMPAT Compatibility of object F being in position h for image v on grid Gx,Gy
% c=compat(v,F,h,Gx,Gy);
import brml.*
%c= sum(sum((v.*placeobject(F,h,Gx,Gy))))./sum(sum(F));
c= sum(sum((v.*placeobject(F,h,Gx,Gy))))./sum(sum(placeobject(F,h,Gx,Gy)));