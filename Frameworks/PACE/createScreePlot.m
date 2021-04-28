%function createScreePlot(yy)
%This function creates the screePlot based on the
%results from FPCA() or FPCder()
%Input yy : returned object from FPCA().
%example:
%yy = FPCA(y,t,p);
%createScreePlot(yy)
%or
%p = setDerOptions('nder',0:2);
%yy = FPCder(y,t,p);
%createScreePlot(yy)     
function createScreePlot(yy)

  ops = getVal(yy,'ops');
  yname = getVal(ops,'yname');
  createSP(getVal(yy,'FVE'), getVal(yy,'no_opt'), yname)

end
