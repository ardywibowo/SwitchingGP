function demoSVDmissing
%DEMOSVDMISSING demo of performing Singular Value Decomposition with missing data
import brml.*
B=randn(10,1); W=randn(1,500); X=B*W; % create data generated from 5 basis vectors
G=real(rand(size(X))<0.2); 
Xmissing=X.*replace(G,0,nan); % create missing data
opts.displayerror=1;
opts.tolerance=10e-7;
[U,S,V]=svdm(Xmissing,5,opts); Xr=U*S*V';
disp(['Mean Square Reconstruction error of all (including non-missing) data = '...
    num2str(mynanmean(mynanmean((X-Xr).^2)))])
maxval=max(max(Xmissing)); % don't plot extreme values
subplot(3,1,1); imagesc(Xmissing); title('X with missing data');
subplot(3,1,2); imagesc(cap(Xr,maxval)); title('SVDm reconstruction');
subplot(3,1,3); imagesc(cap(X,maxval)); title('original X ');