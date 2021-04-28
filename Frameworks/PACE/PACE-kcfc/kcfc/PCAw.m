%===========
%Description:
%===========
%
%            This is the main program to perform Functional Principal Component 
%            Analysis (FPCA) via PACE. The principal component scores can be 
%            estimated through conditional expectation or via classical integration. 
%            For the latter, one can choose a shrinkage method for estimated scores.
%            See Reference [2].  For shrinkage, see Reference [1].
%
%======
%Usage:
%======
%
% function [no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,bw_mu,xcov,bw_xcov,
%           xcovfit, AIC,BIC,FVE,y_pred, y_predOrig, out1,out21,...
%           y,t, regular, rho_opt, mucopy, phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr]...
%        = PCAw(y,t,bwmu,bwmu_gcv, bwxcov,bwxcov_gcv,ntest1,ngrid1,selection_k,
%          FVE_threshold, maxk,control,regular,error, ngrid,method,shrink,...
%          newdata,kernel, numBins, yname, screePlot, designPlot, corrPlot, rho, verbose,method_mu, userrange)
%
%================
%Input Arguments:
%================
%
%Input y:          1*n cell array, y{i} is the vector of measurements for the
%                  ith subject, i=1,...,n.                     
%     
%Input t:          1*n cell array, t{i} is the vector of time points for the
%                  ith subject for which corresponding measurements y{i} are
%                  available, i=1,...,n.
%
%Input bwmu:       scalar bwmu>=0, bandwidth for mean curve mu,
%                  0:  use cross-validation or generalized cross-validation
%                      to choose bandwidth automatically.      [Default]
%                  bwmu(>0): user-specified bandwidth.
%
%Input bwmu_gcv:   For choice bwmu = 0, two types of cross-validation can
%                  be performed: One-curve-leave-out cross-validation (CV)
%                  or  generalized cross-validation (GCV)
%                  0: CV   (may be time-consuming)
%                  1: GCV  (faster)                            [Default]
%                  2: Geometric mean between the minimum bandwidth and the
%                     GCV bandwidth.
%
%Input bwxcov:     1*2 vector, bandwidths for covariance surface used for
%                  smoothing of cov(X(t),X(s))
%                  bwxcov(i): ith coordinate of bandwidth vector, i=1,2.
%                  bwxcov(1)==0 & bwxcov(2)==0: use cross-validation (CV)
%                  or generalized cross-validation (GCV) for automatic
%                  selection.                                  [Default]
%                  bwxcov(1)>0 & bwxcov(2)>0: user-specified bandwidths.
%
%Input bwxcov_gcv: If setting bwxcov = [0 0], automatic bandwidth selection
%                  by CV or GCV choices
%                  0: CV method (may be time-consuming, can be accelerated by 
%                     choosing small values for ntest1 and ngrid1).
%                  1: GCV method (faster)                      [Default]
%                  2: Geometric mean between the minimum bandwidth and GCV
%                     bandwidth
%
%Input ntest1:     integer(<=n), number of curves used for CV when choosing
%                  bandwidths for smoothing the covariance surface. The subjects 
%                  in the test set are randomly selected from n subjects. Small 
%                  ntest1 will accelerate CV at less accuracy. [Default is 30.]
%
%Input ngrid1:     integer, number of support points for the covariance surface 
%                  in the CV procedure (selecting bandwidths of covariance).
%                  Note that input ntest1 and ngrid1 provide options to save
%                  computing time when using CV or GCV.        [Default is 30.]
%
%Input selection_k: the method of choosing the number of principal components K.
%                   'AIC1': use AIC criterion with pseudo-likelihood of
%                           measurements (marginal likelihood).
%                   'AIC2': use AIC criterion with likelihood of measurements
%                           conditional on estimated random coeffecients.
%                   'BIC1': use BIC criterion with pseudo-likelihood of
%                           measurements (marginal likelihood).
%                                                               [Default]
%                   'BIC2': use BIC criterion with likelihood of measurements
%                           conditional on estimated random coeffecients.
%                   'FVE' (fraction of variance explained) : use scree plot
%                           approach to select number of principal
%                           components), see "FVE_threshold" below.
%                    positive integer K: user-specified number of principal 
%                                        components
%
%                   Note: BIC1 and FVE produce the most parsimonious models.
%
%Input FVE_threshold:  a positive number that is between 0 and 1 [Default is 0.85.]
%                      It is used with the option selection_k = 'FVE' to select
%                      the number of principal components that explain at least
%                      "FVE_threshold" of total variation (the fraction
%                      of variance explained).
%
%Input maxk:      integer, the maximum number of principal components to consider
%                 if using automatic methods to choose K, i.e., 'AIC1', 'AIC2',
%                 'BIC1' or 'BIC2' defined by selection_k.      [Default is 20.]
%                 Note: when selection_k = 'FVE' or 'AIC_R', maxk is ignored.
%
%Input control:   'auto', Select K by minimizing AIC or BIC, or find the
%                         first K such that the FVE_threshold is exceeded. [Default]
%                 'look', a scree plot (FVE% Vs No. of PC) will be generated based
%                         on K <= 15. User will be prompted to enter user-specified
%                         K after viewing scree plot. This can be combined
%                         with any setting of selection_k.
%
%Input regular:   0, sparse (or irregular) functional data.      
%                 1, regular data with missing values
%                 2, completely balanced (regular) data.
%                 [], automatically chosen based on the design in t.   [Default]
%
%Input error:     0, no additional measurement error assumed.
%                 1, additional measurement error is assumed.    [Default]
%
%Input ngrid:     integer, number of support points in each direction of 
%                 covariance surface when performing principal component 
%                 analysis ( ngrid > K).                    [Default is 51.]
%Input method:    used for computing random effects \xi_{ik}
%                 'CE': conditional expectation method           [Default]
%                 'IN': classical integration method
%                 Note: 'CE' can be applied for sparse data or regular data, but
%                       'IN' only in the case of regular data.
%
%Input shrink:    indicator of whether applying shrinkage to estimates of random
%                 coefficients (for regular data only)
%                 0:  no shrinkage when method = 'CE' or error = 0 [Default]
%                 1:  shrinkage when method = 'IN' and error = 1, otherwise, this
%                     will be re-set to 0.
%
%Input newdata:   a row vector of user-defined output time grids for
%                 all curves. This corresponds to "out1" in the output argument
%                 If newdata = [], then "out1" corresponds to the set of distinct
%                 time points from the pooled data.
%                 "newdata" is supposed to be a vector in ascending order on
%                  the domain of the functions.                    [Default is []]
%
%Input kernel:    a character string to define the kernel to be used in the
%                 1-D or 2-D smoothing
%                 kernel = 'epan'  ==> Epanechnikov kernel [Default for dense designs with n_i >= 20]
%                          'rect'  ==> Rectangular kernel
%                          'gauss'  ==> Gaussian kernel    [Default for sparse designs, regular designs with
%                                                           missings, dense designs for n_i < 20]
%                 Note: The Gaussian kernel is overall best for sparse designs but is slower than the other kernels 
%                       and if computational speed is of concern then one may wish to use the Epanechnikov kernel 
%                       also for the case of sparse designs.
%
%Input numBins:   0: no binning
%                 a positive interger (>= 10): prebin the data with user-defined
%                 number of bins. When numBins < 10, no binning will be performed.
%		  []:  prebin the data with the following rule    [Default]
%
%                 i) When the input data is regular = 1 or 2
%                    m = max of n_i, where n_i is the number of repeated measurements
%                    for ith subject.
%                 ii) regular = 0
%                    m = median of n_i
%
%                 When m <= 20 subjects, no binning.
%                 When n <= 5000 subjects and m <= 400, no binning.
%                 When n <= 5000 subjects and m > 400, numBins = 400.
%                 When n > 5000 subjects, compute
%
%                 m* = max(20, (5000-n)*19/2250+400)
%
%                 if m > m*, numBins = m*
%                 if m <= m*, no binning
%
%                 This option accelerates analysis, especially for data with
%                 many time measurements.
%
%Input yname:     a character string which denotes the name of the current
%                 function to be estimated.               [Default is []]
%                 It is used in the title part of the scree plot output.
%                 If it is set to [], then yname is set to be the same
%                 name as the first input argument from PCA() or FPCA().
%                 For example, if this analysis is for function y,
%                 then yname = 'y'.
%
%Input screePlot: indicator of whether to create the scree plot
%                 1  a scree plot will be created         
%                 0  no scree plot will be created        [Default]
%Input designPlot: indicator of whether to create the design plot
%                 1  a design plot will be created
%                 0  no design plot will be created       [Default]
%                 Interpretation of design plot: All subareas of the
%                 domain x domain support square of the covariance surface need
%                 to be covered more or less homogeneously with pairs of design points.
%Input corrPlot:  indicator of whether to create the correlation surface plot
%                 1  a correlation surface plot will be created
%                 0  no correlation surface plot will be created       [Default]
%
%Input rho:       truncation threshold for the iterative residual that is used
%                 in the estimation of FPC scores. (see FPCscore.pdf under Help/ for more details)
%                 -1:  compute unadjusted FPC scores (as in previous PACE versions)
%                 >0:  user-defined choice of rho
%                 0:   do not set truncation threshold, but use iterative residuals for sigmanew (see in output description below)
%                 'cv-random':  a character string which specifies to use randomized leave-one-measurement-out CV approach to find
%                        the optimal value of rho.        
%                 Note that this choice contains a random element and therefore the analysis is not exactly
%                  replicable when running  the program twice.
%                 'cv': use non-randomized leave-one-measurement-out CV approach to find the optimal value of rho. [default]
%
%Input verbose:   a character string
%                 'on': display diagnostic messages       [Default]
%                 'off': suppress diagnostic messages
%
%Input xmu:       1*N vector, user-defined mean functions valued at distinct
%                 input time points, in ascending order from all subjects, corresponding to out1. 
%
%
%Input method_mu: method to estimate mu
%                 'RARE': using response-adaptive method to update
%                 estimated       
%                 'PACE': using the original method   [default]
%Input out_percent: a number between 0 and 1 indicates that we leave out out_percent data in the boundary, i.e.,
%                   we shrink the domain. If out_percent > 0.25, reset to 0.25. The dafault is no shrinkage.
%                  When performing local linear smoothing for mu and xcov, 
%                 all the data are used, but the output grid out1 and out21 will be restricted within
%                 the reduced range. This option is used to alleviate the boundary
%                 effect.
%
%Input weight:    option of case weight in the local weighted least square
%                 1:  taking case weight as inverse of number of observation for each subject, i.e., 1/mi
%                 0:  taking case weight as one for all subjects
%========
% Details:
%========
%       1)       bwmu = 0 and bwmu_gcv = 1 ==> use GCV method to choose
%                the bandwidth for mean function
%
%                For Gaussian kernel, the optimal bandwidth from GCV is multiplied
%                by 1.1.
%
%                bwmu = 0 and bwmu_gcv = 0 ==> use CV method to choose the
%                bandwidth for mean function
%
%                bwmu = 0 and bwmu_gcv = 2 ==> use the geometric mean
%                between the minimum bandwidth and GCV bandwidth for mean
%                function
%
%                bwmu > 0 then bwmu_gcv will be ignored, subjective
%                bandwidth choice
%
%       2)       bwxcov = [0 0] and bwxcov_gcv = 1 ==> use GCV method to
%                choose the bandwidth for cov function
%
%                For Gaussian kernel, the optimal bandwidth from GCV is multiplied
%                by 1.1.
%                
%                bwxcov = [0 0] and bwxcov_gcv = 0 ==> use CV method to
%                choose the bandwidth for cov function
%                bwxcov = [0 0] and bwxcov_gcv = 2 ==> use the geometric mean
%                between the minimum bandwidth and the GCV bandwidth for cov
%                function 
%
%                bwxcov(1) > 0 and bwxcov(2) > 0 then bwxcov_gcv will be
%                ignored, subjective bandwidth choice
%
%                If eigenvalues estimation is of primary interest, you may want
%                to undersmooth the covariance surface for better estimates 
%                of eigenvalues.
%
%       3)       If error=1, AIC1, BIC1, AIC2, BIC2, FVE all can be used
%                for choosing K.
%
%                If error=0, only AIC1, BIC1 and FVE can be used by
%                definition of these criteria.
%
%       4)       If some arguments are not needed or you want to use
%                default values, just use [] for these arguments.
%
%                Alternatively, use setOptions() to define nessary input
%                arguments and call PCA through FPCA.
%
%       5)       In general, any output in the form of 1-dimensional
%                vectors will be a row vector and any output related to
%                1-dimensional cell arrays will be a row cell array.
%
%=================
%Output Arguments:
%=================
%Output no_opt:   integer, automatically or subjectively selected value of K, the number of selected components.
%
%Output sigma:    scalar, estimate of measurement error variance if
%                 error=1, while it is [] if error=0.
%
%Output lambda:   1*K vector, estimated eigenvalues (variances of functional principal components scores).
%
%Output phi:      N*K matrix, estimated principal component functions
%                 valued at distinct input time points with ascending
%                 order of all subjects, corresponding to out1
%
%Output eigen:    ngrid*K matrix, estimated principal component functions,
%                 valued at out21, ngrid of the pooled distinct time points
%                 with ascending order of all subjects,
%                 phi is an interpolated version of eigen at out1
%
%Output xi_est:   n*K matrix, predictions for random coeffecients (PC
%                 scores) for n subjects.
%
%Output xi_var:   K*K matrix, Var(PC score)-Var(estimated PC score). The
%                 omega matrix in equation (7) of the paper, which is used
%                 to construct the point-wise C.I. for X_i(t)
%
%Output mu:       1*N vector, estimated mean functions valued at distinct
%                 input time points (newdata = []), in ascending order from
%                 all subjects, corresponding to out1; when newdata is defined,
%                 corresponds to the time points from newdata, same as
%                 out1.
%
%Output bw_mu:    scalar(>0), automatically or subjectively selected
%                 bandwidth for smoothing mean curve.
%
%Output xcov:     ngrid*ngrid matrix, smoothed covariance surface (diagnal
%                 removed), corresponding to out21
%
%Output bw_xcov:  1*2 vector(>0), automatically or subjectively selected
%                 bandwidths for smoothing covariance surface.
%
%Output xcovfit:  ngrid * ngrid matrix, fitted covariance surface, based
%                 on truncated estimate of eigenvalues ("lambda") and
%                 principal component functions ("eigen"), corresponding
%                 to out21
%
%Output xcorr:    ngrid * ngrid matrix, fitted correlation surface, based
%                 on truncated estimate of eigenvalues ("lambda") and
%                 principal component functions ("eigen"), corresponding
%                 to out21
%
%Output AIC:      1*K vector, AIC values obtained when choosing K from
%                 K=1 to K=maxk, where AIC(K) is the minimum. If AIC
%                 method is not applied, it is []
%
%Output BIC:      1*K vector, BIC values obtained when choosing K from
%                 K=1 to K=maxk, where BIC(K) is the minimum. If BIC
%                 method is not applied, it is []
%
%Output FVE:      1*ngrid vector of fraction of variance explained
%
%Output y_pred:   cell array, y_pred{i} is the vector of predictions for
%                 the ith subject evaluated at time points from the output
%                  grid vector "out1".
%
%Output y_predOrig: cell array, y_predOrig{i} is the vector of predictions
%                   for the ith subject at the same time points as the input.
%
%Output out1:     1*N vector, distinct input time points with ascending
%                 order from all subjects (if out_percent is specified, it is restricted to be
%                  within the smaller range) if newdata = []; otherwise, it
%                 is the same as newdata.
%
%Output out21:    1*ngrid vector, a grid of time points for which the
%                 smoothed covariance surface assumes values, i.e.,
%                 ngrids from out1.
%
%Output y:        if no binning is performed, same as the input y
%                 if binning is performed, 1 * n cell array, y{i} is a vector
%                 of measurements after binning for subject i, i = 1,...,n
%                 If out_percent is specified, it is restricted to be
%                 within the smaller range.
%
%Output t:        if no binning is performed, same as the input t
%                 if binning, 1 * n cell array, t{i} is a vector of
%                 time points after binning for subject i, i = 1,...,n
%                 Each value of t{i} corresponds to the midpoints of
%                 each bin.
%                 If out_percent is specified, it is restricted to be
%                 within the smaller range.
%
%Output regular:  if no binning is performed or regular = 2, this is the
%                 same as the input
%                 if binning is performed and regular = 0, it will be
%                 reset to regular = 1. In other words, after binning,
%                 the sparse and irregular case is analyzed as regular
%                 with missings by sampled data.
%
%Output rho_opt:  if rho is set as 'cv' or 'cv-random', then rho_opt is the optimal rho obtained
%                  from the CV method, otherwise, it is the same as the input rho
%                 When rho is 'cv', rho_opt is non-random for a give data set
%                 When rho is 'cv-random', rho_opt can be different due to the randomly leave-out measurements
%                 for the same data set.
%
%Output sigmanew: if rho is set as >0, 0 or 'cv', then sigmanew is the iterative
%                 residual sum of squares (see FPCscore.pdf for more details). It
%                 can be used as an estimate of the variance of the measurement
%                 errors.
%                 if rho is set as -1, then sigmanew is set to the same as output
%                 sigma.
% See also FPCA, setOptions, showOptionNames, example
function [no_opt,sigma,lambda,phi,eigen,xi_est,xi_var,mu,muDense,bw_mu,xcov,bw_xcov, xcovfit, AIC,BIC,FVE,y_pred, y_predOrig, y_predDense, out1,out21,...
	  y,t, regular, rho_opt, sigmanew, mucopy, phicopy, eigencopy, out1copy, out21copy, xcovcopy, xcovfitcopy, xcorr, ops]... 
         =PCAw(y,t,bwmu,bwmu_gcv, bwxcov,bwxcov_gcv,ntest1,ngrid1,selection_k, FVE_threshold, maxk,...
              control,regular,error,ngrid,method,shrink,newdata,kernel, numBins,yname, screePlot, designPlot, corrPlot, rho, verbose, xmu, xcov, method_mu, out_percent,weight)
%clear
%clc

%load PBC
%addpath (genpath('PACE2.15'));
%p = setopts('selection_k','FVE','FVE_threshold',0.9,'kernel','epan',...
%    'rho',-1,'ngrid',m_pt,'regular',2,'method','IN');
%y = albumin; t = Tin;
%bwmu = 0;bwmu_gcv = 1; bwxcov = [0 0]; bwxcov_gcv = 1; ntest1 = 30; 
%ngrid1 = 30; selection_k = 'FVE'; FVE_threshold=0.9; maxk = 20;
%control = 'auto'; regular = 0; error = 1; ngrid = 51; method = 'CE'; 
%shrink = 0; newdata = []; kernel = 'epan'; numBins = []; yname = []; 
%screePlot = 0; designPlot = 0; corrPlot = 0; rho = -1; verbose = 'on'; 
%xmu = []; xcov = []; weight = 1;method_mu=[];
%out_percent=[];

     no_opt =[];sigma=[];lambda=[];phi=[];eigen=[];xi_est=[]; 
     xi_var=[];mu=[];muDense=[];bw_mu=[];bw_xcov=[];xcovfit=[];AIC=[];BIC=[];FVE=[];
     y_pred =[];y_predOrig = [];y_predDense=[];out1 = [];out21 =[];
     rho_opt=[]; sigmanew = []; mucopy=[];phicopy=[];eigencopy=[];out1copy=[];out21copy=[];xcovcopy=[];xcovfitcopy=[];xcorr=[]; ops=[];

invalid = 0;
userrange = [];
[invalid_data]  = CheckData(y, t); %NaN values and number of subjects.
if (invalid_data == 1)
    return
end

if isempty(yname)
    yname = inputname(1);
end

%Set default values for the major input arguments
%when the following input arguments are set to 
%be "[]" by the user

if isempty(bwmu)
    bwmu = 0;    %bandwidth choice for mean function is using CV or GCV
end
if isempty(bwmu_gcv)
    bwmu_gcv = 1; %bandwidth choice for mean function is GCV
end
if isempty(bwxcov)
    bwxcov = [0 0]; %bandwidth choices for covariance function is CV or GCV
end
if isempty(bwxcov_gcv)
    bwxcov_gcv = 1; %bandwidth choices for covariance function is GCV
end
if isempty(ntest1)
    ntest1 = 30;
end
if isempty(ngrid1)
   ngrid1 = 30;
end

if isempty(selection_k) 
     selection_k = 'BIC1';
end

if strcmp(selection_k, 'FVE') && isempty(FVE_threshold)
     FVE_threshold = 0.85;
end

if isempty(error)
    error = 1;       %error assumption with measurement error
end


if isempty(maxk)
    maxk = 20;      %20 PC candidates
end

if isempty(shrink)
     shrink = 0;
end

if isempty(method)
    shrink = 0;
    method = 'CE';   %method to estimate the PC score is through conditional expectation
end


if shrink == 1 && (error ~= 1 || strcmp(method,'IN') ~= 1)
  fprintf(1,'Warning: shrinkage method only had effects when method = "IN" and error = 1 !Reset to shrink = 0 now!\n');
  shrink = 0; 
end

if isempty(control)
    control = 'auto';%automatically choose the number of PC not through interactive
                     %input by the user after visualization
end

if error == 1 
    cut = 1;         %cut off the bounary when estimating sigma
else 
    cut = 0;
end


if isempty(ngrid)
    ngrid = 51;      %number of output time grids is 30
end

if maxk > (ngrid-2)
    fprintf(1,'Warning: maxk can only be less than or equal to ngrid-2! Reset it to be ngrid now!\n');
    maxk = ngrid-2;
end

if isnumeric(selection_k)
  if selection_k > (ngrid-2)
     fprintf(1,'Warning: selection_k can only be less than or equal to ngrid-2! Reset it to be ngrid now!\n');
     selection_k = ngrid-2;
  elseif selection_k <= 0
     fprintf(1,'Warning: selection_k must be a positive integer! Reset it to BIC1 now\n');
     selection_k = 'BIC1';
     FVE_threshold = 0.85;
  end
end

if isempty(screePlot)
   screePlot = 0;
end

if isempty(designPlot)
    designPlot = 0;
end

if isempty(corrPlot)
    corrPlot = 0;
end

if isempty(rho)
    rho = 'cv';
end

if isempty(verbose)
  verbose = 'on';
end

if error == 0 && (strcmp(selection_k,'AIC2') || strcmp(selection_k,'BIC2'))
     fprintf(1,'Warning: When assume no measurement error, cannot use "AIC2" or "BIC2". Reset to "BIC1" now!\n');
     selection_k = 'BIC1';
end

 ncohort=length(t);     % obtain the number of curves or subjects

 ni = zeros(1,ncohort);
 for i = 1:ncohort
   ni(i)= length(t{i});
 end
 if all(ni == 1)
     fprintf(1,'Error: FPCA is aborted because the data do not contain repeated measurements!');
     return;
 end

 if isempty(regular)||~any(regular==0:2)
        regular = isregular(t);
 else
        ireg = isregular(t);
        if ireg < regular
            switch ireg
                case 0
                    fprintf(1,'Warning: the design is sparse but has been specified as regular or regular with missing.  No computation is performed.  Please rerun with default design setting.\n');
                case 1
                    fprintf(1,'Warning: the design is regular with missing but has been specified as regular.  No computation is performed.  Please rerun with default design setting.\n');
            end
            return;
        end
 end


%Prebin the data if conditions are satisfied as specified in the help
%for numBins.
  
if isempty(numBins)
   [newy,newt] = binData(y,t,regular,verbose);
elseif numBins >= 10
   [newy,newt] = binData(y,t,regular,verbose,numBins);
elseif numBins == 0
   newy = [];   %no binning set by user
elseif numBins < 10 || numBins < 0
   newy = [];   %no binning due to number of bins is too small
   fprintf(1,'Warning: number of bins must be at least 10! No binning will be performed!\n');
end

if ~isempty(newy)
    y = newy;
    t = newt;
    if regular == 0
        regular = 1;
    end
end

if isempty(kernel)
    if regular == 2 && length(t{1}) >= 20
       kernel = 'epan'; %kernel: Epanechnikov
    else
       kernel = 'gauss';   %kernel: Gaussian
    end
else
    kernNames = {'rect','gauss','epan','gausvar','quar'};
    if isempty(strmatch(kernel, kernNames, 'exact'))
        fprintf(1,['Warning: kernel ' kernel ' is unrecognizable! Reset to default now.\n']);
        if regular == 2 && length(t{1}) >= 20
            kernel = 'epan';
        else
            kernel = 'gauss';
        end
    end
end

if isempty(method_mu)
    method_mu = 'PACE';
elseif (strcmp(method_mu, 'PACE')~=1 & strcmp(method_mu,'RARE') ~= 1)   % there  should be a  'RARE' option?
    fprintf(1,['Warning: method_mu ' method_mu ' is unrecognizable! Reset to default now.\n']);
    method_mu = 'PACE';
end

ops = struct('bwmu', bwmu,'bwmu_gcv',bwmu_gcv, 'bwxcov',bwxcov,'bwxcov_gcv', bwxcov_gcv,'ntest1',ntest1,...
	     'ngrid1',ngrid1,'selection_k',selection_k, 'FVE_threshold', FVE_threshold, 'maxk', maxk,...
	     'control',control,'regular', regular,'error',error,'ngrid',ngrid,'method', method,'shrink',shrink,...
             'newdata',newdata,'kernel',kernel, 'numBins', numBins,'yname', yname, 'screePlot', screePlot,...
             'designPlot',designPlot,'corrPlot',corrPlot,'rho', rho,'verbose',verbose);

%pool all the subjects and their corresponding time points into 1 x N vectors
%tt: vector to hold time points
%yy: vector to hold the observed measurements
%ind: indices for each subject
tt = cell2mat(t);  % 1 x N vector to hold the observed time points from all subjects
yy = cell2mat(y);  % 1 x N vector to hold the observed measurements from all subjects
%indiv= []; % 1 x N vector to hold the indices of n subjects, e.g., tt(ind == 1), yy(ind == 1) for subject 1


%Initial out1 is based on the unique time points of the pooled data + the unique
%time points of "newdata", the output time grid. When newdata = [], output
%"out1" is equivalent to be the unique sorted pooled time points; otherwise, it 
%corresponds to the unique "newdata".

out1 =unique([tt newdata]);
out21 = linspace(min(out1),max(out1),ngrid);

if designPlot == 1
    createDesignPlot(t, 0, 1, 1, yname);
end

if strcmp(verbose, 'on') == 1
  fprintf(1,'Part I: Obtain smoothed mean curve\n');
end

%when bwmu = 0 and bwmu_gcv = 0, use leave-one-curve-out CV method for bw choice
%when bwmu = 0 and bwmu_gcv = 1, use leave-one-out GCV method for bw choice
%when bwmu > 0, user-defined bw for mean function

if ~isempty(xmu) && (length(xmu) == length(out1))
    mu = xmu;
    muDense = interp1(out1,mu, out21, 'spline');
    bw_mu = [];
else    
    if bwmu == 0
        if (bwmu_gcv == 1) || (bwmu_gcv == 2)
            bw_mu = gcv_lwls(yy,tt,kernel,1,1,0,regular,verbose);   %use GCV method to choose bw for mean function
            if isempty(bw_mu)
                fprintf(1,'Error: FPCA is aborted because the observed data is too sparse to estimate the mean function!');
                mu = []; muDense=[]; no_opt =[];sigma=[];lambda=[];phi=[];eigen=[];xi_est=[];
                xi_var=[]; xcov=[];bw_xcov=[];xcovfit=[];AIC=[];BIC=[];FVE=[];
                y_pred =[];y_predOrig = [];y_predDense=[];out1 = [];out21 =[];
                rho_opt=[]; sigmanew = []; mucopy=[];phicopy=[];eigencopy=[];out1copy=[];out21copy=[];xcovcopy=[];xcovfitcopy=[];
                return;
            end
            bw_mu = adjustBW1(kernel,bw_mu,1,0,regular,verbose);
            if bwmu_gcv == 2  
                minbw = minb(tt,3); %the minimum bw for mean function 
                bw_mu = sqrt(minbw*bw_mu); %the geometric mean between the minimum bw and GCV bw
            end
            if strcmp(verbose, 'on') == 1
               fprintf(1,['Adjusted GCV bandwidth choice for mean function: ' num2str(bw_mu) '\n']);
            end
           
        else                                               %use CV method to choose bw for mean function
            bw_mu = cvfda_lwls(y,t,kernel,1,1,0,regular,verbose);
        end
    elseif bwmu > 0
        bw_mu = bwmu;
    else
        fprintf(1,'Error: Bandwidth choice for the mean function must be positive!\n');
        return;
    end
    
    %define the vector of case weight in the local weighted least square
    if weight == 1
        win1 = [];
        for i = 1:ncohort
            win1 = [win1 repmat(1/ni(i),1,ni(i))];
        end
    else
        win1 = ones(1,length(tt));
    end
    
    [invalid, mu] = lwls(bw_mu,kernel,1,1,0,tt,yy',win1,out1);
    [invalid, muDense] = lwls(bw_mu,kernel,1,1,0,tt,yy',win1,out21);
    
end

if strcmp(verbose, 'on') == 1
   fprintf(1,'Part II: Choose bandwidth of smoothing covariance surface\n');
end
rcov = getRawCov(y,t,out1,mu, regular, 0);           %obtain raw covariance;
if weight == 1
    if regular == 2
        rcov.win = repmat(1/ni(1),1,length(rcov.cxxn));
    elseif regular == 1
          ww = [];
          for j = 1:length(out1)
              Mi = length(find(rcov.tpairn(1,:) == out1(j)));
              ww = [ww repmat(1/Mi,1,Mi)];
          end
          rcov.win = ww;
    else
        ww = 1./ni;
        rcov.win = ww(rcov.indx);
    end
end

if ~isempty(xcov)&&(length(xcov)==ngrid^2)
    if length(bwxcov)>1 && bwxcov(1)>0 && bwxcov(2)>0
        bw_xcov=bwxcov;
    else
        bw_xcov = [bw_mu bw_mu];
    end
else    
    if bwxcov(1)==0 || bwxcov(2)==0
        if (bwxcov_gcv == 1) || (bwxcov_gcv == 2)
            [bw_xcov] = gcv_mullwlsn_20120823(t,ngrid1,regular,error, kernel,rcov,verbose);
            if isempty(bw_xcov)
                fprintf(1,'Error: FPCA is aborted because the observed data is too sparse to estimate the covariance function!');
                no_opt =[];sigma=[];lambda=[];phi=[];eigen=[];xi_est=[];
                xi_var=[]; xcov=[]; bw_xcov=[];xcovfit=[];AIC=[];BIC=[];FVE=[];
                y_pred =[];y_predOrig = [];y_predDense=[];out1 = [];out21 =[];
                rho_opt=[]; sigmanew = []; mucopy=[];phicopy=[];eigencopy=[];out1copy=[];out21copy=[];xcovcopy=[];xcovfitcopy=[];
                return;
            end         
            bw_xcov = adjustBW2(kernel,bw_xcov,1,0,regular,verbose);
            if bwxcov_gcv == 2  
               minbwcov = getMinb(t,unique(cell2mat(t)),regular); %the minimum bw for cov function 
               bw_xcov = sqrt(minbwcov*bw_xcov);  
            end
            if strcmp(verbose, 'on') == 1
                 fprintf(1,['Adjusted GCV bandwidth choice for COV function : (' num2str(bw_xcov(1)) ',' num2str(bw_xcov(2)) ')\n']);
            end
        else
            [bw_xcov] = cv_mullwlsn_120823(y,t,mu,ntest1,ngrid1,regular,error,kernel,rcov,verbose);
        end
    elseif bwxcov > 0
        bw_xcov = bwxcov;
    elseif bwxcov(1) < 0 || bwxcov(2) < 0
        fprintf(1,'Error: Bandwidth choice for the covariance function must be positive!\n');
        return;
    end
    
    if strcmp(verbose, 'on') == 1
        fprintf(1,'Part III: Choose number of principal components functions\n');
    end
    
    if isempty(xcov)
        AB_method = {'full','rand'};
        rcov1 = rcov;
        if error == 1
            tpairn = rcov1.tpairn;
            tneq=find(tpairn(1,:)~=tpairn(2,:));
            cyy = rcov1.cyy;
            rcov1.tpairn = tpairn(:,tneq);
            rcov1.cxxn=cyy(tneq);
            rcov1.win=rcov.win(tneq);
            if regular == 1
                rcov1.count = rcov1.count(tneq);
            end
        end
        
        if regular == 1
            [invalid,xcov]=mullwlsk(bw_xcov,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21,rcov1.count);  %smooth raw covariance;
        else
            [invalid,xcov]=mullwlsk(bw_xcov,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21);  %smooth raw covariance;
        end
    end
end    
xcov = (xcov+xcov')/2;   %transform the smoothed covariance matrix to guarantee it is a symmetric matrix.
tempvar = diag(xcov);
id = find(tempvar < 0);
if ~isempty(id)
    tempvar(id) = min(tempvar(tempvar>0));
    xcov(id,id) = min(tempvar(tempvar>0));
end
xvar = sqrt(tempvar*tempvar');
temp = find(xvar<abs(xcov));
if ~isempty(temp)
    xcov(temp) = sign(xcov(temp)).*xvar(temp);
end

if ~isempty(out_percent)
    if (out_percent > 0.25)   
        fprintf(1,['warning: leaving out ' num2str(out_percent) ' percent of the data in the boundary is too much, we reset out_percent to 0.25 !\n']);
        out_percent = 0.25;
    end
    userrange = zeros(1, 2);
    userrange(1) = quantile([tt newdata], out_percent/2);
    userrange(2) = quantile([tt newdata], 1-out_percent/2);
end

if ~isempty(userrange)
   n = length(y);
   for i = 1:n
       id = (t{i} <= userrange(2) & t{i}>=userrange(1));
       t{i} = t{i}(id);
       y{i} = y{i}(id);
   end
   out1old = out1;
   out21old = out21;
   tt = cell2mat(t);
   out1 =unique([tt newdata]);
   out21 = linspace(min(out1),max(out1),ngrid);
   mu = interp1(out1old, mu, out1);
   muDense = interp1(out21old, muDense, out21);
   xcov = interp2(out21old, out21old, xcov, out21', out21);
   xcov = (xcov+xcov')/2; 
end


if invalid == 0
  [no_opt, FVE] = no_FVE(xcov, FVE_threshold);
else
  no_opt = []; FVE = [];
  return;
end
no_optCopy = no_opt;

pc_options = {'AIC1','AIC2', 'BIC1','BIC2','FVE','user','AIC_R'};
AIC = [];
BIC = [];

if ischar(selection_k)
    k_id = strmatch(selection_k, pc_options,'exact');
    if isempty(k_id)
       fprintf(1,['Warning: Invalid method name for selection_k! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
       k_id = 5;
    end
    if k_id == 1 || k_id == 2
      [no_opt,AIC]=no_AIC(y,t,mu,bw_xcov,ngrid,regular,maxk,AB_method{k_id},method, shrink, out1, out21,kernel, error,cut, rcov,xcov);
    elseif k_id == 3 || k_id == 4
      [no_opt,BIC]=no_BIC(y,t,mu, bw_xcov, ngrid,regular,maxk,AB_method{k_id-2},method, shrink, out1, out21,kernel, error,cut, rcov,xcov);
    elseif k_id == 7
       no_opt = ngrid;
    end
    if isempty(no_opt)
       k_id = 5;
       no_opt = no_optCopy;
    end
elseif isnumeric(selection_k) && selection_k > 0
    no_opt = selection_k;
    k_id = 6;
else
    fprintf(1,['Warning: "selection_k" must be a positive integer! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
    k_id = 5;
end

if k_id ~= 7
    if strcmp(verbose, 'on') == 1
       fprintf(1, ['Best number of principal components selected by ' pc_options{k_id} ': ' num2str(no_opt) '.\n']) ;
    end
    if k_id ~= 5 
         if strcmp(verbose, 'on') == 1
            fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n']);
         end
    else
         if strcmp(verbose, 'on') == 1
           fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation (threshold = ' num2str(FVE_threshold) ').\n']);
         end
    end 

    if no_opt == maxk && k_id < 5
       fprintf(1,['Warning: ' pc_options{k_id}  ' cannot find the best No. of PC for maxk = ' num2str(maxk) '. Increase maxk to get better results.\n']);
    end

    if strcmp(verbose, 'on') == 1
       fprintf(1,['FVE calculated from ' num2str(ngrid) ' possible eigenvalues: \n']);
       disp(FVE);
    end

    if strcmp(control, 'look')
       fve2plot = FVE(1:ceil(ngrid/2))*100;
       figure;
       plot(0:length(fve2plot), [0 fve2plot],'--ro', 'LineWidth',2,...
       'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
       xlabel('\bf{No. of Principal Components}');
       ylabel('\bf{FVE (%)}');
       title(['\bf{Fraction of variance explained by No. of PC (threshold = ' num2str(FVE_threshold) ') for function ' yname '}'])
       hold on
       plot(linspace(0,no_optCopy,30),ones(30,1)*FVE(no_optCopy)*100,'b',...
       ones(30,1)*no_optCopy, linspace(0,FVE(no_optCopy)*100,30),'b');
       text(no_optCopy+0.2, FVE(no_optCopy)*100-10, {[' k = ' num2str(no_optCopy) ', FVE = ' num2str(roundoff(FVE(no_optCopy)*100,3)) '%'] ' (threshold choice)'});
       axis([0 length(FVE)+1 0 101]);
       hold off

       no_opt = input('Enter the number of principal components you want to choose:\nK=');
       if strcmp(verbose, 'on') == 1
          fprintf(1, ['You just chose ' num2str(no_opt) ' principal component(s).\n']);
          fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n\n']);
       end
    end

    %Now output scree plot based on final no_opt
    if screePlot == 1
       createSP(FVE, no_opt, yname);
    end
     
end

if strcmp(verbose, 'on') == 1
   fprintf(1,'Part IV: Perform principal components analysis\n'); 
end
		 
yLength=cellfun(@(x) length(x), y); ynew=y(1,yLength~=0);    tnew=t(1,yLength~=0);       % this line aims at removing subjects with 0 measurement induced by userrange cut.
		 

if error==1
    [invalid, sigma]=pc_covE_120823(t,out1,bw_xcov,ngrid,cut,kernel,rcov);
    if invalid==0
		
      [xi_est, xi_var, lambda, phi, eigen, no_opt, xcovfit,y_predOrig, rho_opt, sigmanew]=pc_est(ynew, tnew, mu, xcov, sigma, no_opt, error, method, shrink, out1, out21,regular,rho, verbose);
    else 
        xi_est=[]; xi_var=[]; lambda=[]; phi=[];
    end
elseif error==0
    %[invalid, xcov, out1, out21]=pc_covNE(y,t,[],mu,bw_xcov,ngrid,regular, out1,kernel, rcov,xcovn);
    sigma=[];
    if invalid==0
      [xi_est, xi_var, lambda, phi, eigen, no_opt, xcovfit,y_predOrig, rho_opt, sigmanew]=pc_est(ynew, tnew, mu, xcov, sigma, no_opt, error, method, shrink, out1, out21,regular, rho, verbose); % this line was modified by Wenwen
    else 
        xi_est=[]; xi_var=[]; lambda=[]; phi=[];
    end
end

if strcmp(method_mu,'RARE') == 1    % this part modified by Wenwen       
       
    pc_x1 = [ones(size(xi_est,1),1),xi_est];      
    lambda_y = lambda(1:no_opt);
		 [X W] = getXW(pc_x1, phi, out1, tnew, lambda_y, sigma, 1); % weighted estimation  % t=> tnew, by Wenwen, size(X)=(1,ncohort-#0LengthSubjects);
    X = (cell2mat(X))';
    W = blkdiag(W{:});
    tty = cell2mat(tnew);  % t=> tnew, by Wenwen
    Y = cell2mat(ynew)'-interp1(out1,mu,tty,'spline')'; % r.h.s y=> ynew
    b = pinv(X'*W*X)*X'*W*Y;
    b = reshape(b,no_opt,no_opt+1)';
        
    phi_YDense = eigen(:,1:no_opt);
    phi_Y=phi(:,1:no_opt);
    beta0Dense = b(1,:)*phi_YDense';
 
    beta0=b(1,:)*phi_Y';
    mu=mu+beta0;
    muDense=muDense+beta0Dense;
    
    rcov = getRawCov(ynew,tnew,out1,mu, regular, 0);           %obtain new raw covariance; 
    if weight == 1
        if regular == 2
            rcov.win = repmat(1/ni(1),1,length(rcov.cxxn));
        elseif regular == 1
              ww = [];
              for j = 1:length(out1)
                  Mi = length(find(rcov.tpairn(1,:) == out1(j)));
                  ww = [ww repmat(1/Mi,1,Mi)];
              end
              rcov.win = ww;
        else
            ww = 1./ni;
            rcov.win = ww(rcov.indx);
        end
    end
    rcov1 = rcov;
        
    if error == 1
            tpairn = rcov1.tpairn;
            tneq=find(tpairn(1,:)~=tpairn(2,:));
            cyy = rcov1.cyy;
            rcov1.tpairn = tpairn(:,tneq);
            rcov1.cxxn=cyy(tneq);
            rcov1.win=rcov.win(tneq);
            if regular == 1
                rcov1.count = rcov1.count(tneq);
            end
    end
        
    if regular == 1
            [invalid,xcov]=mullwlsk(bw_xcov,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21,rcov1.count);  %smooth raw covariance;
    else
            [invalid,xcov]=mullwlsk(bw_xcov,kernel,rcov1.tpairn,rcov1.cxxn',rcov1.win,out21,out21);  %smooth raw covariance;
    end
   
     xcov = (xcov+xcov')/2;   %transform the smoothed covariance matrix to guarantee it is a symmetric matrix.


     if invalid == 0
         [no_opt, FVE] = no_FVE(xcov, FVE_threshold);
       else
          no_opt = []; FVE = [];
       return;
     end
     
     no_optCopy = no_opt;

     pc_options = {'AIC1','AIC2', 'BIC1','BIC2','FVE','user','AIC_R'};
     AIC = [];
     BIC = [];

     if ischar(selection_k)
         k_id = strmatch(selection_k, pc_options,'exact');
         if isempty(k_id)
            fprintf(1,['Warning: Invalid method name for selection_k! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
            k_id = 5;
         end
         
         if k_id == 1 || k_id == 2
           [no_opt,AIC]=no_AIC(ynew,tnew,mu,bw_xcov,ngrid,regular,maxk,AB_method{k_id},method, shrink, out1, out21,kernel, error,cut, rcov,xcov);  % y, t=>ynew, tnew
          elseif k_id == 3 || k_id == 4
           [no_opt,BIC]=no_BIC(ynew,tnew,mu, bw_xcov, ngrid,regular,maxk,AB_method{k_id-2},method, shrink, out1, out21,kernel, error,cut, rcov,xcov);     % y, t=>ynew, tnew
          elseif k_id == 7
            no_opt = ngrid;
         end
         
         if isempty(no_opt)
            k_id = 5;
            no_opt = no_optCopy;
         end
     elseif isnumeric(selection_k) && selection_k > 0
         no_opt = selection_k;
         k_id = 6;
     else
         fprintf(1,['Warning: "selection_k" must be a positive integer! Reset to "FVE" method with threshold = ' num2str(FVE_threshold) '\n']);
         k_id = 5;
     end
     

     if k_id ~= 7
         if strcmp(verbose, 'on') == 1
            fprintf(1, ['Updating mu by RARE method and re-doing FPCA, best number of principal components selected by ' pc_options{k_id} ': ' num2str(no_opt) '.\n']) ;
         end
         if k_id ~= 5 
            if strcmp(verbose, 'on') == 1
                 fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n']);
             end
           else
              if strcmp(verbose, 'on') == 1
                fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation (threshold = ' num2str(FVE_threshold) ').\n']);
              end
          end 

         if no_opt == maxk && k_id < 5
            fprintf(1,['Warning: ' pc_options{k_id}  ' cannot find the best No. of PC for maxk = ' num2str(maxk) '. Increase maxk to get better results.\n']);
         end

        if strcmp(verbose, 'on') == 1
           fprintf(1,['Updating mu by RARE method and re-doing FPCA, FVE calculated from ' num2str(ngrid) ' possible eigenvalues: \n']);
           disp(FVE);
        end

        if strcmp(control, 'look')
           fve2plot = FVE(1:ceil(ngrid/2))*100;
           figure;
           plot(0:length(fve2plot), [0 fve2plot],'--ro', 'LineWidth',2,...
           'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize',5);
           xlabel('\bf{No. of Principal Components}');
           ylabel('\bf{FVE (%)}');
           title(['\bf{Fraction of variance explained by No. of PC (threshold = ' num2str(FVE_threshold) ') for function ' yname ' after RARE update}'])
           hold on
           plot(linspace(0,no_optCopy,30),ones(30,1)*FVE(no_optCopy)*100,'b',...
           ones(30,1)*no_optCopy, linspace(0,FVE(no_optCopy)*100,30),'b');
           text(no_optCopy+0.2, FVE(no_optCopy)*100-10, {[' k = ' num2str(no_optCopy) ', FVE = ' num2str(roundoff(FVE(no_optCopy)*100,3)) '%'] ' (threshold choice)'});
           axis([0 length(FVE)+1 0 101]);
           hold off

           no_opt = input('Enter the number of principal components you want to choose:\nK=');
           
           if strcmp(verbose, 'on') == 1
              fprintf(1, ['You just chose ' num2str(no_opt) ' principal component(s).\n']);
              fprintf(1,['It accounts for ' num2str(roundoff(FVE(no_opt), 4)*100) '%% of total variation.\n\n']);
           end
        end

        %Now output scree plot based on final no_opt
        if screePlot == 1
           createSP(FVE, no_opt, yname);
        end
     
    end



    if error==1
        [invalid, sigma]=pc_covE(tnew,out1,bw_xcov,ngrid,cut,kernel,rcov);  %t=>tnew
        if invalid==0
		    [xi_est, xi_var, lambda, phi, eigen, no_opt, xcovfit,y_predOrig, rho_opt, sigmanew]=pc_est(ynew, tnew, mu, xcov, sigma, no_opt, error, method, shrink, out1, out21,regular,rho, verbose);%y,t=>ynew,tnew;
        else 
            xi_est=[]; xi_var=[]; lambda=[]; phi=[];
        end
     elseif error==0
		sigma=[];
        if invalid==0
		    [xi_est, xi_var, lambda, phi, eigen, no_opt, xcovfit,y_predOrig, rho_opt, sigmanew]=pc_est(ynew, tnew, mu, xcov, sigma, no_opt, error, method, shrink, out1, out21,regular, rho, verbose);%y,t=>ynew,tnew;
        else 
            xi_est=[]; xi_var=[]; lambda=[]; phi=[];
        end
      end


end

%Save these copies for internal use in the FPCreg(), where the time points are guarantteed
%to be evaluated at the distinct time points of pooled t's.
mucopy = mu;
phicopy = phi;
eigencopy = eigen;
out21copy = out21;
out1copy = out1;
xcovcopy = xcov;
xcovfitcopy = xcovfit;

%Map the results to the user-defined output time grid "newdata"
if ~isempty(newdata)
   mu = interp1(out1,mu, newdata, 'spline');
   out21 = linspace(min(newdata),max(newdata),ngrid);
   
   phi = interp1(out1copy, phi, newdata, 'spline');  %eigen function based on newdata
   eigen = interp1(out21copy, eigen, out21,'spline');%eigen function based on ngrid's of newdata
   if no_opt == 1
       if size(phi,2) > 1
         phi = phi';
       end
       if size(eigen,2) > 1
         eigen = eigen';
       end
   end
   XI = repmat(out21,length(out21),1);
   XI = XI(:)';
   YI = repmat(out21,1,length(out21));
   %covariance functions based on newdata
   xcov = reshape(interp2(out21copy,out21copy,xcov,XI,YI,'spline'),length(out21),length(out21));
   xcovfit = reshape(interp2(out21copy,out21copy,xcovfit,XI,YI,'spline'),length(out21),length(out21));
   out1copy = out1;
   out1 = newdata;
end
 
 
if invalid==0 
		y_pred=NaN(ncohort,length(mu));  y_predDense=NaN(ncohort,length(muDense));  %added by Wenwen
		y_pred(yLength~=0,:) = repmat(mu,size(xi_est,1) ,1)+xi_est*phi';   %ncohort=> size(xi_est,1), ypred=> y_pred(yLength~=0,:), by Wenwen;
    y_pred = num2cell(y_pred,2)';
		y_predDense(yLength~=0,:)= repmat(muDense, size(xi_est,1),1)+xi_est*eigen';   %ncohort=>size(xi_est,1) ,y_predDense=>y_predDense(yLength~=0,:),by Wenwen;
    y_predDense = num2cell(y_predDense,2)';
end

if error==0
    sigma=[];
end

xcorr = zeros(length(out21), length(out21));
for i = 1:length(out21)
   for j = 1:length(out21)
       xcorr(i,j) = xcovfit(i,j)/sqrt(xcovfit(i,i)*xcovfit(j,j));
   end
end
%xcorr = max(min(xcorr,1),-1);
%		xi_est_copy=xi_est;xi_var_copy=xi_var;   % added by Wenwen, xi_var is a cell;
%		xi_est=NaN(ncohort,size(xi_est,2));xi_var=cell(1,ncohort); xi_est(yLength~=0,:)=xi_est_copy;xi_var(1,yLength~=0)=xi_var_copy;    % added by Wenwen

if corrPlot == 1
    if no_opt == 1
        fprintf(1,'Warning: Correlation surface is not available when only one principal component is used.');
    else        
        figure;
        mesh(out21,out21,xcorr);
        xlabel('\bf{t}');
        ylabel('\bf{t}');
        zlabel('\bf{Correlation}');
        title(['\bf{Fitted correlation surface for function ' yname '}']);
    end
end
end

function [regular] = isregular(t)
    tt = cell2mat(t);
    f = length(tt)/length(unique(tt))/length(t);
    if f==1
        regular = 2;
    elseif f>.75
        regular = 1;
    else
        regular = 0;
    end
end
 

