% ============
% Description:
% ============
% 
% This is the main function to implement Functional Quadratic Regression (QuadReg), where the
% predictor is a function X(t_x) and the response Y is a scalar. 
% 
% Reference: Yao, F. and M\"uller, H.G. (2009). Functional Quadratic Regression.
% Biometrika.
% 
% 
% This model can be used for both sparsely and densely sampled trajectories. 
% Gaussian assumptions are required for the coefficient estimates to be 
% consistent in the sparse case. Please refer to the above paper for more details.  
% 
% Includes the following steps:
% 
% 1) FPCA using the PACE method for X(t_x).
% 
% 2) Computation of the quadratic regression model components
% 
% 3) Prediction of response for given predictor values
% 
% 4) Computation of (generalized) R-squared.
% 
% ========
% Usage:
% ========
% 
% [res, xx] = FPCQuadReg(x, t_x, y, param_X, K, isNewSub)
% 
% =======
% Input:
% =======
% 
% x    :   1*n cell array for predictor function x, where x{i} is the row 
%          vector of measurements for the ith subject, i=1,...,n. It may 
%          contain data for new subjects whose response is to be predicted,
%          while not contributing to model fitting; this is controlled by 
%          "isNewSub" which is either a vector consisting of 0's and 1's, 
%          according to whether a subject is used for prediction (0) or 
%          estimation (1), or is controlled by a positive integer nn. In 
%          this case, nn is the number of subjects to be used for 
%          estimation and n-nn is the number of remaining subjects to be 
%          used for prediction, corresponding to the last n-nn data rows. 
%          When "isNewSub" is set to [], all n subjects are used for 
%          estimation and no prediction will be calculated; see "isNewSub" 
%          for more details.
% 
% t_x  :   1*n cell array, t_x{i} is the row vector of time points for the ith
%          subject at which corresponding measurements x{i} are taken,
%          i=1,...,n. It contains subjects that are used for prediction.
%          See above for 2 different cases of "isNewSub" and the definition of
%          "isNewSub" for more details.
% 
% y    :   i) When no prediction is requested, that is, isNewSub = [] or
%             isNewSub = 0:
%             1*n vector for scalar response y, then y(i) is the response value
%             for the ith subject, i = 1,...,n.
% 
%         ii) When prediction is requested, that is isNewSub is either a vector of
%             0's (subjects used for estimation) and 1's (subjects used for prediction but not estimation)
%             or a positive integer:
%             a vector of dim 1*nn (scalar response case). In either case, nn is the number of subjects 
%             used for estimation.
% 
%        iii) When prediction is requested and y is 1*n, it will be
%             truncated to 1*nn according to the isNewSub definition below.
%             (see "isNewSub" for more details).
%
% param_X: an object that is returned by setOptions() that sets the input
%          arguments for FPCA() of the X (predictor) functions (for default, set param_X = []).
%          The default method for choosing the number of principal components
%          is 'AIC_R', which uses the regression AIC criterion, selecting the number
%          of principal components based on the linear relationship between predictor
%          and response (Note that if FIT = -1, the number of principal
%          components of X and Y are jointly selected by 'AIC_R'). 
%          If the main goal of the analysis is the estimation of the regression parameter function or surface, 
%          suggested settings are FIT = 0 and 'BIC1' for 'selection_k' in param_X.
%          For optimal prediction, it is better to use the default options.
%          For other default values, see setOptions().
% 
% K:     positive integer, not required. Number of principal components
%          of predictor x used in regression, must be smaller than or equal to the maximum
%          number of principal components available from FPCA.  Default K = [], in this case
%          K is the number of components selected by FPCA.m for predictor functions.
% 
% 
% isNewSub: i) 1*n vector of 0s or 1s, where
% 
%              1 : the data for the corresponding subject, i.e.,
%                  x(isNewSub == 1), t_x(isNewSub==1), t_y(isNewSub == 1)
%                  are used for prediction only;
% 
%                 The count is n-nn for subjects with isNewSub = 1.
% 
%              0 : the data for the corresponding subject, i.e.,
%                  x(isNewSub == 0), t_x(isNewSub == 0), y and
%                  t_y(isNewSub == 0) are used for estimation only.
% 
%                  The count is nn for subjects with isNewSub = 0
% 
%           This option is convenient for computing leave-one-out
%           prediction if desired.
% 
% 
%          ii) If it is a positive integer, say 5, then the last 5 subjects (in the order from
%              top to bottom within the array x) and their values for
%              t_x, t_y are used for prediction. In other words, when choosing
%              this option, one would append the ``new'' subjects for which one desires prediction of the response
%              to occur at the end of x, t_x and t_y. Then the first nn rows of the arrays x, t_x,t_y and the entire 
%              array y (of length nn) are used for estimation and the remainder (last n-nn rows) 
%              of x,t_x and t_y (of length n-nn) will be used for prediction.
%              Example: isNewSub = 5; n = length(x)
%                       x(1:(n-5)), t_x(1:(n-5)), y, t_y(1:(n-5)) are used for estimation.
%                       x((n-4):n), t_x((n-4):n) and t_y((n-4):n) are used for prediction.
% 
%           This option is convenient for obtaining predictions for a set of new subjects.
% 
%           iii) set isNewSub = [] for the default value, which is no prediction.
% 
%        Note:
% 
%        o   nn > n - nn, the number of subjects used for estimation should generally be larger than
%            number of subjects used for prediction.
% 
%        o   When no prediction is requested, x,t_x, y, will be of
%            length n = nn.
% 
%        o   When prediction is requested, x, t_x will be of length n,
%            which includes nn subjects for estimation and n-nn subjects for
%            prediction. Here y is always of length nn (< n).
% 
%        o   When prediction is requested and y is of length n < nn, a warning
%            message will be given and y will be truncated to length nn. This
%            assumes that only nn among the n available data for y are used for
%            estimation and the remaining n-nn observations will be ignored.
% 
% Details: i) There are no default values for the first 3 arguments, that
%             is, they are always part of the input arguments.
%         ii) Any unspecified or optional arguments can be set to "[]" for
%             default values;
%        iii) FPCA() calls PCA(), so setOptions() sets the input
%             arguments for PCA() and the returned object contains all
%             values returned by PCA();
%         iv) Names of objects can be retrieved through names() function i.e.,
%             names(xx), names(res) and the actual values can be
%             retrieved through the getVal() function, example: getVal(res, 'BETA'),
%             getVal(res, 'newy') etc.
%          v) When isNewSub is set to be [], no prediction will be performed,
%             and the number of subjects in x is the same as the number of subjects
%             in y;
%             when isNewSub is either a vector of 0's and 1's or a positive
%             integer, the number of subjects in x is larger than the number of subjects
%             of y, since the former contains new subjects for prediction only,
%             whose response will not be used in the model estimation process.
%          Vi)Default: no prediction will be performed.
% 
% =======
% Output:
% =======
%   res:  an aggregated object that contains coeff, coeffFun, design, newx, new_tx, newy, 
%        Q, fitted_y, K.
% 
%         1) coeff: an object that contains alpha, beta and gamma.
%                  i)  alpha: the estimate of the intercept.
%                  ii) beta: a K_x*1 vector of estimated coefficients for the 
%                      linear term in the eigenbasis representation, where K_x
%                      is the number of included FPC components.
%                  iii) gamma: a {K*(K+1)/2}*1 vector of estimated
%                  coefficients for the quadratic terms in the eigenbasis
%                  representation. The sequence is gamma_{11}, gamma_{21}, gamma_{22},
%                   gamma_{31},... ,gamma_{KK}, as described in the
%                   reference.
%         2) coeffFun: an object that contains betafun, gammafun, and grid.
%                  i) betafun: a 1*ngrid vector; estimated regression
%                  function for the linear term, evaluated at given grid points. 
%                  ii) gammafun: a ngrid*ngrid vector; estimated regression
%                  surface for the quadratic term. 
%                  iii) grid: 1*ngrid vector, contains ngrid of distinct
%                  time points for the x function, where ngrid is defined
%                  in setOptions().
%         3) design: n*{(K+1)*(K+2)/2} design matrix for the regression model. 
%          4) newx: 1*numNewSub cell array that contains measurements for new x
%                  (predictor) functions.
% 
%         5) new_tx: 1*numNewSub cell array that contains time points corresponding
%                    to newx.
% 
%         6) newy:   1*numNewSub vector that contains predicted measurements for
%                    corresponding newx.

%         7) Q:      Quasi R-square; a measure of the fraction of the variation
%                    of the responses that is explained by the regression.
%                       Q = 1-sum((Y_i-Yhat_i)^2)/sum((Y_i-mean(Y))^2).
% 
%         9) fitted_y: 1*nn (same length of y) vector containing fitted
%                      measurements for corresponding y values that were used for
%                      estimation.
%        
%         10) K:   number of principal components used in regression for predictor X.
%    
% 
%  xx:   an aggregated object that contains the returned values from
%        FPCA(x,t_x,param_X).
%        See PCA() or pcaHELP.txt for more details. Use getVal() to retrieve
%        the individual values.
% 
% 
%    o    See exampleQuadReg.m for an example of the Functional Quadratic
%         Regression.

% See also PCA, FPCA, FPCreg

function [res, xx] = FPCQuadReg(x, t_x, y, param_X, K, isNewSub)


   [x, t_x, y, newx, new_tx, newy, invalid] = pre1(x, t_x, y, isNewSub);

   if invalid == 1
       return;
   end

   if isempty(param_X)
       param_X = setOptions('selection_k','BIC1');   %set default for param_X
   elseif strcmp(getVal(param_X,'selection_k'),'BIC1') == 0
       fprintf(1,'Reminder: Suggested method of choosing the number of principal components for predictor X is BIC1!\n');
   end
   if ~isempty(K)
       if K > 0
           param_X.selection_k = K;
       else
           fprintf(1,'Warning: K must be a positive integer, reset selection_k to BIC1 now\n');
           param_X.selection_k = 'BIC1';
       end      
   end
       
   verbose_x = getVal(param_X,'verbose');
   if strcmp(verbose_x, 'on') == 1
       fprintf(1, 'Obtain functional object for x:\n');
   end
   xx = FPCA(x, t_x, param_X);   %perform PCA on x
   K_x= getVal(xx,'no_opt');
   xi = getVal(xx,'xi_est');  
   if isempty(K)
       K=K_x;
   end
   
   n=size(xi,1);
   xi_K=xi(:,1:K);
   qdm1=[ones(n,1),qdmscores1(xi_K)];
   qcoef1=pinv(qdm1'*qdm1)*qdm1'*y';
   m=length(qcoef1);
   
   design=qdm1;
   alpha=qcoef1(1);
   beta=qcoef1(2:(K+1));
   gamma=qcoef1((K+2):m);
   coeffNames = {'alpha','beta','gamma'};
   coeff = {alpha,beta,gamma,coeffNames};
   
   phi=getVal(xx,'eigencopy');
   grid=getVal(xx,'out21');
   phi_K=phi(:,1:K);
   betafun=(phi_K*beta)';
   S=size(phi_K,1);
   gammafun=zeros(S,S);
   C_kl=0;
   for k=1:K
       for l=1:k
           C_kl=C_kl+1;
           gammafun=gammafun+gamma(C_kl)*phi_K(:,k)*phi_K(:,l)';
       end
   end
   coeffFunNames={'betafun','gammafun','grid'};
   coeffFun={betafun,gammafun,grid,coeffFunNames};
   
   if ~isempty(newx)
       [ypred, newxscore] =  FPCApred(xx, newx, new_tx);
       clear ypred;
       newxscore_K = newxscore(:,1:K);
       newy = predQ(newxscore_K,coeff); 
   end
   
       fitted_y = predQ(xi_K,coeff);
       Q = 1-sum((y-fitted_y).^2)/sum((y-mean(y)).^2);
   
   resNames={'coeff', 'coeffFun', 'design', 'newx', 'new_tx','newy', 'Q', 'fitted_y', 'K'};
   res = {coeff, coeffFun, design, newx, new_tx, newy, Q, fitted_y, K, resNames};

   end
