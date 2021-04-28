% ============
% Description:
% ============
% 
% This is the primary function for performing Functional Linear Regression, where the
% predictor is a function X(t_x) and the response can be either a function
% Y(t_y) or a scalar.
% 
% Reference: Yao, F., M\"{u}ller, H.G., Wang, J.L. (2005). Functional
% linear regression analysis for longitudinal data.  Annals of
% Statistics 33, 2873--2903. (Reference [3])
% 
% This function includes the following steps:
% 
% 1) FPCA using the PACE method for X(t_x) and/or Y(t_y) (if the responses are functions/trajectories)
% 
% 2) Computation of the regression parameter function (regression parameter surface if the 
% responses are functions)
% 
% 3) Prediction of responses or response functions, given predictor values
% 
% 4) Computation of the values of functional R-square and Quasi R-square (a simple
%    predictive measure based on the concept of proportion of variance explained
%    by the regression).
% 
% ========
% Usage:
% ========
% 
% [res, xx, yy] = FPCreg(x, t_x, y, t_y,  param_X, param_Y, FIT, K_x, K_y, isNewSub, alpha)
% 
% =======
% Input:
% =======
% 
% x    :   1*n cell array for predictor function x, where x{i} is the row vector of
%          measurements for the ith subject, i=1,...,n. It may contain data for subjects
%          that are used for prediction; this is controlled by "isNewSub" which is either a vector consisting of
%          0's and 1's according to whether subject is used for prediction (0) or estimation (1), or is 
%          a positive integer nn. In this case, nn is the number of subjects to be used for estimation and 
%          n-nn is the number of remaining subjects to be used for prediction, corresponding to 
%          the last n-nn data rows. When "isNewSub" is set to [], all n subjects
%          are used for estimation and no prediction will be calculated; see "isNewSub" for details.
% 
% 
% t_x  :   1*n cell array, t_x{i} is the row vector of time points for the ith
%          subject at which corresponding measurements x{i} are taken,
%          i=1,...,n. It contains subjects that are used for prediction.
%          See above for 2 different cases of "isNewSub" and the definition of
%          "isNewSub" for more details.
% 
% y    :   i) When no prediction is requested, that is, isNewSub = [] or
%             isNewSub = 0:
%             1*n cell array for response function y, where y{i} is the response
%             row vector for the ith subject, i = 1,...,n, 
%             or a 1*n vector for scalar response y, in which case y(i) is the response value
%             for the ith subject, i = 1,...,n.
% 
%         ii) When prediction is requested, that is isNewSub is either a vector of
%             0's (subjects used for estimation) and 1's (subjects used for prediction but not estimation)
%             or a positive integer:
%             Cell array of dim 1*nn (functional response case) or a vector
%             of dim 1*nn (scalar response case), in either case, nn is the number of subjects 
%             used for estimation.
% 
%        iii) When prediction is requested and y is still 1*n, it will be
%             truncated to 1*nn according to the isNewSub definition below.
% 
%          (see "isNewSub" for more details).
% 
% t_y  :   1*n cell array, where t_y{i} is the row vector of time points for the ith
%          subject at which corresponding measurements y{i} are taken (for each response trajectory),
%          i = 1,...,n, or [] if y is a scalar response.  See above
% 	 and the definition of "isNewSub" for more details.
% 
% param_X: an object that is returned by setOptions(), that sets the input
%          arguments for FPCA() of the X (predictor) functions (default setting is param_X = []).
%          The default method for choosing the number of principal components
%          is 'AIC_R', which uses the regression AIC criterion, selecting the number
%          of principal components based on the linear relationship between predictor
%          and response (if FIT = -1, the number of principal
%          components of X and Y are jointly selected by 'AIC_R'). 
%          If the primary target of the analysis is the estimation of the regression parameter function or surface,
%          suggested settings are FIT = 0 and 'BIC1' for 'selection_k' in param_X.
%          For optimal prediction, it is better to use the default options.
%          For other default values, see setOptions().
% 
% param_Y: an object that is returned by setOptions(), that sets the input
%          arguments for FPCA() of the Y (response) functions (default setting is param_Y = []).
%          The default method of choosing the number of principal
%          components for the response trajectories is 'BIC1'. See setOptions() for more details.
%          When y is a scalar, this object will be ignored.
% 
%  FIT:    an indicator with values 0 or -1:
%             FIT =  0: Fitting a functional linear regression (default) through decomposition into 
%                       simple linear regressions between the FPC scores of Y (response) and X 
%                       (predictor). This approach is denoted as 'YPC' in the following.
%             FIT = -1: Prediction with a functional linear regression
%                       between the actual response observations Y and the FPC scores 
%                       of X, using weighted** least squares with the covariance surface of Y as the
%                       weight matrix, this approach is denoted as 'YCO' in the following.
%             Note: Option FIT = -1 is recommended for optimal
%                prediction in the case of sparsely observed functional responses (i.e.
%                Y is sparse and functional), but it may not yield a good parameter surface estimate. 
% 
% K_x:     positive integer, normally not needed. Number of principal components
%          of functional predictor X used in functional regression, must be smaller 
%          than or equal to the maximum number of principal components available from FPCA. 
%          Default K_x = [], which means that K_x is chosen as the number of components 
%          selected by FPCA.m for predictor functions.
% 
% K_y:     positive integer, normally not needed. Has no effect for the
%          case of a scalar response. The number of principal components for the response Y,
%          used for functional regression, must be smaller than or equal to the maximum
%          number of principal components available from FPCA. Default K_y = [],
%          in which case K_y is chosen as the number of components selected by FPCA.m 
%          for the response functions.
% 
% isNewSub: i) 1*n vector of 0s or 1s, where n is the total count of subjects and 
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
%                  The count is nn for subjects with isNewSub = 0.
% 
%           This option is convenient for computing leave-one-out
%           prediction if desired.
% 
% 
%          ii) If it is a positive integer, say 5, then the last 5 subjects (in the order from
%              top to bottom within the array x) and their values for
%              t_x, t_y are used for prediction. In other words, when choosing
%              this option, one would append the ``new'' subjects for which one desires 
%              prediction of the response to occur at the end of x, t_x and t_y. Then the first 
%              nn rows of the arrays x, t_x,t_y and the entire array y (of length nn) 
%              are used for estimation and the remainder (last n-nn rows) 
%              of x,t_x and t_y (of length n-nn) will be used for prediction.
%              Example: isNewSub = 5; n = length(x)
%                       x(1:(n-5)), t_x(1:(n-5)), y, t_y(1:(n-5)) are used for estimation.
%                       x((n-4):n), t_x((n-4):n) and t_y((n-4):n) are used for prediction.
% .
%           This option is convenient for obtaining predictions for a set of new subjects.
% 
%           iii) set isNewSub = [] for the default value, which is no prediction.
% 
%        Note:
% 
%        o   nn > n - nn, the number of subjects used for estimation, should generally be larger than
%            the number of subjects used for prediction.
% 
%        o   When no prediction is requested, x,t_x, y, t_y will be of
%            length n = nn.
% 
%        o   When prediction is requested, x, t_x, t_y will be of length n,
%            which includes nn subjects for estimation and n-nn subjects for
%            prediction. Here y is always of length nn (< n).
% 
%        o   When prediction is requested and y is of length n < nn, a warning
%            message will be given and y will be truncated to length nn. This
%            assumes that only nn among the n available data for y are used for
%            estimation and the remaining n-nn observations will be ignored.
%
% alpha:   Optional, the level of the confidence bands.  alpha = 0.05 if is left empty. 
%          No confidence bands will be produced if the alpha is left out or not within the range of (0, 1).
%
% Details: i) There are no default values for the first 3 arguments, that
%             is, they are always part of the input arguments. Note, t_y can be
%             set to [] when y is a scalar response;
%         ii) Any unspecified  or optional arguments can be set to "[]" for
%             default values;
%        iii) FPCA() calls PCA(), so setOptions() sets the input
%             arguments for PCA() and the returned object contains all
%             values returned by PCA();
%         iv) Names of objects can be retrieved through names() function i.e.,
%             names(xx), names(yy), names(res) and the actual values can be
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
% 
% =======
% Output:
% =======
%   res:  an aggregated object that contains BETA, b, newx, newtx, newy,
%         new_ty, r2, Q, fitted_y, K_x, K_y.
% 
%         1) BETA:  an object that contains beta, grid_x and grid_y.
% 
%                   i) beta: a ngrid_x*ngrid_y matrix; estimated regression parameter
%                      surface or function, evaluated at given grid points.
%                  ii) grid_x : 1*ngrid_x vector, contains ngrid of distinct time
%                               points for the x function (predictor function), where
%                               ngrid is defined in setOptions().
%                 iii) grid_y : 1*ngrid_y vector, contains ngrid of distinct time
%                               points for th y function (response function), where ngrid
%                               is defined in setOptions().
%                               when y is a scalar, grid_y is [].
%                 only calculated when FIT = 0; when FIT = -1, BETA = [].
% 
%         2)  b:  a K_x*K_y matrix of estimated coefficients for the regression
%                 parameter function in the eigenbasis representation, where K_x
%                 is the number of selected FPC for X and K_y is the number of selected
%                 FPC for Y; alternatively interpreted as the regression slope coefficients of
%                 the response FPC scores on predictor FPC scores.
% 
%         3) newx:   1*numNewSub cell array that contains measurements for new x
%                    (predictor) functions.
% 
%         4) new_tx: 1*numNewSub cell array that contains time points corresponding
%                    to newx.
% 
%         5) newy:   1*numNewSub cell array that contains predicted measurements for
%                    corresponding newx.
% 
%         6) new_ty: 1*numNewSub cell array that contains time points corresponding
%                    to newy.
% 
%         7) newcb:  when y is functional, 2*numNewSub cell array, newcb{1,i} is the
%                    lower confidence bound for the i-th new subject, while newcb{2,i} 
%                    is the upper confidence bound;
%                    when y is scalar, 2*numNewSub vector, newcb(1,i) is the
%                    lower confidence bound for the i-th new subject, while
%                    newcb(2,i) is the upper confidence bound.
%
%         8) r2:     functional R-square.
% 
%         9) Q:      Quasi R-square; a measure of the fraction of the variation
%                    of the responses that is explained by the regression.
%                    If Y is a scalar,
%                       Q = 1-sum((Y_i-Yhat_i)^2)/sum((Y_i-mean(Y))^2).
%                    If Y is a function,
%                       Q=1-sum((Y_i-Yhat_i)'*(Y_i-Yhat_i)/n_i)/sum((Y_i-mean(Y))'*(Y_i-mean(Y))/n_i).
%                    (See Reference [20] Eq. (20)).
% 
%         10) fitted_y: when y is functional, 1*nn (same length of y) cell array containing fitted
%                      measurements for corresponding y that were used for
%                      estimation; when y is scalar, it is a 1*nn vector.
%        
%         11) K_x:   number of principal components used in regression for predictor X.
%         12) K_y:   number of principal components used in regression for response Y.
% 
%  xx:   an aggregated object that contains the returned values from
%        FPCA(x,t_x,param_X),
%        see PCA() or pcaHELP.txt for more details. Use getVal() to retrieve
%        the individual values.
% 
%  yy:   an aggregated object that contains the returned values from
%        FPCA(y,t_y,param_Y) or returns "[]" if y is a scalar response
%        (see PCA() or pcaHELP.txt for more details). Use getVal() to retrieve
%        the individual values.
% 
%  See
%    o    example0.m for an example of the sparse irregular data case,
%         where both predictors and responses are functional.
%    o    example2.m for an example of the regular data case,
%         where both predictors and responses are functional.
%    o    example_scal for an example of the regular data case,
%         where the predictor is functional and the response is scalar.
% 
% See also PCA, FPCA, FPCdiag, Rtest.
% 
% Use Rtest.m to conduct a hypothesis testing of Quasi R-square Q: Q=0, based on 
% bootstrapping.
% 
% Diagnostics of the fitted functional linear regression model is available
% only for the case FIT = 0. For more information about diagnostics,
% consult the help file of FPCdiag.

function [res, xx, yy] = FPCreg(x, t_x, y, t_y, param_X, param_Y, FIT, K_x, K_y, isNewSub, alpha)

    if nargin < 11
        alpha = -1;
    elseif isempty(alpha)
        alpha = .05;
    end
    
   [x, t_x, y, t_y, isYFun, newx, new_tx, newy, new_ty, invalid] = pre(x, t_x, y, t_y, isNewSub);
   
   if invalid == 1
       return;
   end
   
   if isempty(FIT) || FIT == 0
       method = 'YPC';
   elseif FIT == -1
       method = 'YCO';
   end

   if isempty(param_X)
       param_X = setOptions('selection_k','AIC_R');   %set default for param_X
   elseif strcmp(getVal(param_X,'selection_k'),'AIC_R') == 0
       fprintf(1,'Reminder: Suggested method of choosing the number of principal components for predictor X is AIC_R!\n');
   end
   if ~isempty(K_x)
       if K_x > 0
           param_X.selection_k = K_x;
       else
           fprintf(1,'Warning: K_x must be a positive integer, reset selection_k to AIC_R now\n');
           param_Y.selection_k = 'AIC_R';
           K_x = [];
       end      
   end
       
   verbose_x = getVal(param_X,'verbose');
   if strcmp(verbose_x, 'on') == 1
       fprintf(1, 'Obtain functional object for x:\n');
   end
   xx = FPCA(x, t_x, param_X);   %perform PCA on x
   
   if isYFun == 1
       
       if isempty(param_Y)
           param_Y = setOptions('selection_k','BIC1');   %set default for param_Y
       elseif strcmp(getVal(param_Y,'selection_k'),'BIC1') == 0
           fprintf(1,'Reminder: Suggested method of choosing the number of principal components for response Y is BIC1!\n');
       end
       if ~isempty(K_y)
           if K_y > 0
               param_Y.selection_k = K_y;
           else
               fprintf(1,'Warning: K_y must be a positive integer, reset selection_k to BIC1 now\n');
               param_Y.selection_k = 'BIC1';
               K_y = [];
           end
       end

       verbose_y = getVal(param_Y,'verbose');
       if strcmp(verbose_y, 'on') == 1
           fprintf(1,'Obtain functional object for y:\n');
       end
       yy = FPCA(y, t_y, param_Y);   %perform PCA on y
       y = getVal(yy, 'y');               
       t_y = getVal(yy,'t');

   else

      yy = [];        
      verbose_y = [];
      
   end
   

   % select the numbers of principal components of X and Y included in
   % regression
   Kopt_x = getVal(xx,'no_opt');
   if isempty(K_x)
       if strcmp(getVal(param_X,'selection_k'),'AIC_R') == 1
           [K_x, K_y, xx, yy] = getNumK(xx, yy, y, t_y, param_X, param_Y, isYFun, method);
       else
           K_x = Kopt_x;          %set default value for K_x
       end
   elseif K_x ~= Kopt_x
       K_x = Kopt_x;
   end
   
   % error checking on K_y   
   if isYFun == 1
       K_y = getVal(yy,'no_opt');
   else
       if ~isempty(K_y)
           fprintf(1, 'Warning: K_y is not needed when Y is a scalar!\n');
           K_y = [];
       end
   end
  
   
   if strcmp(verbose_x, 'on') == 1 || strcmp(verbose_y, 'on') == 1
       fprintf(1,'Calculate regression function:\n');
   end
    
   [BETA,b] = getBeta(xx, yy, y, t_y, isYFun, K_x, K_y, method);
   
   if ~isempty(newx)
       [ypred, newxscore] =  FPCApred(xx, newx, new_tx);
       clear ypred;
       newxscore = newxscore(:,1:K_x);
       [newy newcb] = predict(newxscore, new_tx, xx, y, yy, new_ty, b, isYFun, K_x, K_y, method, alpha); 
   end
   
   pc_x = getVal(xx,'xi_est');
   pc_x = pc_x(:,1:K_x);
   out1_y = getVal(yy,'out1');
   nsub_y = length(y);
   if isYFun==1       
%        fitted_y = predict(pc_x,y,yy,t_y,b,isYFun,K_x,K_y,method);
       [fitted_y fitted_cb]= predict(pc_x,t_x,xx,y,yy,mat2cell(repmat(out1_y,nsub_y,1),ones(1,nsub_y))',b,isYFun,K_x,K_y,method,alpha);
       fitted_y1 = cell(1,nsub_y);
       for i = 1:nsub_y
           fitted_y1{i} = interp1(out1_y,fitted_y{i},t_y{i},'spline');
       end
       Q = getQ(fitted_y1, yy, y, t_y, isYFun);
   else
       [fitted_y fitted_cb]= predict(pc_x,t_x,xx,y,yy,[],b,isYFun,K_x,K_y,method,alpha);
       Q = getQ(fitted_y, [], y, t_y, isYFun);
       fitted_y1 = fitted_y;
   end
   
   r2 = getR2(b, xx, yy, y, K_x, K_y, isYFun);
   
   resNames = {'BETA', 'b', 'newx', 'new_tx', 'newy', 'new_ty', 'newcb', 'r2','Q', 'fitted_y','fitted_cb','K_x','K_y','method','isYFun','fitted_y1'};
   res = {BETA, b, newx, new_tx, newy, new_ty, newcb, r2,Q, fitted_y, fitted_cb, K_x,K_y, method, isYFun, fitted_y1, resNames};

end
