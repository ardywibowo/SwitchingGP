function [Xtildes, UiTildes, A, Bjs, C, sigma2, sTrans, pStGivenX1Ts, logLikelihoods, initD, initPStGivenX1Ts] = inputGroupSARLearn(Xis, Uis, Lx, Lu, S, Tskip, maxIterations, imputeParams, PStGivenX1TInit, Dinit)
%INPUTGROUPSARLEARN EM training of a Switching AR Population model
% [Xtildes, UiTilde, A, Bjs, C, sigma2, sTrans, pStGivenX1Ts, logLikelihoods] = inputGroupSARLearn(Xis, Uis, Lx, Lu, S, Tskip, maxIterations)
%
% Inputs:
% Xis               : Cell vector. Each cell contains state timepoints for a different
%                     users in row vector.
% Uis               : Cell vector. Each cell contains input timepoints for a different
%                     user. Each user has different cells corresponding to different
%                     input variables in row vector.
% Lx and Lu         : Order of the AR model for state and inputs respectively.
% S                 : Number of AR models.
% Tskip             : Forces the switches to make a transition only at times t for mod(t,Tskip)==0
% maxIterations     : Maximum number of iterations
% 
% Outputs:
% Xtildes           : Imputed values of Xis
% UiTildes          : Imputed values of Uis
% A                 : Learned AR coefficients for x. Shared between each
%                     subject.
% Bjs               : Cell vector. Each cell contains learned AR coefficients for
%                     different input types Uis, shared between each subject.
% C                 : Column vector of constant coefficients
% The AR coefficients are in reverse order X(t-L) : X(t-1)
% Each column represents a different switching state
% sigma2            : Learned innovation noise
% sTrans            : Learned transition distribution
% pStGivenX1Ts      : Smoothed posterior p(S(t)|X(1:T))
% logLikelihoods    : Loglikelihood measurements for each EM loop

%% Notation:
% Mostly in camelcase except when referring to math notation
% When referring to math, variable names can be read as if it were math notation
% i: Index for user
% j: Index for input variable type
% Prefix p refers to probability
% Suffix s refers to the collection of objects

import brml.*
numSubjects = length(Xis);
numInputs = length(Uis{1});
numCoefs = Lx + Lu * numInputs + 1;

%% Initialization

D = condp(randn(numCoefs, S));
if nargin >= 10
	disp('Oh');
	D = Dinit;
end

initD = D;

A = D(1 : Lx, :);
Bjs = cell(numInputs, 1);
for j = 1 : numInputs
	Bjs{j} = D(Lx + 1 + (j-1)*Lu : Lx + j*Lu, :);
end
C = D(Lx + numInputs*Lu + 1 : end, :);

% Initialize imputed values
hasMissing = zeros(numSubjects, 1);
Xtildes = Xis;
UiTildes = Uis;
Xmissings = cell(numSubjects, 1);
UiMissings = cell(numSubjects, 1);

% Initialize transition probabilities
sTrans          = cell(numSubjects, 1);
sPriors         = cell(numSubjects, 1);
pStGivenX1Ts    = cell(numSubjects, 1);
pStStpGivenX1Ts = cell(numSubjects, 1);
sigma2          = cell(numSubjects, 1);

for i = 1 : numSubjects
	% Check for missing values
	Xi = Xis{i};
	Ui = Uis{i};
	
	Ti = length(Xi);
	UiMatrix = zeros(Ti, numInputs);
	for j = 1 : numInputs
		UiMatrix(:, j) = Ui{j};
	end
	nanUi = isnan(UiMatrix);
	UiMissings{i} = nanUi;
	
	nanXi = isnan(Xi);
	Xmissings{i} = nanXi;
	hasMissing(i) = sum(nanUi(:)) + sum(nanXi) > 0;
	
	% Initialize switches randomly
	pStGivenX1T = zeros(S, Ti);
	for t = 1 : Ti
		if mod(t-1, Tskip) == 0 || Tskip == 0
			randomState = randi([1 S]);
			pStGivenX1T(randomState, t) = 1;
		else
			pStGivenX1T(:, t) = pStGivenX1T(:, t-1);
		end
	end
	pStGivenX1Ts{i} = pStGivenX1T;
	
	% Initial Imputation
	if hasMissing(i)		
		[Xtildes{i}, UiTildes{i}] = imputeMean(Xi, Ui);
	end
	
	% Initialize switching distributions
	sTrans{i} = condp(ones(S, S)); % switch transition
	sPriors{i} = condp(ones(S, 1)); % switch prior
	
	% Initialize variances
	sigma2{i} = var(Xtildes{i}) * ones(1, S); % Separate variances
end

if nargin >= 9
	disp('Yea');
	pStGivenX1Ts = PStGivenX1TInit;	
end
initPStGivenX1Ts = pStGivenX1Ts;

% Initialize plotting parameter
logLikelihoods = zeros(maxIterations, 1);

%% Expectation Maximization Loop

for emLoopCount = 1 : maxIterations
	%% Perform expectation separately
	for i = 1 : numSubjects
		Xi = Xtildes{i};
		Ui = UiTildes{i};
		
		% Inference using HMM structure:
		[logAlpha, logLikelihood] = HMMforwardInputSAR(Xi, Ui, sTrans{i}, sPriors{i}, A, Bjs, C, sigma2{i}, Tskip);			
		logBeta = HMMbackwardInputSAR(Xi, Ui, sTrans{i}, A, Bjs, C, sigma2{i}, Tskip);
		[pStGivenX1T, pStStpGivenX1T] = HMMsmoothInputSAR(logAlpha, logBeta, A, Bjs, C, sigma2{i}, sTrans{i}, Xi, Ui, Tskip);
		
		logLikelihoods(emLoopCount) = logLikelihoods(emLoopCount) + logLikelihood;
		
		pStGivenX1Ts{i} = pStGivenX1T;
		pStStpGivenX1Ts{i} = pStStpGivenX1T;
	end
	
	%% Perform maximization
	% For each switching state
	for s = 1 : S
		SumXiViHat = zeros(numCoefs, 1); 
		SumViHatViHat = zeros(numCoefs);
		
		% For each user
		for i = 1 : numSubjects
			Xi = Xtildes{i};
			Ui = UiTildes{i};
			
			pStGivenX1T = pStGivenX1Ts{i};
			
			sigmaSum = 0; 
			sigmaNum = 0;
			
			% For each time point
			Ti = length(Xi);
			for t = 1 : Ti
				[mean, XiHat, UiHat] = sarMean(Xi, Ui, t, A, Bjs, C);
				mean = mean(s);
				ViHat = [XiHat; UiHat; 1];
				
				SumXiViHat = SumXiViHat + pStGivenX1T(s, t) * Xi(t) * ViHat ./ sigma2{i}(s);
				SumViHatViHat = SumViHatViHat + pStGivenX1T(s, t) * (ViHat * ViHat') ./ sigma2{i}(s);
				sigmaSum = sigmaSum + pStGivenX1T(s, t) * (Xi(t) - mean).^2;
				sigmaNum = sigmaNum + pStGivenX1T(s, t);
			end
			
			% Update covariance
			sigma2{i}(s) = sigmaSum / sigmaNum;
		end
		
		% Update AR Coefficients
		D = pinv(SumViHatViHat) * SumXiViHat;
		A(:, s) = D(1 : Lx);
		for j = 1 : numInputs
			Bjs{j}(:, s) = D(Lx+1 + (j-1)*Lu : Lx + j*Lu);
		end
		C(:, s) = D(end);
	end
	
	for i = 1 : numSubjects
		pStStpGivenX1T = pStStpGivenX1Ts{i};
		
		Xi = Xtildes{i};
		Ti = length(Xi);
		t = 1 : Ti-1;
		if Tskip ~= 0
			t = t(mod(t+1, Tskip) == 0);
		end
		
		% Update transition distribution
		sTrans{i} = condp(sum(pStStpGivenX1T(:, :, t), 3)');
	end
	
	%% Impute missing values on x and u for each user
	for i = 1 : numSubjects
		if hasMissing(i)
			oldXi = Xtildes{i};
			oldUi = UiTildes{i};
			pStGivenX1T = pStGivenX1Ts{i};
			
			[Xi, Ui] = imputeLeastSquares(oldXi, oldUi, Xmissings{i}, UiMissings{i}, A, Bjs, C, Lx, Lu, pStGivenX1T, imputeParams);
			
			Xtildes{i} = Xi;
			UiTildes{i} = Ui;
		end
	end
	
	disp(emLoopCount);
end

end