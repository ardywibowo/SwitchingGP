tic;
clear; close all; clc;

addpath(genpath('../'));

%% Things to try
% 1. Averaging each context * subject instance individually
% 2. Fit Gaussian Process baseline on all population
% 3. 

%% Info
% 1 WALKING
% 2 WALKING_UPSTAIRS
% 3 WALKING_DOWNSTAIRS
% 4 SITTING
% 5 STANDING
% 6 LAYING

%% Import Data
disp('Importing Data');

freq = 1000/50; % ms

signals_train = importdata('../Data/UCI HAR Dataset/train/X_train.txt');
signals_test = importdata('../Data/UCI HAR Dataset/test/X_test.txt');

subject_train = importdata('../Data/UCI HAR Dataset/train/subject_train.txt');
subject_test = importdata('../Data/UCI HAR Dataset/test/subject_test.txt');

raw_context_train = importdata('../Data/UCI HAR Dataset/train/y_train.txt');
raw_context_test = importdata('../Data/UCI HAR Dataset/test/y_test.txt');

raw_train_directory = dir('../Data/UCI HAR Dataset/train/Inertial Signals/*.txt');
raw_test_directory = dir('../Data/UCI HAR Dataset/test/Inertial Signals/*.txt');

signals_train_raw = [];
for i = 1 : length(raw_train_directory)
	curr_raw = importdata(['../Data/UCI HAR Dataset/train/Inertial Signals/', raw_train_directory(1).name]);
	signals_train_raw = cat(3, signals_train_raw, curr_raw);
end
	
signals_test_raw = [];
for i = 1 : length(raw_test_directory)
	curr_raw = importdata(['../Data/UCI HAR Dataset/test/Inertial Signals/', raw_test_directory(1).name]);
	signals_test_raw = cat(3, signals_test_raw, curr_raw);
end

clear curr_raw i raw_train_directory raw_test_directory;

% For now
clear signals_test_raw signals_train_raw;

%% Divide Data Based on Different Subjects
disp('Formatting Data by Subjects');

num_subjects = max(subject_train);
num_contexts = max(raw_context_train);

train_subjects = 1 : num_subjects;
signals_training = cell(num_subjects, 1);
context_train = cell(num_subjects, 1);
for i = 1 : length(signals_training)
	signals_training{i} = signals_train(subject_train == i, :);
	context_train{i} = raw_context_train(subject_train == i);
end

empty_cells = cellfun(@isempty, signals_training);
signals_training(empty_cells) = [];
context_train(empty_cells) = [];
train_subjects(empty_cells) = [];

test_subjects = 1 : num_subjects;
signals_testing = cell(num_subjects, 1);
context_test = cell(num_subjects, 1);
for i = 1 : length(signals_testing)
	signals_testing{i} = signals_test(subject_test == i, :);
	context_test{i} = raw_context_test(subject_test == i);
end

empty_cells = cellfun(@isempty, signals_testing);
signals_testing(empty_cells) = [];
context_test(empty_cells) = [];
test_subjects(empty_cells) = [];

signals_train = signals_training;
signals_test = signals_testing;
clear empty_cells signals_training signals_testing raw_context_train raw_context_test i subject_train subject_test;

%% Divide by contexts
disp('Dividing by Contexts');

signals_training = cell(num_contexts, 1);
signals_testing = cell(num_contexts, 1);

for j = 1 : num_contexts
	signals_training{j} = cell(length(signals_train), 1);
	signals_testing{j} = cell(length(signals_test), 1);
	
	for i = 1 : length(signals_train)
		current_context = context_train{i} == j;
		disjoints = (current_context ~= 0)';
		
		dis1 = strfind([0 disjoints], [0 1]);
		dis2 = strfind([disjoints 0], [1 0]);
		context_signals = arrayfun(@(x,y) signals_train{i}(x:y, :), dis1, dis2, 'un', 0);
		
		signals_training{j}{i} = context_signals';
	end
	
	for i = 1 : length(signals_test)
		current_context = context_test{i} == j;
		disjoints = (current_context ~= 0)';

		dis1 = strfind([0 disjoints], [0 1]);
		dis2 = strfind([disjoints 0], [1 0]);
		context_signals = arrayfun(@(x,y) signals_test{i}(x:y, :), dis1, dis2, 'un', 0);
		
		signals_testing{j}{i} = context_signals';
	end
end

signals_train = signals_training;
signals_test = signals_testing;

clear i j signals_training signals_testing dis1 dis2 disjoints context_signals current_context;

%% Format and Reduce Data Dimensionality
disp('Applying PCA');

% Timepoints
time_train = cell(num_contexts, 1);
time_test = cell(num_contexts, 1);
pca_coeffs = cell(num_contexts, 1);

signal_means = cell(num_contexts, 1);

pca_rank = 10;

% Subtract means
for j = 1 : num_contexts
	
	signal_means{j} = cell(length(signals_train{j}), 1);
	for i = 1 : length(signals_train{j})
		signal_cat = cell2mat(signals_train{j}{i});
		
		signal_means{j}{i} = mean(signal_cat);
		for k = 1 : length(signals_train{j}{i})
			signals_train{j}{i}{k} = signals_train{j}{i}{k} - signal_means{j}{i};
		end
	end
end

% Combine means
total_means = cell(num_contexts, 1);
for j = 1 : num_contexts
	all_signals = cat(1, signals_train{j}{:});
	total_means{j} = mean(cat(1, all_signals{:}));
end

for j = 1 : num_contexts
	indiv_cat = cat(1, signals_train{j}{:});
	signals_cat = cat(1, indiv_cat{:});
	[coeff, score] = pca(signals_cat, 'NumComponents', pca_rank);
	
	% TODO Check means of individuals to see if they are somewhat equal
	pca_coeffs{j} = coeff;
	
	time_train{j} = cell(length(signals_train{j}), 1);
	curr_index = 0;
	for i = 1 : length(signals_train{j})
		time_train{j}{i} = cell(length(signals_train{j}{i}), 1);
		for k = 1 : length(signals_train{j}{i})
			time_train{j}{i}{k} = 1 : size(signals_train{j}{i}{k}, 1);
			signals_train{j}{i}{k} = score(curr_index+1 : curr_index + size(signals_train{j}{i}{k}, 1), :);
			curr_index = curr_index + size(signals_train{j}{i}{k}, 1);
		end
	end
	
	time_test{j} = cell(length(signals_test{j}), 1);
	for i = 1 : length(signals_test{j})
		time_test{j}{i} = cell(length(signals_test{j}{i}), 1);
		for k = 1 : length(signals_test{j}{i})
			time_test{j}{i}{k} = 1 : size(signals_test{j}{i}{k}, 1);
		end
	end
end

feature_model = struct;
feature_model.means = total_means;
feature_model.pca_coeffs = pca_coeffs;

save('../Results/MTGP/feature_model.mat', 'feature_model');

clear i j k coeff curr_index std_x mean_x cat_x indiv_cat signal_means pca_coeffs score signal_cat signals_cat total_means signals_training signals_testing;

%% Train HMM
disp('Training HMM');

p_trans = zeros(num_contexts);
for i = 1 : length(context_train)
	subject_context = context_train{i};
	transitions = subject_context(diff([0; subject_context]) ~= 0);
	
	for t = 1 : length(transitions)-1
		p_trans(transitions(t+1), transitions(t)) = p_trans(transitions(t+1), transitions(t)) + 1;
	end
end

% Add prior and normalize
% p_trans = p_trans + 50 * ones(num_contexts);
p_trans = p_trans ./ sum(p_trans);

save('../Results/MTGP/hmm_model.mat', 'p_trans');

clear i t transitions subject_context;

% Remark: Transitions very biased and useless for real life.

%% Train Gamma Process
disp('Training Gamma Process');

gamma_model = cell(num_contexts, 1);
for j = 1 : num_contexts
	time = cat(1, time_train{j}{:});
	time_lengths = cellfun('length', time);
	
	v = log(mean(time_lengths)) -  mean(log(time_lengths));
	k_hat = (3 - v + sqrt((v-3)^2+24*v)) / (12 * v);
	beta_hat = mean(time_lengths) / k_hat;
	
	gamma_model{j} = [k_hat, beta_hat];
end
save('../Results/MTGP/gamma_model.mat', 'gamma_model');

% Remarks: Notice that walking upstairs/downstairs is smallest duration
% highest duration is walking and standing

clear j time k_hat beta_hat v;

%% Train GP
disp('Training GP');

% Format Data
signals_training = cell(num_contexts, 1);
signals_testing= cell(num_contexts, 1);
for j = 1 : num_contexts
	signals_training{j} = cell(length(signals_train{j}), 1);
	for i = 1 : length(signals_train{j})
		signals_training{j}{i} = cat(1, signals_train{j}{i}{:});
		time_train{j}{i} = 1 : size(signals_training{j}{i}, 1);
	end
	
	signals_testing{j} = cell(length(signals_test{j}), 1);
	for i = 1 : length(signals_test{j})
		signals_testing{j}{i} = cat(1, signals_test{j}{i}{:});
		time_test{j}{i} = 1 : size(signals_testing{j}{i}, 1);
	end
end

signals_train = signals_training;
signals_test = signals_testing;

clear signals_training signals_testing j i;

%% GP Parameters
% Constants
num_features = size(signals_train{1}{1}, 2);
rank_approx = num_features;
num_obs = ones(num_features, 1);

% Parameters
num_iter = 200;
cov_baseline = {'covSEard'};
cov_contexts = {'covMatern5iso'};

%% Begin Training

task_index = cell(num_contexts, 1);
time_index = cell(num_contexts, 1);
signals_train_in = cell(num_contexts, 1);

for j = 1 : num_contexts
	task_index{j} = cell(length(signals_train{j}), 1);
	time_index{j} = cell(length(signals_train{j}), 1);
	signals_train_in{j} = cell(length(signals_train{j}), 1);
	for i = 1 : length(signals_train{j})
		task_index{j}{i} = kron(1 : size(signals_train{j}{i}, 2), ones(1, size(signals_train{j}{i}, 1)))';
		time_index{j}{i} = repmat(1 : size(signals_train{j}{i}, 1), 1, size(signals_train{j}{i}, 2))';
		
		time_train{j}{i} = time_train{j}{i}';
		
		signals_train_in{j}{i} = signals_train{j}{i}(:);
	end
end

cov_func = [cov_baseline, cov_contexts];
data = {cov_func, time_train, signals_train_in, num_features, rank_approx, num_obs, task_index, time_index};
[log_theta, deriv_range] = init_jointmtgp_default(signals_train_in, cov_func, num_features, rank_approx);
[log_theta, nl] = learn_jointmtgp(log_theta, deriv_range, data, num_iter);

save('../Results/MTGP/gp_model.mat', 'log_theta');


%% Plot and Test Model
disp('Plotting Results');

dataRatio = 4;

Ypred_all = cell(num_contexts, 1);
Vpred_all = cell(num_contexts, 1);
differences = cell(num_contexts, 1);
for j = 1 : num_contexts
	time = time_train{j};
	
	% Predict for each subject
	Ypred_all{j} = cell(length(time), 1);
	Vpred_all{j} = cell(length(time), 1);
	differences{j} = cell(length(time), 1);
	for i = 1 : length(time)
		% True signal
		time_pred = time{i};
		signal_true = signals_train{j}{i}(time_pred, :);
		
		% Observed signal
		time_obs = downsample(time{i}, dataRatio);
		signal_obs = signals_train{j}{i}(time_obs, :);
		
		% Format observed signal for prediction
		task_index_obs = kron(1 : size(signal_obs, 2), ones(1, size(signal_obs, 1)))';
		time_index_obs = repmat(1 : size(signal_obs, 1), 1, size(signal_obs, 2))';
		signal_obs_in = signal_obs(:);
		
		% Predict trajectory
		D = size(time{1}, 2);
		
		cov_pred = {'covSum', {'covSEard', 'covMatern5iso'}};
		length_Lf = rank_approx * (2*num_features - rank_approx +1) / 2;
		length_baseline = eval(feval(cov_func{1}));
		length_contexts = eval(feval(cov_func{2}));
		
		theta_lf = log_theta(1 : length_Lf);
		theta_baseline = log_theta(length_Lf+1 : length_Lf+length_baseline);
		theta_contexts = log_theta(length_Lf+length_baseline+length_baseline + (j-1)*length_contexts + 1 : ...
			length_Lf+length_baseline+length_baseline + j*length_contexts);
		theta_noise = log_theta(length_Lf+length_baseline + num_contexts*length_contexts + 1 : end);
				
		pred_coeffs = [theta_lf; theta_baseline; theta_contexts; theta_noise];
		data_obs = {cov_pred, time_obs, signal_obs_in, num_features, rank_approx, num_obs, task_index_obs, time_index_obs};
		[signal_pred, signal_variance] = predict_popmtgp_all_tasks(pred_coeffs, data_obs, time_pred);
		
		% Prediction and confidence itervals
		Vpred = signal_variance;
		Ypred_all{j}{i} = signal_pred;
		Vpred_all{j}{i} = signal_variance;
		
		Ypred_up = signal_pred + 1.96 * sqrt(Vpred);
		Ypred_down = signal_pred - 1.96 * sqrt(Vpred);
		
		% Compute Differences
		difference = abs(signal_pred - signal_true);
		differences{j}{i} = difference;
		
		% Plot results
		figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
		suptitle(['Predictions Context ' num2str(j), ' Trial ' num2str(i)]);
		for p = 1 : num_features
			subplot(5, 2, p);
			hold on; box on;
			
			% Prediction and confidence intervals
			conf_color = [0.95 0.95 0.95];
			ciplot(Ypred_up(:, p), Ypred_down(:, p), time_pred * freq, conf_color);
			plot(time_pred * freq, signal_pred(:, p));
			
			% Data points
			scatter(time_obs * freq, signal_obs(:, p), 60, 'd', 'red', 'filled');
			scatter(time_pred * freq, signal_true(:, p), 'o', 'blue');
			
			xlabel('Time (ms)'), ylabel(['Feature ' num2str(p)]);
		end
		saveFigure('../Results/MTGP/', ['Context ' num2str(j), ' Trial ' num2str(i)]);
	end
end

differences_all = cat(1, differences{:});
differences_all = cat(1, differences_all{:});
mse_all = mean(differences_all(:).^2);
abs_all = mean(differences_all(:));

save('../Results/MTGP/errors.mat', 'mse_all', 'abs_all');
save('../Results/MTGP/predictions.mat', 'Ypred_all', 'Vpred_all');


% Predict with multiple steps


toc;

%% Matern Stuff

% % Train each context separately
% 
% % Distance matrix
% dist_train = cell(num_contexts, 1);
% dist_test = cell(num_contexts, 1);
% for j = 1 : num_contexts
% 	dist_train{j} = {};
% 	for i = 1 : length(time_train{j})
% 		dist_train{j}{i} = squareform(pdist(time_train{j}{i}'));
% 	end
% 	
% 	dist_test{j} = {};
% 	for i = 1 : length(time_test{j})
% 		dist_test{j}{i} = squareform(pdist(time_test{j}{i}'));
% 	end
% end

% function ell = ell_matern(dist, x, range, nu, sigma_l)
% 
% cov_matern = matern_cov(dist, range, nu);
% 
% 
% 
% cov_all
% 
% sigma_c = chol(cov_matern, 'lower')';
% 
% out = sigma_c \ x;
% 
% quad_form = sum(out.^2);
% det_part = 2 * sum(log(diag(sigma_c)));
% 
% ell = 0.5*det_part + 0.5*quad_form;
% 
% end

function cov_mat = matern_cov(dist, range, nu)

num_samples = length(dist);
cov_m = cell(num_samples, 1);
for i = 1 : num_samples
	cov_m{i} = matern(dist{i}, range, nu);
end
cov_mat = spdiags(blkdiag(cov_m{:}));

end







