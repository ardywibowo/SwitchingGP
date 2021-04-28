clear model_new


% model_new = SetSelectedModel(model, selected_set);
%[a, b] = ContextPredictionValidation(model_new, signals_train, context_train, pred_lag, max_duration)
%crossvalidation
%train the model and store the trained model, for example 5 models
%and then calculate the correctness for each selected set
%[c, d] = ContextPredictionValidation(model_new, signals_test, context_test, pred_lag, max_duration)

%% Begin feature selection

%initial a feature cell
feature_set_num = 5;
feature_cell = cell(1, feature_set_num);
for m = 1:feature_set_num
    feature_cell{m} = [2*m-1, 2*m];
end
selected_feature_set_num = 2;
pred_lag  = 5;
max_duration = 5;

selected_set_index = FeatureSetSelect(model, feature_cell, selected_feature_set_num, pred_lag, max_duration);


function selected_set_index = FeatureSetSelect(model, feature_cell, selected_feature_set_num, pred_lag, max_duration)%for cross validation, substitute model by "model cell"
%The function select feature sets among a group of feature sets, with
%greedy method. 
%feature_cell is a set of feature set
%model, trained by training data
%selected_feature_set_num, how many feature set we want to select
%pred_lag, max_duration are context prediction related parameters

%selected_set_index includes the indexes of feature_cell

foldnum = 1;


selected_set_index = [];
remain_set_index = 1:feature_set_num;

for n = 1:selected_feature_set_num
    correct_cumulate = zeros(1, feature_set_num);
    for j = 1:foldnum
        for m = remain_set_index

            selected_set_union = Cellindex2Unionset([selected_set_index, m], feature_cell);
            %model_new = SetSelectedModel(model{j}, selected_set_union);
            model_new = SubsetModel(model, selected_set_union);
            %correct_context = ContextPredictionValidation(model_new, signals_test{j}, context_test{j}, pred_lag, max_duration);%%update for crossvalidation
            correct_context = ContextPredictionValidation(model_new,...
                signals_test, context_test, pred_lag, max_duration)
            correct_cumulate(m) = correct_cumulate(m)+correct_context;
        end
    end
%     update the remain set index
    [~, max_index] = max(correct_cumulate);
    remain_set_index = setdiff(remain_set_index, max_index);
    selected_set_index = union(selected_set_index, max_index);
end

end


function set = Cellindex2Unionset(index_set, feature_cell)

set = [];

for m = index_set
    set = union(set, feature_cell{m});
end

end


    
 
function model_new = SubsetModel(model, selected_set)
%selected_set:       is a vector of selected feature index
%model:              is the trined model

%model_new:          modified from model by attributes model.gp_model and
%model.feature_model
    model_new = model;


    %feature model initial
    feature_model = model.feature_model;
    pca_coeffs = feature_model.pca_coeffs;
    pca_coeffs_new = pca_coeffs;

    %gp model initial
    gp_model = model.gp_model;
    gp_coeffs = gp_model.coeffs;
    gp_model_new = model_new.gp_model;
    gp_coeffs_new = gp_coeffs;
    
    
    
    num_context = length(gp_coeffs);
    rank_approx = gp_model.rank_approx;
    num_features = rank_approx;
    nlf = rank_approx * (2*num_features - rank_approx +1) / 2; 
    
    for j = 1:num_context
        log_theta = gp_coeffs{j};
        vlf = log_theta(1:nlf);
        Lf = vec2lowtri_inchol(vlf, num_features, rank_approx);
        Kf = Lf*Lf';
        Kf_new = Kf(selected_set, selected_set);
        Lf_new = chol(Kf_new, 'lower');
        vlf_new = lowtri2vec_inchol(Lf_new, length(selected_set), length(selected_set));
        log_theta_new = [vlf_new; log_theta(nlf+1:end)];
        gp_coeffs_new{j} = log_theta_new;

        %pca feature update
        pca_coeffs_new{j} = pca_coeffs{j}(:, selected_set);
    end
    
    %modify model_new
    gp_model_new.coeffs = gp_coeffs_new;
    gp_model_new.rank_approx = length(selected_set);
    model_new.gp_model = gp_model_new;
    
    feature_model.pca_coeffs = pca_coeffs_new;
    model_new.feature_model = feature_model;
    
end



function [correct_context, total_context] = ContextPredictionValidation(model, signals_test, context_test, pred_lag, max_duration)
%Input trained model and test data set. Data structure shown in
%runContestPrediction.m

%return the correctness ratio 

total_context = 0;
correct_context = 0;


for i = 1 : length(signals_test)
	signal = signals_test{i};
	time = 1 : size(signal, 1);
	true_context = context_test{i};
	
	[predicted_context, ~] = predictContext(signal, time, model, pred_lag, max_duration, i, true_context);
	total_context = total_context + length(true_context);
    correct_context = sum(predicted_context == true_context);
end
end