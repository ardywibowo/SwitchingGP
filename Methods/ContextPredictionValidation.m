function [correct_context, total_context] = ContextPredictionValidation(model, signals_test, context_test, pred_lag, max_duration)
%Input trained model and test data set. Data structure shown in
%runContestPrediction.m

%return the correct detected context number and total context number

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