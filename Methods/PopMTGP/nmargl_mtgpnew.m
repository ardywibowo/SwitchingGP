function [total_nl, total_gradnl] = nmargl_mtgpnew(logtheta, logtheta_all, cov_func, times, signals, ...
				      num_features, rank_approx, task_index, time_index, deriv_range)
%NMARGL_MTGPNEW Summary of this function goes here
%   Detailed explanation goes here

total_nl = 0;
total_gradnl = 0;
for i = 1 : length(signals)
	signal = signals{i};
	time = times{i};
	num_obs = ones(length(signal), 1);
	
	[nl, gradnl] = nmargl_mtgp(logtheta, logtheta_all, cov_func, time, signal, num_features, rank_approx, num_obs, task_index{i}, time_index{i}, deriv_range);
	total_nl = total_nl + nl;
	total_gradnl = total_gradnl + gradnl;
end

end

