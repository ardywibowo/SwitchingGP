function [ys, ys_L, ys_U] = convert_to_errorbar_format(errors, len_x, nalg)
  means = zeros(len_x,nalg);
  maxs  = zeros(len_x,nalg);
  mins  = zeros(len_x,nalg);
  for i = 1:len_x
    yysi = errors{i};
    for ialg = 1:nalg
      means(i,ialg) = mean(yysi(:,ialg));
      maxs (i,ialg) = max (yysi(:,ialg));
      mins (i,ialg) = min (yysi(:,ialg));
    end
  end
  ys   = means;
  ys_L = means - mins;
  ys_U = maxs - means;
end
