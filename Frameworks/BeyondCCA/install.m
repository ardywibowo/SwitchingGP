
cd algorithms;            addpath(pwd); cd ..
cd algorithms/helpers;    addpath(pwd); cd ../..
cd data;                  addpath(pwd); cd ..
cd evaluation;            addpath(pwd); cd ..
cd helpers;               addpath(pwd); cd ..
cd sample;                addpath(pwd); cd ..
cd scripts;               addpath(pwd); cd ..
cd scripts/plots;         addpath(pwd); cd ../..

try
  cd algorithms/helpers; make_jd;   cd ../..
catch
  disp('Someting went wrong with mex for jd...')
  cd ../..
end

try 
  cd algorithms/helpers; make_nojd; cd ../..
catch
  disp('Something went wrong with mex for nojd...')
  cd ../..
end
