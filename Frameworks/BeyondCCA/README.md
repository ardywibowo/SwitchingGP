### Beyond CCA: Moment Matching for Multi-View Models

####About

This project contains implementation of several moment matching based algorithms for the 
estimation in the multi-view models closely related to canonical correlation analysis.
Moreover, this code reproduces the experiments from the following paper where these 
models and algorithms where introduced:
* A. Podosinnikova, F. Bach, S. Lacoste-Julien. [Beyond CCA: Moment Matching for Multi-View Models](https://arxiv.org/abs/1602.09013). ICML, 2016.

Please cite this paper if you use this code for your research.


These algorithms designed for the situations, where classical canonical correlation analysis
(CCA) is normally used, but are better adapted to handling, e.g. non-negative count data
or a combination of non-negative count data and continuous data in different views. The proposed algorithms have good identifiability properties, as opposed to classical CCA, and therefore are better suitabile for interpretation purposes.

The code is mostly written in Matlab, but there are several C++ functions called through mex-files. Even if your Matlab does not recognize a C++ compiler, compiled mex files for Mac 
and Linux are provided.


####Quick start
* make sure your Matlab recognizes a C++ compiler: ```mex -setup```
* save all required paths and build mex files: ```install.m```
* reproduce experiments: ```reproduce_figure_*```
* synthetic data and plots are saved to the ```experiments/``` directory
* when finished, remove all paths with ```deinstall.m```

The functions to run the estimation algorithms proposed in this paper can be found in 
the ```algorithms/``` directory:
* ```dcca_gencov.m``` - Discrete CCA with generalized covariance matrices
* ```dcca.m``` - Discrete CCA with order-3 cumulants-based tensors
* ```ncca.m``` - Non-Gaussian CCA with generalized covariance matrices


Note that **discrete CCA** is designed for the cases when both views of the data are non-negative count data; **non-Gaussian CCA** is designed for the cases when both views of the data are continuous (i.e. the setting of the classical CCA, but the new model is better suited for interpretability due to the good identifiability properties); **mixed CCA** is designed for the cases when one view of the data is non-negative count data and the other one is continuous (is currently not included in this package, please contact me if you need this code).




####Questions and feedback

I would be grateful for any feedback: success and failure stories of the methods, bugs, suggestions.. I would also be happy to answer questions regarding the code and related models, so don not hesitate to contact me. My name is Anastasia Podosinnikova and you can find my current contact details in my Github profile or on my webpage.
