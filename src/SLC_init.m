function conf = SLC_init(conf)
%%  SLC Sparse Lapacian Consensus
%
%   CONF = SLC_INIT(CONF) sets the default configuration for SLC.
%
%   gamma: Percentage of inliers in the samples. This is an inital value
%       for EM iteration, and it is not important. Default value is 0.9.
%
%   delta: Paramerter of Edge Weight, W(x, y) = exp(-||x-y||^2/delta^2),
%       Default value is 3. Closely related to the speed of convergence.
%
%   lambda: Represents the trade-off between the goodness of data fit 
%       and regularization. Default value is 1. Control convergence results.  
%
%   theta: Define how could be an inlier. If the posterior probability of
%       a sample is an inlier is larger than theta, then it is regarded as an
%       inlier. Default value is 0.75.
%
%   a: Paramerter of the model of outliers. We assume the outliers obey
%       uniform distribution, and the volume of outlier's variation space is a.
%       Default Value is 10.
%
%   MaxIter: Maximum iterition times. Defualt value is 500.
%
%   ecr: The minimum limitation of energy change rate in the iteration
%       process. Default value is 1e-5.
%
%   minP: The posterior probability Matrix P may be singular for matrix
%       inversion. We set the minimum value of P as minP. Default value is
%       1e-5.
%
%   Author: Yifan Xia (xiayifan@whu.edu.cn)


if ~isfield(conf,'MaxIter'), conf.MaxIter = 500; end;
if ~isfield(conf,'gamma'), conf.gamma = 0.9; end;
if ~isfield(conf,'delta'), conf.delta = 1; end; % 
if ~isfield(conf,'lambda'), conf.lambda = 0.01; end; %   
if ~isfield(conf,'theta'), conf.theta = 0.85; end; % 0.85
if ~isfield(conf,'a'), conf.a = 1.0; end;
if ~isfield(conf,'ecr'), conf.ecr = 1e-5; end;
if ~isfield(conf,'minP'), conf.minP = 1e-5; end;
