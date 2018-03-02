function [y,truth] = syntheticData(N)
%SYNTHETICDATA(N) Create Synthetic Log-Returns with Generative Model
%   [Y,TRUTH] = SYNTHETICDATA(N) samples log return sequence Y from
%   generative model. All latent variables and parameters are returned in
%   the TRUTH structure.
%
%   Revision History:
%   1.0     15-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.0 $  $Date: 2017/06/15 19:23:54 $
%   
%   Permission is hereby granted, free of charge, to any person obtaining a 
%   copy of this software and associated documentation files (the 
%   "Software"), to deal in the Software without restriction, including 
%   without limitation the rights to use, copy, modify, merge, publish, 
%   distribute, sublicense, and/or sell copies of the Software, and to 
%   permit persons to whom the Software is furnished to do so, subject to 
%   the following conditions:
% 
%   The above copyright notice and this permission notice shall be included 
%   in all copies or substantial portions of the Software.
% 
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
%   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
%   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
%   CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
%   TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
%   SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%Set Hyperparameters
truth.alpha = 1.1;
truth.gamma = 2;
truth.shape = 1;
truth.rate = 0.2;
truth.seed = 4;

%Set Random Seed
rng(truth.seed);

%Sample Yule-Simon Process
truth.x = ones(1,N);
nk = 1;
for kk = 2:N
    pk = nk/(nk+truth.alpha);
    u = rand;
    st = u>pk;
    if st==0
        truth.x(kk) = truth.x(kk-1);
        nk = nk+1;
    else
        truth.x(kk) = truth.x(kk-1)+1;
        nk = 1;
    end
end

%Sample Chinese Restaurant Process
K = truth.x(end);
truth.z = ones(1,K);
for ii = 2:K
    cur = truth.z(1:ii-1);
    n = hist(cur,1:max(cur));
    p = [n,truth.gamma];
    truth.z(ii) = discreteSelect(p,1);
end

%Sample Volatility
truth.lambda = gamrnd(truth.shape,1/truth.rate,1,max(truth.z));

%Sample Measurements
sigmas = 1./sqrt(truth.lambda(truth.z(truth.x)));
y = sigmas.*randn(1,N);

