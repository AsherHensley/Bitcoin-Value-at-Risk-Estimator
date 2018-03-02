function LL = logLikelihood(x,Chain)
%LOGLIKELIHOOD Compute Log-Likelihood for Given State of Markov Chain
%   LL = LOGLIKELIHOOD(X,CHAIN) returns the log-likelihood LL of the given
%   state of the Markov chain defined by the CHAIN stucture conditioned on
%   the observations X.
%
%   Revisions:
%   1.1     12-Jan-2018     Hensley     Replaced CRP likelihood calculation
%                                       with call to crpLikelihood.m
%   1.0     22-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.1 $  $Date: 2018/01/12 19:23:54 $
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

%p(x|y,z,lambda)
sigma = Chain.lambda(Chain.z(Chain.y)).^(-1/2);
L1 = sum(log(normpdf(x,0,sigma)));

%p(y|alpha)
n = hist(Chain.y,unique(Chain.y));
L2 = sum(log(Chain.alpha0*beta(n,Chain.alpha0+1)));

%p(z|gamma)
L3 = crpLikelihood(Chain.z,Chain.gamma0,'log');

%p(lambda|a,b)
F = @(z,a,b)b^a/gamma(a)*z.^(a-1).*exp(-b*z);
L4 = sum(log(F(Chain.lambda,Chain.shape0,Chain.rate0)));

%Total
LL = [L1; L2; L3; L4];

%Debug
if any(isnan(LL)) || any(isinf(LL))
    error('Calculation problem in loglLikelihood.m')
end

