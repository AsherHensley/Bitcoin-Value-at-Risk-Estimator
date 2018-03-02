function [Chain,History] = gibbsUpdate(x,S,n)
%GIBBSUPDATE Run Gibbs Sampler for Switching Variance Process
%   [CHAIN,HISTORY] = GIBBSUPDATE(X,N) runs the Gibbs sampling
%   algorithm for the switching variance process for N iterations based on 
%   the observations in the vector X. The variable S can either be the 
%   CONFIG or CHAIN structure. The current state of the Markov 
%   chain is returned in the CHAIN structure along with the HISTORY 
%   structure documenting the trajectory of the Gibbs sampling 
%   algorithm during the N iterations for the latent variables.
%
%   NOTES: The mathematical representation of the generative model uses x
%   to denote latent regime indicators and y to denote observations. This 
%   program and it's supporting routines use the reverse representation.
%
%   Revision History:
%   1.5     25-Jan-2018     Hensley     Fixed zero-liklelihood problem with
%                                       CRP gamma sampler.
%   1.4     24-Jan-2018     Hensley     Force observations to row vector
%   1.3     14-Jan-2018     Hensley     Added separate gamma/shape samplers
%   1.2     11-Jan-2018     Hensley     Removed Metropolis-Hastings
%                                       sampler. Integrated gibbsInit.m.
%   1.1     19-Jul-2017     Hensley     Added Metropolis-Hastings rejection
%                                       count to History struct.
%   1.0     15-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.5 $  $Date: 2018/02/28 19:23:54 $
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

%Force to Row
x = x(:)';

%Varargin
if isfield(S,'z')
    Chain = S;
else
    Chain = gibbsInit(x,S);
end
    
%Setup History Struct
History.y = zeros(length(x),n+1);
History.z = zeros(length(x),n+1);
History.zt = cell(1,n+1);
History.lambda = zeros(length(x),n+1);
History.LL = zeros(4,n+1);
History.shape = zeros(1,n+1);
History.rate = zeros(1,n+1);
History.alpha = zeros(1,n+1);
History.gamma = zeros(1,n+1);

%Update with Current State
History.y(:,1) = Chain.y';
History.z(:,1) = Chain.z(Chain.y');
History.zt{1} = Chain.z;
History.lambda(:,1) = Chain.lambda(Chain.z(Chain.y'));
History.LL(:,1) = logLikelihood(x,Chain);
History.shape(1) = Chain.shape0;
History.rate(1) = Chain.rate0;
History.alpha(1) = Chain.alpha0;
History.gamma(1) = Chain.gamma0;

%Run
hw = waitbar(0,'Running MCMC Update 0 of 0');
ns = num2str(n);
for it = 2:n+1
    
    %Update Yule-Simon Partitions
    Chain = sampleYuleSimonPartitions(x,Chain);

    %Update Partition Parameter Assignments
    Chain = sampleCrpPosterior(x,Chain);
    
    %Sample Alpha Posterior
    Chain = sampleAlphaPosterior(Chain);
    
    %Sample Rate Posterior
    Chain = sampleRatePosterior(Chain);
    
    %Sample Gamma Posterior
    gsweep = 0.25:0.25:20;
    p = zeros(1,length(gsweep));
    for kk = 1:length(gsweep)
        p(kk) = crpLikelihood(Chain.z,gsweep(kk),'log');
    end
    p0 = (length(p):-1:1)/sum(length(p):-1:1);
    p = p + log(p0);
    p = exp(p - min(p) - 745);
    if all(p==0), p = p0; disp('WARNING: CRP Likelihoods are Zero'),end
    Chain.gamma0 = gsweep(discreteSelect(p,1));
    
    %Sample Shape Posterior
    sup = 0.05:0.05:2;
    p = zeros(1,length(sup));
    for kk = 1:length(sup)
        p(kk) = prod(gampdf(Chain.lambda,sup(kk),1/Chain.rate0));
    end
    p0 = (length(p):-1:1)/sum(length(p):-1:1);
    p = p.*p0;
    if all(p==0), p = p0; disp('WARNING: Shape Likelihoods are Zero'),end
    Chain.shape0 = sup(discreteSelect(p,1));

    %Update History
    History.y(:,it) = Chain.y';
    History.z(:,it) = Chain.z(Chain.y');
    History.zt{it} = Chain.z;
    History.lambda(:,it) = Chain.lambda(Chain.z(Chain.y'));
    History.LL(:,it) = logLikelihood(x,Chain);
    History.shape(it) = Chain.shape0;
    History.rate(it) = Chain.rate0;
    History.alpha(it) = Chain.alpha0;
    History.gamma(it) = Chain.gamma0;
    
    %Update Waitbar
    waitbar(it/n,hw,['Running MCMC Update ',num2str(it),' of ' ns]);
    
end
delete(hw)



