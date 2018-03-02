function [VaR,ES] = valueAtRisk(History,nsamp,dn,q,xi,varargin)
%VALUEATRISK Value-at-Risk Predictor
%   [VAR,ES] = VAR(HISTORY,NSAMP,DN,Q,XI) value-at-risk VAR and expected
%   shortfall predictions ES given the sampled states of the Markov chain
%   in the HISTORY structure. Predicitions are made based on NSAMP Gibbs
%   samples with a downsample factor of DN. VaR predicitons are computed
%   based on probability of exceedance Q using an initial guess of XI for
%   the optimization problem.
%
%   [VAR,ES] = VAR(...,type) TYPE can be used as an optional input to set
%   the way VAR and ES predictions are computed. TYPE is string which can
%   be set to 'multi-sample' or 'single-sample'. The multi-sample
%   configuration (default) solves for a single VAR and ES prediction using
%   all of the NSAMP Gibbs samples. The single-sample configuration returns
%   NSAMP VAR and ES predictions (one for each sample).
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 0.1 $  $Date: 2018/02/28 19:23:54 $
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

%Setup
ntotal = size(History.LL,2);
F = @(z,L)normcdf(z,0,L.^(-1/2));
E = @(theta,L) -1./sqrt(2*pi*L).*exp(-theta^2*L/2);

%Run
ptr = 1;
h = cell(1,nsamp);
g = cell(1,nsamp);
s = fliplr(ntotal:-dn:1);
svec = s(end-nsamp+1:end);
for kk = svec
    
    %Latent Variables
    yt = History.y(end,kk);
    nyt = sum(History.y(:,kk)==yt);
    lambda_zyt = History.lambda(end,kk);
    
    %Hyperparameters (Sampled)
    alpha = History.alpha(kk);
    shape = History.shape(kk);
    rate = History.rate(kk);
    GAMMA = History.gamma(kk);
    
    %Get Lambdas
    [~,I] = unique(History.y(:,kk),'first');
    lambdas = History.lambda(I,kk);
    [uzt,J] = unique(History.zt{kk},'first');
    lambdas = lambdas(J);
    
    %Get CRP Table Counts
    nzk = hist(History.zt{kk},uzt);
    
    %VaR Objective Function
    h{ptr} = @(theta) ...
        nyt/(nyt+alpha) * F(theta,lambda_zyt) + alpha/(nyt+alpha) * ( ...
        nzk/(yt+GAMMA) * F(theta,lambdas) + ...
        GAMMA/(yt+GAMMA) * tcdf(sqrt(shape/rate)*theta,2*shape));
    
    %Expected Shortfall
    g{ptr} = @(theta) nan;
    if shape>1/2
        g{ptr} = @(theta) 1/q*( nyt/(nyt+alpha)*E(theta,lambda_zyt) + ...
            alpha/(nyt+alpha) * ( nzk/(yt+GAMMA)*E(theta,lambdas) + ...
            GAMMA/(yt+GAMMA) * (-1) * rate^shape*gamma(shape-0.5) / ...
            (sqrt(2*pi)*gamma(shape) * (theta^2/2+rate)^(shape-0.5)) ) );
    end

    %Update Pointer
    ptr = ptr+1;
    
end

%Clear Temp Variables
clear yt nyt alpha lambda_zyt shape rate lambdas nzk

%Get Type
type = 'multi-sample';
if ~isempty(varargin)
    type = varargin{1};
end

%Configure Output
switch type
    
    case 'multi-sample' 
        
        %VaR
        fmin = [];
        for kk = 1:nsamp
            idx = num2str(kk);
            fmin = [fmin,'+h{',idx,'}(theta)'];
        end
        fmin(1) = [];
        fmin = eval(['@(theta) 1/nsamp * (',fmin,') - q']);
        
        %Solve for VaR
        VaR = fzero(fmin,xi);
        
        %Compute ES
        ES = zeros(1,nsamp);
        for kk = 1:nsamp
            ES(kk) = g{kk}(VaR);
        end
        ES = nanmean(ES);
        
    case 'single-sample'
        
        %Solve for VaR/ES Separately for Each Gibbs Sample
        VaR = zeros(1,nsamp);
        ES = zeros(1,nsamp);
        for kk = 1:nsamp
            fmin = @(theta)h{kk}(theta)-q;
            VaR(kk) = fzero(fmin,xi);
            ES(kk) = g{kk}(VaR(kk));
        end
        
    otherwise
        error('I dont know that one')
        
end




