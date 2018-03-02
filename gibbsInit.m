function Chain = gibbsInit(x,Config)
%GIBBSINIT Initialize State of Markov Chain for Gibbs Sampler
%   CHAIN = GIBBSINIT(X,CONFIG) initializes the state of the Markov Chain
%   for the Gibbs sampling algorithm and returns the result in the CHAIN
%   structure based on the observations vector X an configuration struct
%   CONFIG.
%
%   Revision History:
%   1.0     13-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.0 $  $Date: 2017/06/13 19:23:54 $
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

%Set Default Parameters
Chain.a0 = Config.a0;
Chain.b0 = Config.b0;
Chain.alpha0 = Config.alpha0;
Chain.gamma0 = Config.gamma0;
Chain.shape0 = Config.shape0;
Chain.rate0 = Config.rate0;
Chain.y = [1,x(2:end)*0];
Chain.z = 1;
shape = Chain.shape0;
rate = Chain.rate0;
Chain.lambda = sampleLambdaPosterior(x(1),shape,rate);

%Init Chain
k = 1;
wtest = [];
for t = 2:length(x)
    
    %Transition Probabilities
    w = forwardTransitionWeights(x(t),Chain,t);
    wtest = [wtest; w];
    
    %Sample Transition
    u = rand;
    stateTransition = ~(u<w(1));
    if stateTransition
        k = k+1;
        Chain.z(k) = sampleCrpForward(x(t),Chain);
    end
    Chain.y(t) = k;
    
    %Update Lambda
    mask = Chain.z(Chain.y(1:t))==Chain.z(k);
    newLambda = sampleLambdaPosterior(x(mask),shape,rate);
    if Chain.z(k)>length(Chain.lambda)
        Chain.lambda = [Chain.lambda,newLambda];  %New Table
    else
        Chain.lambda(Chain.z(k)) = newLambda;
    end
    
end


function w = forwardTransitionWeights(x,Chain,ptr)
%FORWARDTRANSITIONWEIGHTS Yule-Simon Forward Process Transition Weights
%   W = FORWARDTRANSITIONWEIGHTS(X,CHAIN,PTR) returns the Yule-Simon
%   forward transition weights based on the scalar observation X, state of
%   the Markov Chain CHAIN, and the current time index pointer PTR.
%
%   Revision History:
%   1.0     13-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.0 $  $Date: 2017/06/13 19:23:54 $
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

%Unpack
ylast = Chain.y(ptr-1);
n = sum(Chain.y(1:ptr-1)==ylast);
alpha = Chain.alpha0;
lambda = Chain.lambda(Chain.z(ylast));
shape = Chain.shape0;
rate = Chain.rate0;

%Compute
prior = [n./(n+alpha),alpha./(n+alpha)];
L0 = normpdf(x,0,1/sqrt(lambda));
L1 = studentsPdf(x,0,shape/rate,2*shape);
w = [L0,L1].*prior;
w = w./sum(w);


function newTable = sampleCrpForward(xt,Chain)
%SAMPLECRPFORWARD Draw Samples From Chinese Restaurant Process
%   NEWTABLE = SAMPLECRPFORWARD(XT,CHAIN) returns NEWTABLE assignment
%   based on the observed data XT and the state of the Markov chain CHAIN
%
%   Revision History:
%   1.0     13-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.0 $  $Date: 2017/06/13 19:23:54 $
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

counts = hist(Chain.z,unique(Chain.z));
prior = [counts,Chain.gamma0]/(sum(counts)+Chain.gamma0);
L0 = normpdf(xt,0,1./sqrt(Chain.lambda));
L1 = studentsPdf(xt,0,Chain.shape0/Chain.rate0,2*Chain.shape0);
w = [L0,L1].*prior;
w = w./sum(w);
newTable = discreteSelect(w,1);



