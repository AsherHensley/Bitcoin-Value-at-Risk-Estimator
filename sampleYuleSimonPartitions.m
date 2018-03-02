function Chain = sampleYuleSimonPartitions(x,Chain)
%SAMPLEYULESIMONPARTITIONS Gibbs Sampler for Latent Yule Simon Partitions
%   CHAIN = SAMPLEYULESIMONPARTITIONS(X,CHAIN) updates the Yule-Simon 
%   partitions in the CHAIN struct for the Gibbs sampler based on the 
%   observations X.
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

%Sample First Point
Chain = sampleCurrentPoint(x,1,Chain);

%Sample Internal Points
N = length(Chain.y);
for t = 2:N-1
    C = Chain.y(t)~=Chain.y(t-1) || Chain.y(t)~=Chain.y(t+1);
    if C
        Chain = sampleCurrentPoint(x,t,Chain);
    end
end

%Sample Last Point
Chain = sampleCurrentPoint(x,N,Chain);


function Chain = sampleCurrentPoint(x,t,Chain)
%SAMPLECURRENTPOINT Gibbs Sampler for a Given Point in the Sequence
%   CHAIN = SAMPLECURRENTPOINT(X,T,CHAIN) updates the Yule-Simon 
%   partitions in the CHAIN struct for the Gibbs sampler based on the 
%   observation X(T).
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

%Get Transition Weights
[w,type] = getTransitionWeights(x,t,Chain);

%Return If Interior Point
if strcmp(type,'interiorPoint')
    return
end

%Update Markov Chain
Chain = updateMarkovChain(w,type,x,t,Chain);


function [w,type] = getTransitionWeights(x,t,Chain)
%GETTRANSITIONWEIGHTS Markov Chain Transition Probabilities
%   [W,TYPE] = GETTRANSITIONWEIGHTS(X,T,CHAIN) returns the transition
%   weights W with the corresponding TYPE based on the current state of the
%   Markov Chain CHAIN and the observation X(T).
%
%   Revision History:
%   1.1     28-Jul-2017     Hensley     Optimized histogram calculation
%   1.0     13-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.1 $  $Date: 2017/07/28 19:23:54 $
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

%Parameters
alpha = Chain.alpha0;
sigma = 1./sqrt(Chain.lambda);
shape = Chain.shape0;
rate = Chain.rate0;
y = Chain.y;
z = Chain.z;

%In-Line Functions
F = @(xt,idx)normpdf(xt,0,sigma(z(idx))); 
H = @(xt)alpha/(1+alpha)*studentsPdf(xt,0,shape/rate,2*shape);

%Counts
n = diff(find([1,diff(y),1]));

N = length(x);
L = length(n);
j = y(t);

%Init
w = nan;
type = 'interiorPoint';

%Determine Case
if t==1 %First Sample
    
    if y(t+1)==1
        w(1) = (n(1)-1)/(n(1)+alpha)*F(x(1),1);
        w(2) = H(x(1));
        type = 't1_noBoundary';
        
    else
        w(1) = n(2)/(n(2)+alpha+1)*F(x(1),2);
        w(2) = H(x(1));
        type = 't1_rightBoundary';
        
    end
    
elseif t==N %Last Sample
    
    if y(t-1)==L
        w(1) = (n(L)-1)/(n(L)+alpha)*F(x(N),L);
        w(2) = H(x(N));
        type = 'tN_noBoundary';
     
    else
        w(1) = n(L-1)/(n(L-1)+alpha+1)*F(x(N),L-1);
        w(2) = H(x(N));
        type = 'tN_leftBoundary';

    end
    
elseif y(t-1)~=j && y(t+1)==j %Left Boundary
    
    w(1) = (n(j)-1)/(n(j)+alpha)*F(x(t),j);
    w(2) = n(j-1)/(n(j-1)+alpha+1)*F(x(t),j-1);
    w(3) = H(x(t));
    type = 'leftBoundary';
    
elseif y(t-1)==j && y(t+1)~=j %Right Boundary
    
    w(1) = (n(j)-1)/(n(j)+alpha)*F(x(t),j);
    w(2) = n(j+1)/(n(j+1)+alpha+1)*F(x(t),j+1);
    w(3) = H(x(t));
    type = 'rightBoundary';
    
elseif y(t-1)~=j && y(t+1)~=j %Double Boundary
    
    w(1) = n(j-1)/(n(j-1)+alpha+1)*F(x(t),j-1);
    w(2) = n(j+1)/(n(j+1)+alpha+1)*F(x(t),j+1);
    w(3) = H(x(t));
    type = 'doubleBoundary';

end


function Chain = updateMarkovChain(w,type,x,t,Chain)
%UPDATEMARKOVCHAIN Update State of Markov Chain
%   CHAIN = UPDATEMARKOVCHAIN(W,TYPE,X,T,CHAIN) updates the Markov Chain
%   struct CHAIN based on the transition weights W, the sampling case TYPE,
%   and the observation X(T).
%
%   Revision History:
%   1.1     23-Jun-2017     Hensley     Fixed t1_rightBoundary bug
%   1.0     13-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.1 $  $Date: 2017/06/23 19:23:54 $
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
j = Chain.y(t);
u = discreteSelect(w,1);

%Update Chain Struct
switch type
    
    case 't1_noBoundary'
        
        if u==1 %No Change
            Chain.y(1) = 1;
            
        elseif u==2 %Add New Partition
            Chain.y = Chain.y+1;
            Chain.y(1) = 1;
            new = sampleCrpStatic(x(1),Chain);
            Chain.z = [new,Chain.z];
            
        else
            error([type ' error'])
            
        end
      
    case 't1_rightBoundary'
        
        if u==1 %Merge Right
            Chain.y(1) = 2;
            Chain.y = Chain.y-1; 
            Chain.z(1) = [];
            
        elseif u==2 %Add New Partition
            Chain.y(1) = 1;
            new = sampleCrpStatic(x(1),Chain);
            Chain.z(1) = new;
            
        else
            error([type ' error'])
            
        end
        
    case 'leftBoundary'
        
        if u==1 %no change
            Chain.y(t) = j;
            
        elseif u==2 %merge left
            Chain.y(t) = j-1;
            
        elseif u==3 %add new partition
            Chain.y(t+1:end) = Chain.y(t+1:end)+1;
            new = sampleCrpStatic(x(t),Chain);
            Chain.z = [Chain.z(1:j-1),new,Chain.z(j:end)];
            
        else
            error([type ' error'])
            
        end
        
    case 'rightBoundary'
        
        if u==1 %no change
            Chain.y(t) = j;
            
        elseif u==2 %merge right
            Chain.y(t) = j+1;
            
        elseif u==3 %add new partition
            Chain.y(t) = j+1;
            Chain.y(t+1:end) = Chain.y(t+1:end)+1;
            new = sampleCrpStatic(x(t),Chain);
            Chain.z = [Chain.z(1:j),new,Chain.z(j+1:end)];
            
        else
            error([type ' error'])
            
        end
        
    case 'doubleBoundary'
        
        if u==1 %merge left
            Chain.y(t) = j-1;
            Chain.y(t+1:end) = Chain.y(t+1:end)-1;
            Chain.z(j) = [];
            
        elseif u==2 %merge right
            Chain.y(t+1:end) = Chain.y(t+1:end)-1;
            Chain.z(j) = [];
            
        elseif u==3 %add new partition
            new = sampleCrpStatic(x(t),Chain);
            Chain.z(j) = new;
            
        else
            error([type ' error'])
        end
      
    case 'tN_noBoundary'
        
        if u==1 %no change
            Chain.y(t) = j;
        
        elseif u==2 %add new partition
            Chain.y(t) = j+1;
            new = sampleCrpStatic(x(t),Chain);
            Chain.z = [Chain.z,new];
            
        else
            error([type ' error'])
        end

    case 'tN_leftBoundary'
        
        if u==1 %merge left
            Chain.y(t) = j-1;
            Chain.z(j) = [];
            
        elseif u==2 %add new partition
            new = sampleCrpStatic(x(t),Chain);
            Chain.z(j) = new;

        else
            error([type ' error'])
            
        end
        
end


function newTable = sampleCrpStatic(xt,Chain)
%SAMPLECRPSTATIC Draw Samples From Static Chinese Restaurant Process
%   NEWTABLE = SAMPLECRPFORWARD(XT,CHAIN) returns NEWTABLE assignment for
%   based on the observed data XT and the state of the Markov chain CHAIN
%   without the possibility of creating new tables.
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

counts = hist(Chain.z,1:length(Chain.lambda));
prior = counts/sum(counts);
L = normpdf(xt,0,1./sqrt(Chain.lambda));
w = L.*prior;
w = w/sum(w);
newTable = discreteSelect(w,1);



