function Chain = sampleCrpPosterior(x,Chain)
%SAMPLECRPPOSTERIOR Gibbs Sampler for Chinese Restaurant Process
%   CHAIN = SAMPLECRPPOSTERIOR(X,CHAIN) samples the table assignments for
%   the each of the Yule-Simon partitions in the CHAIN structure
%   conditioned on the current state of the Markov chain and the
%   observations in the X vector.
%
%   Copyright 2017 Asher Hensley
%   Revisions:
%   1.2     18-Jul-2017     Hensley     Fixed Gaussian likelihood overflow
%   1.1     17-Jun-2017     Hensley     Improved speed of Gaussian
%                                       likelihood.
%   1.0     13-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.2 $  $Date: 2017/07/18 19:23:54 $
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
N = length(Chain.z);
shape = Chain.shape0;
rate = Chain.rate0;

%Gaussian In-Line Likelihood
L = @(z,s)prod((2*pi*s^2)^(-1/2)*exp(-0.5/s^2*z.^2))+1e-300;

%Remove Empty Tables
[n,tbl] = hist(Chain.z,1:length(Chain.lambda));
Chain = removeEmptyTables(Chain,n,tbl);

%Resample Table Assignments
for kk = 1:length(Chain.z)
    
    %Get Current Counts
    [n,tbl] = hist(Chain.z,1:length(Chain.lambda));
    
    %Remove Current Assignment
    n(Chain.z(kk)) = n(Chain.z(kk))-1;
    
    %Form Weights
    mask = Chain.y==kk;
    w = 0*tbl;
    for ii = 1:tbl(end)
        w(ii) = L(x(mask),1/sqrt(Chain.lambda(ii)));
        w(ii) = n(ii)/(N-1+Chain.gamma0)*w(ii);
    end
    wnew = prod(studentsPdf(x(mask),0,shape/rate,2*shape));
    wnew = Chain.gamma0/(N-1+Chain.gamma0)*wnew;
    w = [w,wnew];
        
    %Clip/Scale as Needed
    wtemp = min(w,exp(709));
    if isinf(sum(wtemp))
        wtemp = wtemp*1e-6;
    end
    if ~isequal(w,wtemp)
        warning('Overflow in sampleCrpPosterior.m')
    end
    w = wtemp;
    
    %Update State Indicator (i.e. Table Number)
    Chain.z(kk) = discreteSelect(w,1);
    
    %Error Chk
    if Chain.z(kk)==0
        error('Bad CRP Table Assignment')
    end
    
    %Update State Parameters (i.e. Dish)
    if Chain.z(kk)==tbl(end)+1
        lnew = sampleLambdaPosterior(x(mask),shape,rate);
        Chain.lambda = [Chain.lambda,lnew];
        n = [n,1];
        tbl = [tbl,tbl(end)+1];
    end

    %Remove Empty Tables
    Chain = removeEmptyTables(Chain,n,tbl);
    
end

%Keep Table Ids in Ascending Order
Chain = sortTables(Chain);

%Sample Dishes
Chain = sampleCrpDishes(x,Chain);


function Chain = removeEmptyTables(Chain,n,tbl)
%REMOVEEMPTYTABLES Remove Empty Tables and Reset Table indices
%   CHAIN = REMOVEEMPTYTABLES(X,N,TBL) removes empty tables from the CHAIN
%   structure using the table counts N and table indices TBL.
%
%   Revisions:
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

Chain.lambda(n==0) = [];
tblMask = tbl(n>0);
for kk = 1:length(tblMask)
    Chain.z(Chain.z==tblMask(kk)) = kk;
end


function Chain = sortTables(Chain)
%SORTTABLES Keep Table IDs in Ascending Order
%   CHAIN = SORTTABLES(CHAIN) Modifies tables and parameters to keep table
%   indices in ascending order.
%
%   Revisions:
%   1.0     22-Jun-2017     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.0 $  $Date: 2017/06/22 19:23:54 $
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

if any(sort(Chain.z)~=Chain.z)
    ptr = 1;
    mapping.old = [];
    mapping.new = [];
    for ii = 1:length(Chain.z)
        if ~any(mapping.old==Chain.z(ii))
            mapping.old(ptr) = Chain.z(ii);
            mapping.new(ptr) = ptr;
            ptr = ptr+1;
        end
    end
    temp = Chain.z;
    for ii = 1:length(mapping.new)
        mask = Chain.z==mapping.old(ii);
        temp(mask) = mapping.new(ii);
    end
    
    %Debug
    if temp(1)~=1
        error('CRP Sorting Error')
    end
    
    Chain.z = temp;
    Chain.lambda = Chain.lambda(mapping.old);
end


function Chain = sampleCrpDishes(x,Chain)
%SAMPLECRPPARAMETERS Sample the CRP Parameters for Each Table
%   CHAIN = SAMPLECRPPARAMETERS(X,CHAIN) samples the parameters (i.e.
%   dishes) for each existing table in the current CRP configuration given
%   the current state of the Markov chain CHAINand the observed data X.
%
%   Revisions:
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

%Setup
shape = Chain.shape0;
rate = Chain.rate0;
tbl = unique(Chain.z);

%Run
for kk = 1:length(tbl)
    mask = Chain.z(Chain.y)==tbl(kk);
    Chain.lambda(kk) = sampleLambdaPosterior(x(mask),shape,rate);
end



