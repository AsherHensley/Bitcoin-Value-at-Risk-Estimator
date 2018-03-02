function L = crpLikelihood(z,gamma0,type)
%CRPLIKELIHOOD Compute Likelihood for Given CRP Seating Arrangement
%   L = CRPLIKELIHOOD(Z,GAMMA0,TYPE) returns the likelihood for the current
%   seating arrangement Z given concentration parameter GAMMA0. TYPE can
%   either be 'log' for log-likelihood or 'linear' for the standard
%   likelihood.
%
%   Revisions:
%   1.0     12-Jan-2018     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.0 $  $Date: 2018/01/12 19:23:54 $
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

%Loop Through Customer Arrivals
L = z*0;
L(1) = 1;
n = 1;
for kk = 2:length(L)
    PDF = [n,gamma0]/(kk-1+gamma0);
    L(kk) = PDF(z(kk));
    if z(kk)==length(PDF)
        n = [n,1];
    else
        n(z(kk)) = n(z(kk))+1;
    end
end

%Compute Output
switch type
    case 'log'
        L = sum(log(L));
    case 'linear'
        L = prod(L);
    otherwise
        error('Unknown Likelihood Type')
end


