function p = studentsPdf(x,mu,lambda,nu)
%STUDENTSPDF Student's-t Likelihood
%   P = STRUDENTSPDF(X,MU,LAMBDA,NU) evaluates the Student's-t distribution
%   with mean MU, precision LAMBDA, and degrees of freedom NU at the
%   support values in the vector X.
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

%Student's-t PDF
p = gamma(nu/2+1/2)/gamma(nu/2)*sqrt(lambda/pi/nu)...
    *(1+1/nu*lambda*(x-mu).^2).^(-nu/2-1/2);


