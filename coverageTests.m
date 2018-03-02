function coverageTests(x,p)
%COVERAGETESTS() Run/Print Unconditional and Conditional Coverage Tests
%   COVERAGETESTS(X,P) Runs and prints results for unconditional and
%   conditional coverage hypothesis tests based on binary hit (VaR
%   exceedance) sequence X and null probability P.
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

%UC Test
ph = mean(x);
n = length(x);
k = sum(x);
L_uc = 2 * (log(ph^k*(1-ph)^(n-k)) - log(p^k*(1-p)^(n-k)));
pval_uc = 1-chi2cdf(L_uc,1);

%Print UC Test Results
disp(['# Violations = ' num2str(sum(x)) '/' num2str(n)])
disp(['p_MLE = ' num2str(mean(ph))])
disp(' ')
disp('UC TEST:')
disp(['Test Statistic = ' num2str(L_uc)])
disp(['p-value = ' num2str(pval_uc)])
disp(' ')

%IND Test
T00 = sum(x(1:end-1)==0 & x(2:end)==0);
T01 = sum(x(1:end-1)==0 & x(2:end)==1);
T10 = sum(x(1:end-1)==1 & x(2:end)==0);
T11 = sum(x(1:end-1)==1 & x(2:end)==1);
pi01 = T01/(T00+T01);
pi11 = T11/(T10+T11);
pi0 = (T01+T11)/n;
LA = (1-pi01)^T00 * pi01^T01 * (1-pi11)^T10 * pi11^T11;
L0 = (1-pi0)^(T00+T10) * pi0^(T01+T11);
L_ind = 2 * (log(LA) - log(L0));
pval_ind = 1-chi2cdf(L_ind,1);

%Print IND Test Results
disp('IND TEST:')
disp(['Test Statistic = ' num2str(L_ind)])
disp(['p-value = ' num2str(pval_ind)])
disp(' ')

%CC Test
L_cc = L_uc + L_ind;
pval_cc = 1-chi2cdf(L_cc,2);

%Print CC Test Results
disp('CC TEST:')
disp(['Test Statistic = ' num2str(L_cc)])
disp(['p-value = ' num2str(pval_cc)])
disp(' ')

