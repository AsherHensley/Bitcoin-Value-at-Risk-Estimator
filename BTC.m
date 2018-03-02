function [r,t] = BTC()
%BTC() Read Bitcoin Data BTC.xlsx
%   [R,T] = BTC() returns Bitcoin log return vector R with fraction year 
%   suport vector T.
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

%Import
[num,txt] = xlsread('BTC.xlsx');

%Closing Prices
p = num(:,strcmp(txt,'Close Price'));
t = num(:,strcmp(txt,'Datenum'));

%Down-Sampling
dn = 1;
p = p(1:dn:end);
t = t(1:dn:end);

%Log Returns
r = diff(log(p));
t = t(2:end);

%Output
r = r(end-1999:end);
t = datevec(t(end-1999:end)) * [1,1/12,1/365,0,0,0]';


