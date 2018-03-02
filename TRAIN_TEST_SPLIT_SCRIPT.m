%TRAIN_TEST_SPLIT_SCRIPT() Plot Train/Test Split for SYN & BTC Experiments
%   TRAIN_TEST_SPLIT_SCRIPT() plots train/test splits for SYN and BTC
%   experiments by reading TRAIN_SCRIPT_synthetic.mat and 
%   TRAIN_SCRIPT_bitcoin.mat
%
%   Revision History:
%   1.0     28-Feb-2018     Hensley     Initial release
%
%   MIT License
% 
%   Copyright (c) 2018 Asher A. Hensley Research
%   $Revision: 1.0 $  $Date: 2018/02/28 19:23:54 $
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

%Clean Up Workspace
close all
clear
clc

%Load
d1 = load('TRAIN_SCRIPT_synthetic.mat');
d2 = load('TRAIN_SCRIPT_bitcoin.mat');

%Plot
figure
subplot(211)
plot(d1.x/max(d1.x)*100,'color',[1,1,1]*0.5)
set(gca,'yticklabel',{},'fontsize',12)
xlabel('Sample'),ylabel('Log-Return'),title('Synthetic Data')
axis tight
ax = axis;
hold on
plot(1500*[1,1],ax(3:4),'k','linewidth',2)
grid on
text(1300,80,'\bf TRAIN')
text(1530,80,'\bf TEST')

subplot(212)
plot(d2.t,d2.x/max(d2.x)*100,'color',[1,1,1]*0.5)
set(gca,'yticklabel',{},'fontsize',12,'xtick',[2013:2018])
xlabel('Year'),ylabel('Log-Return'),title('Bitcoin/USD Exchange Rate')
axis tight
ax = axis;
hold on
plot(d2.t(1500)*[1,1],ax(3:4),'k','linewidth',2)
grid on
text(d2.t(1300),80,'\bf TRAIN')
text(d2.t(1530),80,'\bf TEST')



