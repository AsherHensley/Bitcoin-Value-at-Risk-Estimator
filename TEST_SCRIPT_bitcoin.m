%TEST_SCRIPT_BITCOIN() Bitcoin Data Test Script
%   TRAIN_SCRIPT_SYNTHETIC() runs prediction phase of Bitcoin data
%   experiment. Results are saved to the file TEST_SCRIPT_bitcoin.mat
%   and diagnostic plots/results are created.
%
%   NOTES: This script may take several hours to run. After each
%   prediction, the current state of the algorithm is saved to temp.mat. By
%   doing this, this script can be stopped at anytime and run later
%   starting from the stopping point. To run this script using the state
%   saved to the temp.mat file, set: from_file=true. 
%
%   Revision History:
%   1.0     28-Feb-2018     Hensley     Initial release
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

%Clean Up Workspace
clear
close all
clc

%Load Training Data
load TRAIN_SCRIPT_bitcoin.mat

%Run Flag
from_file = false;

%Setup
dn = 1;
q = [0.01,0.05];
ngibbs = 100;
N = length(x(winsize+1:end));

%Setup
if ~from_file

    %Init Random Number Generator
    seed = 15;
    rng(seed);
    
    %Init Variables,etc.
    start = 1;
    VaR = zeros(length(q),N);
    ES = zeros(length(q),N);
    xin = xtrain;
    VaR_guess = 0;
    
else
    
    %Load
    load temp.mat
    
    %Init Random Number Generator
    rng(rng_state);
    
    %Init Start Point
    start = tau;
    
end

%Run
for tau = start:N

    %Save State of Random Number Generator
    rng_state = rng;
    
    %Save Results to temp.mat
    save temp.mat
    
    %Prediction
    tic
    for kk = 1:length(q)
        [VaR(kk,tau),ES(kk,tau)] = valueAtRisk(History,ngibbs,dn,q(kk),VaR_guess);
        VaR_guess = VaR(kk,1);
    end
    
    %Update Markov Chain
    xin = x(1:winsize+tau);
    Chain.y = [Chain.y,Chain.y(end)];
    [Chain,History] = gibbsUpdate(xin,Chain,ngibbs);
    
    %Cool Down 
    pause(10)
    
    %Print Out
    dt = toc;
    disp([num2str(tau/N*100),'% Complete - Time Elapsed: ',num2str(dt),' seconds'])

end

%Save
save TEST_SCRIPT_bitcoin.mat

%Plot VaR Results
c1 = [1,1,1]*0.8;
c5 = [1,1,1]*0.6;
figure,
area(t(winsize+1:end),VaR(1,:)','facecolor',c1,'edgecolor',c1)
hold on
area(t(winsize+1:end),VaR(2,:)','facecolor',c5,'edgecolor',c5)
stem(t(winsize+1:end),x(winsize+1:end),'k','marker','none')
set(gca,'fontsize',12,'yticklabel',{})
xlabel('Sample')
ylabel('Observations')
axis tight
grid on
legend('1% Quantile','5% Quantile','Data','location','northwest')

%Print VaR Exceedance/Hypothesis Test Results
xtest = x(winsize+1:winsize+tau);
ttst = t(winsize+1:winsize+tau);
e1 = xtest(1:tau)'<VaR(1,1:tau);
e2 = xtest(1:tau)'<VaR(2,1:tau);
disp(' ')
disp('--------------')
disp('1% Test:')
disp('--------------')
coverageTests(e1,q(1));
disp('--------------')
disp('5% Test:')
disp('--------------')
coverageTests(e2,q(2));

