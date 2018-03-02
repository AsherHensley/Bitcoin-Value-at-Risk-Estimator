function trainingDiagnostics(History,ttrain,xtrain,nburn)
%TRAININGDIAGNOSTICS() Make Training Diagnostics Plots
%   TRAININGDIAGNOSTICS(HISTORY,TTRAIN,XTRAIN,NBURN) makes training 
%   diagnostic plots for Gibbs sampler inference and convergence using
%   HISTORY structure, training data vectors TTRAIN and XTRAIN, and the
%   number of burned Gibbs sample NBURN.
%
%   Revisions:
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

%Plot Hyperparameter Samples
figure,
subplot(221)
plot(History.rate)
set(gca,'fontsize',12)
xlabel('MCMC Iteration')
ylabel('Rate')
grid on
axis tight

subplot(222)
plot(History.alpha)
set(gca,'fontsize',12)
xlabel('MCMC Iteration')
ylabel('Alpha')
grid on
axis tight

subplot(223)
plot(History.gamma)
set(gca,'fontsize',12)
xlabel('MCMC Iteration')
ylabel('Gamma')
grid on
axis tight

subplot(224)
plot(History.shape)
set(gca,'fontsize',12)
xlabel('MCMC Iteration')
ylabel('Shape')
grid on
axis tight

%Plot Hyperparameter Histograms
figure,
subplot(221)
hist(History.rate(nburn+1:end),30)
set(gca,'fontsize',11)
ylabel('Histogram')
xlabel('Rate')
grid on

subplot(222)
hist(History.alpha(nburn+1:end),30)
set(gca,'fontsize',11)
ylabel('Histogram')
xlabel('Alpha')
grid on

subplot(223)
hist(History.gamma(nburn+1:end),0.25:0.25:20)
set(gca,'fontsize',11)
ylabel('Histogram')
xlabel('Gamma')
grid on
xlim([1,25])

subplot(224)
hist(History.shape(nburn+1:end),0.05:0.05:2)
set(gca,'fontsize',11)
ylabel('Histogram')
xlabel('Shape')
grid on
axis tight

%Plot Training Data, Log Posterior, and Number of Regimes
figure,
subplot(2,2,1:2)
plot(ttrain,xtrain,'k')
grid on
axis tight
ax = axis;
xlim([ax(1:2)+[-0.1,0.1]])
set(gca,'fontsize',11)
title('Training Data')
xlabel('Sample')
ylabel('Daily Log Returns')

subplot(223)
plot(sum(History.LL),'k')
set(gca,'fontsize',11)
xlabel('MCMC Iteration')
ylabel('Log Posterior')
axis tight
ax = axis;
xlim([-500,ax(2)+500])
grid on

subplot(224)
plot(History.y(end,:),'k')
set(gca,'fontsize',11)
xlabel('MCMC Iteration')
ylabel('# Volatility Clusters')
axis tight
ax = axis;
xlim([-500,ax(2)+500])
grid on

%Plot Regime Change Probabilities
figure,
subplot(211),
plot(ttrain,xtrain,'k'),axis tight
set(gca,'fontsize',11)
xlabel('Sample')
ylabel('Log Returns')
grid on
ha(1) = gca;

subplot(212),
plot(ttrain,[0;mean(diff(History.y(:,nburn+1:end)),2)],'k')
set(gca,'fontsize',11)
xlabel('Sample')
ylabel('Pr(Regime Change)')
grid on
ha(2) = gca; 
axis tight
linkaxes(ha,'x')



