%TRAIN_SCRIPT_BITCOIN() Bitcoin Data Training Script
%   TRAIN_SCRIPT_SYNTHETIC() runs training phase of Bitcoin data
%   experiment. Results are saved to the file TRAIN_SCRIPT_bitcoin.mat
%   and diagnostic plots are created.
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

%Clean Up
clear
close all
clc

%Read BTC Data
[x,t] = BTC();

%Scale 
scale = 10;
x = scale * x;

%Training Data
winsize = 1500; 
xtrain = x(1:winsize);
ttrain = t(1:winsize);

%Random Number Generator
seed = 3;
rng(seed);

%Gibbs Sampler Setup
Config.a0 = 1;
Config.b0 = 0.1;
Config.alpha0 = 3;              
Config.gamma0 = 20;
Config.shape0 = 2;
Config.rate0 = 1;

%Training
nburn = 2000;          
ngibbs = 8000;
tic
[Chain,History] = gibbsUpdate(xtrain,Config,nburn+ngibbs);
toc

%Save Results to File
save TRAIN_SCRIPT_bitcoin.mat

%Make Diagnostic Plots
trainingDiagnostics(History,ttrain,xtrain/scale,nburn)


