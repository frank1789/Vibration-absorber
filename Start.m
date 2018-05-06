%% Automatic control project
close all;
clear;
clc;
addpath(genpath('/Users/francesco/Documents/MATLAB/Automatic Control/toolboxes'));
yalmip('clear');

% Data definition
% mass [kg]
m1 = 100;
m2 = 15;
% spring modulus [N/m]
k1 = 15e3;
k2 = 2e3;
% damping modulus [Ns/m]
c2 = 30;

%% Component equation of motion
M = [m1 0; 0 m2]; % Defining Mass matrix
Cdamp = [c2 -c2; -c2 c2]; % Damping matrix
K = [k1+k2 -k2; -k2 k2]; % Stiffness matrix

%% Specify the state-space model
sys = ss([-2 -1;1 -2],[1 1;2 -1],[1 0],[0 1]);

%% Convert State-Space Model to Transfer Function
% Convert this model to a transfer function.
tf(sys)

bodeplot(sys,'g--')