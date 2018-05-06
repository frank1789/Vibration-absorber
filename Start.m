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
% Considering the system as time invariant, these equations can also be
% represented in the standard state - space form as,
% x_dot(t) = Ax(t) +Bu(t)
% y (t) = Cx(t)+Du(t)
% where:
% A = System matrix
% B = Input matrix
% C = Output matrix
% D = Disturbance matrix
% x = State Space Vector
% u = Control input Vector

u = [1 0].';

A = [zeros(2) eye(2); -M\K -M\Cdamp];   % left divide for the inverse
B = [zeros(2,1); M\u];                  % single input
C = [1 0 0 0];                          % multiple output

sys = ss(A,B,C,[]);
%% Convert State-Space Model to Transfer Function
tf(sys)

figure(101)
subplot(2,1,1)
step(sys)
subplot(2,1,2)
impulse(sys)

figure(102)
opts = bodeoptions;
opts.Grid = 'on';
opts.MagUnits = 'abs';
bodeplot(tf(sys),opts,'g')

%% Lyapunot stability as an LMI problem
% clear the internal memory of YALMIP
yalmip('clear')

% set the parameters of the LMI solver
opts = sdpsettings;     % set the solver settings
opts.solver = 'Lmilab'; % chose the solver
opts.verbose = 1;       % show iterations 0 = don't show, 1 = show

% set the variable to optimize
P = sdpvar(4, 4, 'symmetric');
obj = trace(P);                             % set the objective function
constr = [ A'* P + P * A < 0; P > eye(length(A))];  % set the constraints

% solve the problem
yalmipdiagnostics = solvesdp(constr,obj,opts);

% display the results during iteration
msg = ['feasibility = ', num2str(yalmipdiagnostics.problem)];
disp(msg);
P = double(P); % extract the result
disp(P);

% check that is a proper Lyapunov function
[V,D,W] = eig(A' * P + P * A);
disp(V); % V of right eigenvectors
disp(D); % D of eigenvalues
