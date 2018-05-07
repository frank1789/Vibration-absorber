%% Automatic control project
close all
clear
clc
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
C = [1 0 0 0];                          % single output
D = 0;
sys = ss(A,B,C,D);

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
opts.solver = 'sdpt3';  % chose the solver
opts.verbose = 1;       % show iterations 0 = don't show, 1 = show

% set the variable to optimize
P = sdpvar(4, 4, 'symmetric');
obj = trace(P);                             % set the objective function
constr = [ A'* P + P * A < 0; P > eye(length(A))];  % set the constraints

% solve the problem
yalmipdiagnostics = optimize(constr,obj,opts);

% display the results during iteration
if yalmipdiagnostics.problem == 0
    disp('Solver thinks it is feasible')
elseif yalmipdiagnostics.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end

P = double(P); % extract the result
disp(P);

% check that is a proper Lyapunov function
[V,E] = eig(A' * P + P * A);
disp(V); % V of right eigenvectors
disp(E); % E of eigenvalues

% check the constraint
check(constr);
save statespace.mat A B C D

%% 3) estimate numerically L2 gain
clear
load statespace.mat;
yalmip('clear')
% decision variable
P = sdpvar(length(A));
gamma = sdpvar(1);      % scalar

n = length(A);
[ny, nu] = size(D);

% define inequality constraint
M = [P * A          P * B           zeros(n,ny);
    zeros(nu,n) -gamma/2*eye(nu)    zeros(nu,ny);
    C                 D             -gamma/2*eye(ny)'];

constr = [M + M' < 0, P > 0];

opts = sdpsettings;
opts.solver = 'sdpt3';

diagnostic = optimize(constr,gamma,opts);
if diagnostic.problem == 0
    disp('Solver thinks it is feasible')
elseif diagnostic.problem == 1
    disp('Solver thinks it is infeasible')
else
    disp('Something else happened')
end

% evaluate solution variables (if any)
Psol=double(P);
gammasol=double(gamma);
disp(Psol);
disp(gammasol);

% compare to alternative Hinf norm computation
sys = pck(A,B,C,D);
out = hinfnorm(sys);
disp([out(2) gammasol])