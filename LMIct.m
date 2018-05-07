%% LMI control toolbox
clear
clc
load statespace.mat

% Initialize the LMI system
setlmis([]);
% Specify the variables of the LMI
P     = lmivar(1,[length(A) 1]);
gamma = lmivar(1,[1 0]);
% Describe the LMI constraints
% 1st LMI (P>0)
Ppos = newlmi;
lmiterm([-Ppos 1 1 P],1,1);

% 2nd LMI (L2 gain)
L2lmi = newlmi;
lmiterm([L2lmi 1 1 P],1,A,'s');
lmiterm([L2lmi 1 2 P],1,B);
lmiterm([L2lmi 2 2 gamma],-1/2,eye(size(B,2)),'s');
lmiterm([L2lmi 3 1 0],C);
lmiterm([L2lmi 3 2 0],D);
lmiterm([L2lmi 3 3 gamma],-1/2,eye(size(C,1)),'s');
% Assign a name to the LMI system
mylmisys=getlmis;
% Solve the LMI
% Choose the function to be minimized (gamma)
n = decnbr(mylmisys);
cost = zeros(n,1);
for j=1:n
    cost(j)=defcx(mylmisys,j,gamma);
end
% Run the LMI solver
[copt,xopt]=mincx(mylmisys,cost);
% Decode the solution (if any)
if not(isempty(copt))
    Psol     = dec2mat(mylmisys,xopt,P);
    gammasol = dec2mat(mylmisys,xopt,gamma);
    % compare to alternative Hinf norm computation
    sys = pck(A,B,C,D);
    out = hinfnorm(sys);
    disp([out(2) gammasol])
end