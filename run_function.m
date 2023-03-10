%--------------------------------------------------------------------------
% Project: System identification of discretized hybrid systems
%--------------------------------------------------------------------------

function [t1, j1, errornorm] = run_function(example_number)

%% Setup simulations
global s gamma_c1 gamma_d1 dimpsi dimt dimth A B H J u C D error psi theta kj1

%pick example 
[s, gamma_c1, gamma_d1, psi0, A, B, H, J, u, C, D, TSPAN, JSPAN] ...
    = pickexample(example_number);

dimpsi = size(A,1); dimu = size(u,1); dimt = dimpsi + dimu;

%% Simulation 1: generate y_s and psi_s
% initial conditions
t0 = 0; k0 = 0; j0 = 0;
x0 = [psi0, t0, k0, j0];

% priority: rule = 1 (jumps), 2 (flows)
rule = 1;
options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[t1,j1,x1] = HyEQsolver( @f_est1,@g_est1,@C_est1,@D_est1,...
    x0,TSPAN,JSPAN,rule,options,'ode23t');

kj1 = x1(:,end-1:end);

%% Simulation 2: compute theta trajectory
for i = 1:size(x1,1)-1
    if kj1(i+1,2) == kj1(i,2)
        psi(i,:) = [x1(i,1:end-3)'; u(i); zeros(dimt,1)];
    elseif kj1(i+1,2) ~= kj1(i,2)
        psi(i,:) = [zeros(dimt,1); x1(i,1:end-3)'; u(i)];
    end
end 

theta = [A B H J]; dimth = size(theta);
hattheta0 = theta;
error = hattheta0-theta;
for i = 1:dimth(1)
    hatthetarows(1+(i-1)*dimth(2):dimth(2)+(i-1)*dimth(2)) = hattheta0(i,:);
end

% initial conditions
t0 = 0; k0 = 0; j0 = 0;
x0 = [hatthetarows, t0+1, k0, j0];

% priority: rule = 1 (jumps), 2 (flows)
rule = 1;
options = odeset('RelTol',1e-6,'MaxStep',.1);

% simulate
[t2,j2,x2] = HyEQsolver( @f_est2,@g_est2,@C_est2,@D_est2,...
    x0,TSPAN,[0 size(kj1,1)],rule,options,'ode23t');

thetarows = x2(end,1:end-3);
for i = 1:dimth(1)
    hattheta(i,:) = thetarows(1+(i-1)*dimth(2):dimth(2)+(i-1)*dimth(2));
end

% calculate error norm
for i = 1:size(t2,1)
    errornorm(i) = norm(x2(i,1:end-3));
end

end