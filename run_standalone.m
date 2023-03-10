% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: run_standalone.m
%--------------------------------------------------------------------------
% Project: Discretized hybrid gradient descent algorithm (DiscHyGradAlg)
% Description: Standalone script to simulate DiscHyGradAlg for chosen
% example, returns figures for sampled signal \psi_s and trajectory of
% error norm \tilde\theta_s
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision:

%% Setup simulations

clc; clear all; close all;

global s gamma_c1 gamma_d1 dimpsi dimt dimth A B H J u C D error psi theta kj1

%pick example 
[s, gamma_c1, gamma_d1, psi0, A, B, H, J, u, C, D, TSPAN, JSPAN] ...
    = pickexample(1);

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

%% Plots
% plot solution
modificatorF2{1} = 'b';
modificatorF2{2} = 'LineWidth';
modificatorF2{3} = 1.5;

modificatorJ2{1} = 'r:';
modificatorJ2{2} = 'LineWidth';
modificatorJ2{3} = 1.5;
modificatorJ2{4} = 'Marker';
modificatorJ2{5} = '.';
modificatorJ2{6} = 'MarkerEdgeColor';
modificatorJ2{7} = 'r';
modificatorJ2{8} = 'MarkerFaceColor';
modificatorJ2{9} = 'r';
modificatorJ2{10} = 'MarkerSize';
modificatorJ2{11} = 1.5;

figure(1) % Plot for psi_s
tlt = tiledlayout(2, 1);
varm = {'$x_1$', '$x_2$'};
for i = 1:dimpsi
    nexttile(i)
    plotHarc(t1,j1,x1(:,i),[],modificatorF2,modificatorJ2);
    grid on
    ylabel(varm{i},'interpreter','latex','FontSize',18);
    if i == 1
        title('Trajectory of $\psi_s$', 'interpreter','latex');
    end
end
xlabel('t [s]');
pos = get(gcf, 'Position');
tlt.Padding = "none";
tlt.TileSpacing = "compact";
set(gcf, 'Position',  [pos(1), pos(2), 1*pos(3), 1*pos(4)])
set(gca, 'LooseInset', get(gca,'TightInset'))

for i = 1:size(t2,1)
    errornorm(i) = norm(x2(i,1:end-3));
end

modificatorF{1} = 'b';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;

modificatorJ{1} = 'b:';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 2;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '.';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'b';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'b';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 2;

figure(2) % Plot for evolution of \hat\theta_s
tlt = tiledlayout(2, 1);
nexttile(1)
plotHarc(t1(1:end-1),j1(1:end-1),errornorm,[],modificatorF,modificatorJ);
title('Norm of error $\tilde\theta_s$ vs Time','interpreter','latex');
xlabel('$t$ [s]','interpreter','latex','FontSize',18);
ylabel('$|\tilde\theta_s|$','interpreter','latex','FontSize',18);
grid on; box on;
pos = get(gcf, 'Position');
tlt.Padding = "none";
tlt.TileSpacing = "compact";
set(gcf, 'Position',  [pos(1), pos(2), 1*pos(3), 1*pos(4)])
set(gca, 'LooseInset', get(gca,'TightInset'))

%% Calculations for Gamma, mu and other parameters to check for stability

x_cont = []; t_cont = []; j_cont = []; ic = 0;
x_disc = []; t_disc = []; j_disc = []; id = 0;
for i = 1:size(x1,1)-1
    if x1(i,end) == x1(i+1,end)
        ic = ic + 1;
        x_cont(:,ic) = x1(i,1:dimpsi);
        t_cont(:,ic) = x1(i,end-2);
        j_cont(:,ic) = x1(i,end);
    elseif x1(i,end) ~= x1(i+1,end)
        id = id + 1;
        x_disc(:,id) = x1(i,1:dimpsi);
        t_disc(:,id) = x1(i,end-2);
        j_disc(:,id) = x1(i,end);
        j_d_ind(:,id) = i;
    end
end
if id ~= 0
    [~,j_diff_ind] = max(diff(j_d_ind));
    vmin = j_d_ind(j_diff_ind)+1; vmax = j_d_ind(j_diff_ind+1);
else
    vmin = 1; vmax = size(x1,1)-1;
end

Gamma = (x1(vmax,end-1)-x1(vmin,end-1))/s + x1(vmax+1,end)-x1(vmin+1,end);
mu = min(eig(s*gamma_c1*x1(vmin:vmax,1:dimt)'*x1(vmin:vmax,1:dimt)+(gamma_d1*x1([vmin+1,vmax+1],1:dimt)'*x1([vmin+1,vmax+1],1:dimt))/(1+gamma_d1*x1([vmin+1,vmax+1],1:dimt)'*x1([vmin+1,vmax+1],1:dimt))));
sigma = mu/(1+sqrt(2*(Gamma+2)^3))^2;
kappa = sqrt(1/(1-sigma));
lambda = 1/(2*(Gamma+1))*log(1/(1-sigma));

eigc_tot = eig(x_cont*x_cont'); eigd_tot = eig(x_disc*x_disc');
eigc = eig(x1(vmin:vmax,1:dimpsi)'*x1(vmin:vmax,1:dimpsi));
eigd = eig(x1([vmin-1,vmax],1:dimpsi)'*x1([vmin-1,vmax],1:dimpsi));
eigt = eig(x1(vmin:vmax,1:dimpsi)'*x1(vmin:vmax,1:dimpsi) + x1([vmin-1,vmax],1:dimpsi)'*x1([vmin-1,vmax],1:dimpsi));
eigt2 = eig(x1(1:end,1:dimpsi)'*x1(1:end,1:dimpsi));
psi_eig = eig(x1(1:end-2,1:dimpsi)'*x1(1:end-2,1:dimpsi));
psi_m2 = max(x1(1:end-2,1:dimpsi))*max(x1(1:end-2,1:dimpsi))';
gamma_c_max = 1/(s*psi_m2);