% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: plots.m
%--------------------------------------------------------------------------
% Project: Discretized hybrid gradient descent algorithm (DiscHyGradAlg)
% Description: Create figures for trajectory of error norm to compare
% discretized continuous, discretized, and discretized hybrid GD algorithms
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision:

clc; clear all; close all;

[t1c, j1c, errornormc] = run_function(1);
[t1d, j1d, errornormd] = run_function(2);
[t1cd, j1cd, errornormcd] = run_function(3);

modificatorF{1} = 'c';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1;

modificatorF2{1} = 'r';
modificatorF2{2} = 'LineWidth';
modificatorF2{3} = 1;

modificatorF3{1} = 'b';
modificatorF3{2} = 'LineWidth';
modificatorF3{3} = 1;

modificatorJ{1} = 'c:';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 2;
modificatorJ{4} = 'Marker';
modificatorJ{5} = '.';
modificatorJ{6} = 'MarkerEdgeColor';
modificatorJ{7} = 'c';
modificatorJ{8} = 'MarkerFaceColor';
modificatorJ{9} = 'c';
modificatorJ{10} = 'MarkerSize';
modificatorJ{11} = 2;

modificatorJ2{1} = 'r:';
modificatorJ2{2} = 'LineWidth';
modificatorJ2{3} = 2;
modificatorJ2{4} = 'Marker';
modificatorJ2{5} = '.';
modificatorJ2{6} = 'MarkerEdgeColor';
modificatorJ2{7} = 'r';
modificatorJ2{8} = 'MarkerFaceColor';
modificatorJ2{9} = 'r';
modificatorJ2{10} = 'MarkerSize';
modificatorJ2{11} = 2;

modificatorJ3{1} = 'b:';
modificatorJ3{2} = 'LineWidth';
modificatorJ3{3} = 2;
modificatorJ3{4} = 'Marker';
modificatorJ3{5} = '.';
modificatorJ3{6} = 'MarkerEdgeColor';
modificatorJ3{7} = 'b';
modificatorJ3{8} = 'MarkerFaceColor';
modificatorJ3{9} = 'b';
modificatorJ3{10} = 'MarkerSize';
modificatorJ3{11} = 2;

modificatorJ4{1} = 'g:';
modificatorJ4{2} = 'LineWidth';
modificatorJ4{3} = 2;
modificatorJ4{4} = 'Marker';
modificatorJ4{5} = '.';
modificatorJ4{6} = 'MarkerEdgeColor';
modificatorJ4{7} = 'g';
modificatorJ4{8} = 'MarkerFaceColor';
modificatorJ4{9} = 'g';
modificatorJ4{10} = 'MarkerSize';
modificatorJ4{11} = 2;

figure(3)
%clf;
tlt = tiledlayout(2,1);
nexttile(1)
hold on; plot(0,'c'); plot(0,'r'); plot(0,'b'); plot(0,'g'); hold off;
legend('Discretized Continuous GD','Discrete GD','Discretized Hybrid GD','AutoUpdate','off'); hold on;
plotHarc(t1c(1:end-1),j1c(1:end-1),errornormc,[],modificatorF,modificatorJ); hold on
plotHarc(t1d(1:end-1),j1d(1:end-1),errornormd,[],modificatorF2,modificatorJ2); hold on
plotHarc(t1cd(1:end-1),j1cd(1:end-1),errornormcd,[],modificatorF3,modificatorJ3); hold on
xlabel('$t$ [s]','interpreter','latex','FontSize',18);
ylabel('$|\tilde\theta_s|$','interpreter','latex','FontSize',18);
grid on; box on;
tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1*pos(3), 1*pos(4)])