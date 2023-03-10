% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: pickexample.m
%--------------------------------------------------------------------------
% Project: Discretized hybrid gradient descent algorithm (DiscHyGradAlg)
% Description: function to pick example to simulate, input case number to
% output example parameters
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision:

function [s, gamma_c1, gamma_d1, psi0, A, B, H, J, u, C, D, TSPAN, JSPAN] = pickexample(x)
%% Example System Choice

switch x
    case 1 %bouncing ball, continuous GD
        s = 0.01; psi0 = [3 0]; gamma_c1 = 1/(s*2000); gamma_d1 = 0;
        TSPAN = [0 5]; JSPAN = [0 TSPAN(end)/s];
        A = [0 1; 0 0]; B = [0; -9.8];
        H = [0 0; 0 -1]; J = [0; 0];
        u = @(t) 1;
        C = @(psi,t) psi(1) >= 0;
        D = @(psi,t) (psi(1) <= 0 && psi(2) <= 0);

    case 2 %bouncing ball, discrete GD
        s = 0.01; psi0 = [3 0]; gamma_c1 = 0; gamma_d1 = 1/(500*s);
        TSPAN = [0 5]; JSPAN = [0 TSPAN(end)/s];
        A = [0 1; 0 0]; B = [0; -9.8];
        H = [0 0; 0 -1]; J = [0; 0];
        u = @(t) 1;
        C = @(psi,t) psi(1) >= 0;
        D = @(psi,t) (psi(1) <= 0 && psi(2) <= 0);

    case 3 %bouncing ball, hybrid GD
        s = 0.01; psi0 = [3 0]; gamma_c1 = 1/(s*2000); gamma_d1 = 1/(500*s);
        TSPAN = [0 5]; JSPAN = [0 TSPAN(end)/s];
        A = [0 1; 0 0]; B = [0; -9.8];
        H = [0 0; 0 -1]; J = [0; 0];
        u = @(t) 1;
        C = @(psi,t) psi(1) >= 0;
        D = @(psi,t) (psi(1) <= 0 && psi(2) <= 0);
end