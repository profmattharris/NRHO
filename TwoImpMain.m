clear all; close all; clc

%--------------------------------------------------------------------------
% This script specifies an initial guess type and point in NRHO and then
% calls a function to generate locally optimal solutions.
% Inputs:
%   1) guess_flag (1 for Hohmann guess, 2 for two-body optimal guess)
%   2) rev (vector of numbers between 1 and 2 indicating points in NRHO)
%
% There are no outputs because this is a script. 
%--------------------------------------------------------------------------

%--------%
% INPUTS %
%--------%
flag = 2;
rev  = [1:.1:1.9];


%----------------%
% FUNCTION CALLS %
%----------------%
for i = 1:length(rev)
    [dv2b(i),dv3b(i),tflight(i)] = TwoImpOptimizer( rev(i), 2 );
end


%--------%
% GRAPHS %
%--------%
figure, plot(rev-1,1000*dv2b,'linewidth',2), grid on, hold on
plot(rev-1,1000*dv3b,'linewidth',2)
xlabel('Rev Fraction'), ylabel('delta-v (m/s)')
legend('2BP','3BP')

figure, plot(rev-1,tflight/3600,'linewidth',2), grid on
xlabel('Rev Fraction'), ylabel('Transfer Time (hr)')

figure
yyaxis left
plot(rev-1,1000*dv2b,'linewidth',2), hold on
plot(rev-1,1000*dv3b,'linewidth',2), grid on
xlabel('Rev Fraction')
ylabel('Delta-v (m/s)')
yyaxis right
plot(rev-1,tflight/3600,'linewidth',2)
ylabel('Flight Time (hr)')