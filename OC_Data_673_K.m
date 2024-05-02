%% CE40254
% Methanol-to-aromatics kinetic model at 673 K to determine rate constants 
% for olefin oligomerisation-cracking
% Written by Ben Mayoh
% Last modified: 01/04/24

clc
clear all
close all

%% Ethene

% Duration (s)
ti_C2 = 0; % Intial time
tf_C2 = 81; % Final time
tspan_C2 = [ti_C2, tf_C2]; % Time span

% Feed Mass Fractions
yi_C2_CH3OH = 0; % Methanol
yi_C2_C2H4 = 0.4; % Ethene
yi_C2_C3H6 = 0; % Propene
yi_C2_N2 = 1 - sum(yi_C2_CH3OH + yi_C2_C2H4 + yi_C2_C3H6); % Nitrogen
yi_C2 = [yi_C2_CH3OH, 0, 0, 0, 0, yi_C2_C2H4, 0, yi_C2_C3H6, 0, 0, 0, 0,...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, yi_C2_N2];...
    % All species

% ODE Solver
[t_C2, y_C2] = ode15s(@Kinetic_Model, tspan_C2, yi_C2);

% Modelled Data for Ethene at 673 K
y_C2_C2H4 = y_C2(:,6); % Ethene
y_C2_C3H6 = y_C2(:,8); % Propene
y_C2_C4H8 = y_C2(:,10); % Butene
y_C2_C5H10 = y_C2(:,12); % Pentene
y_C2_C6H12 = y_C2(:,15); % Hexene
y_C2_C7_plus = y_C2(:,18) + y_C2(:,21) + y_C2(:,24) + y_C2(:,27)...
    + y_C2(:,29) + y_C2(:,30); % C7+

% Experimental Data for ethene at 773 K
t_exp_C2 = [8.03, 11.92, 20.64, 31.12, 67.75]; % Contact time (s)
y_exp_C2_C2H4 = [39.341, 38.441, 38.195, 36.620, 26.424]/100; % Ethene
y_exp_C2_C3H6 = [0.054, 0.170, 0.357, 0.814, 3.812]/100; % Propene
y_exp_C2_C4H8 = [0.275, 0.341, 0.672, 0.940, 4.554]/100; % Butene
y_exp_C2_C5H10 = [0.040, 0.089, 0.273, 0.665, 2.721]/100; % Pentene
y_exp_C2_C6H12 = [0.112, 0.096, 0.232, 0.218, 0.906]/100; % Hexene
y_exp_C2_C7_plus = [0.000, 0.000, 0.000, 0.000, 0.198]/100; % C7+

% Plots
figure
plot(t_C2, y_C2_C2H4, 'b', t_exp_C2, y_exp_C2_C2H4, 'bx')
title('Ethene Mass Fraction with Ethene Feed at 673 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('Mod', 'Exp')

figure
hold on
plot(t_C2, y_C2_C3H6, 'r')
plot(t_exp_C2, y_exp_C2_C3H6, 'rx')
plot(t_C2, y_C2_C4H8, 'k')
plot(t_exp_C2, y_exp_C2_C4H8, 'kx')
plot(t_C2, y_C2_C5H10, 'c')
plot(t_exp_C2, y_exp_C2_C5H10, 'cx')
plot(t_C2, y_C2_C6H12, 'y')
plot(t_exp_C2, y_exp_C2_C6H12, 'yx')
plot(t_C2, y_C2_C7_plus, 'm')
plot(t_exp_C2, y_exp_C2_C7_plus, 'mx')
title('Olefin Mass Fractions with Ethene Feed at 673 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C3 - Mod', 'C3 - Exp', 'C4 - Mod', 'C4 - Exp', 'C5 - Mod',...
    'C5 - Exp', 'C6 - Mod', 'C6 - Exp', 'C7+ - Mod', 'C7+ - Exp')

% Error Calculation
sol = ode15s(@Kinetic_Model, tspan_C2, yi_C2); % ODE system solution

t0_mod_C2 = deval(sol, t_exp_C2(1))'; % Model solution at contact time 1
t1_mod_C2 = deval(sol, t_exp_C2(2))'; % Model solution at contact time 2
t2_mod_C2 = deval(sol, t_exp_C2(3))'; % Model solution at contact time 3
t3_mod_C2 = deval(sol, t_exp_C2(4))'; % Model solution at contact time 4
t4_mod_C2 = deval(sol, t_exp_C2(5))'; % Model solution at contact time 5

y_t0_mod_C2 = [t0_mod_C2(6), t0_mod_C2(8), t0_mod_C2(10), t0_mod_C2(12),...
    t0_mod_C2(15), (t0_mod_C2(18) + t0_mod_C2(21) + t0_mod_C2(24)...
    + t0_mod_C2(27) + t0_mod_C2(29) + t0_mod_C2(30))]; % Model mass...
    % fractions at contact time 1
y_t1_mod_C2 = [t1_mod_C2(6), t1_mod_C2(8), t1_mod_C2(10), t1_mod_C2(12),...
    t1_mod_C2(15), (t1_mod_C2(18) + t1_mod_C2(21) + t1_mod_C2(24)...
    + t1_mod_C2(27) + t1_mod_C2(29) + t1_mod_C2(30))]; % Model mass..
    % fractions at contact time 2
y_t2_mod_C2 = [t2_mod_C2(6), t2_mod_C2(8), t2_mod_C2(10), t2_mod_C2(12),...
    t2_mod_C2(15), (t2_mod_C2(18) + t2_mod_C2(21) + t2_mod_C2(24)...
    + t2_mod_C2(27) + t2_mod_C2(29) + t2_mod_C2(30))]; % Model mass...
    % fractions at contact time 3
y_t3_mod_C2 = [t3_mod_C2(6), t3_mod_C2(8), t3_mod_C2(10), t3_mod_C2(12),...
    t3_mod_C2(15), (t3_mod_C2(18) + t3_mod_C2(21) + t3_mod_C2(24)...
    + t3_mod_C2(27) + t3_mod_C2(29) + t3_mod_C2(30))]; % Model mass...
    % fractions at contact time 4
y_t4_mod_C2 = [t4_mod_C2(6), t4_mod_C2(8), t4_mod_C2(10), t4_mod_C2(12),...
    t4_mod_C2(15), (t4_mod_C2(18) + t4_mod_C2(21) + t4_mod_C2(24)...
    + t4_mod_C2(27) + t4_mod_C2(29) + t4_mod_C2(30))]; % Model mass...
    % fractions at contact time 5

mat_mod_C2 = [y_t0_mod_C2; y_t1_mod_C2; y_t2_mod_C2;...
    y_t3_mod_C2; y_t4_mod_C2]; % Modelled mass fraction values matrix
mat_exp_C2 = [y_exp_C2_C2H4', y_exp_C2_C3H6', y_exp_C2_C4H8',...
    y_exp_C2_C5H10', y_exp_C2_C6H12', y_exp_C2_C7_plus'];...
    % Experimental mass fraction values matrix

SSE_C2 = (mat_mod_C2 - mat_exp_C2).^2; % Squared errors

WSSE_C2_C2H4 = (SSE_C2(:,1))/(sum(mat_exp_C2(:,1))); % Weight for ethene
WSSE_C2_C3H6 = (SSE_C2(:,2))/(sum(mat_exp_C2(:,2))); % Weight for propene
WSSE_C2_C4H8 = (SSE_C2(:,3))/(sum(mat_exp_C2(:,3))); % Weight for butene
WSSE_C2_C5H10 = (SSE_C2(:,4))/(sum(mat_exp_C2(:,4))); % Weight for pentene
WSSE_C2_C6H12 = (SSE_C2(:,5))/(sum(mat_exp_C2(:,5))); % Weight for hexene
WSSE_C2_C7_plus = (SSE_C2(:,6))/(sum(mat_exp_C2(:,6))); % Weight for C7+

mat_WSSE_C2 = [WSSE_C2_C2H4, WSSE_C2_C3H6, WSSE_C2_C4H8, WSSE_C2_C5H10,...
    WSSE_C2_C6H12, WSSE_C2_C7_plus]; % Weighted sum of squared errors...
    % matrix

WSSE_C2 = sum(sum(mat_WSSE_C2))*100; % Weighted sum of sqaured errors...
    % (multiplied by 100 to avoid rounding errors)

fprintf('Ethene WSSE = %.3f\n', WSSE_C2) % Error outputted

% Writing to Excel File
data_C2_mod = table(t_C2, y_C2_C2H4, y_C2_C3H6, y_C2_C4H8, y_C2_C5H10,...
    y_C2_C6H12, y_C2_C7_plus);
data_C2_exp = table(t_exp_C2', y_exp_C2_C2H4', y_exp_C2_C3H6',...
    y_exp_C2_C4H8', y_exp_C2_C5H10', y_exp_C2_C6H12',...
    y_exp_C2_C7_plus');

writetable(data_C2_mod, 'Oligomerisation-Cracking_Data.xlsx', 'Sheet',...
    'C2 Mod 673 K')
writetable(data_C2_exp,'Oligomerisation-Cracking_Data.xlsx', 'Sheet',...
    'C2 Exp 673 K')