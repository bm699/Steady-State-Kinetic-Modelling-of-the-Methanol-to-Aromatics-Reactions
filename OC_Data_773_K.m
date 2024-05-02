%% CE40254
% Methanol-to-aromatics kinetic model at 773 K to determine rate constants 
% for olefin oligomerisation-cracking
% Written by Ben Mayoh
% Last modified: 01/04/24

clc
clear all
close all

%% Ethene

% Duration (s)
ti_C2 = 0; % Intial time
tf_C2 = 57; % Final time
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

% Modelled Data for ethene at 773 K
y_C2_C2H4 = y_C2(:,6); % Ethene
y_C2_C3H6 = y_C2(:,8); % Propene
y_C2_C4H8 = y_C2(:,10); % Butene
y_C2_C5H10 = y_C2(:,12); % Pentene
y_C2_C6H12 = y_C2(:,15); % Hexene
y_C2_C7_plus = y_C2(:,18) + y_C2(:,21) + y_C2(:,24) + y_C2(:,27)...
    + y_C2(:,29) + y_C2(:,30); % C7+

% Experimental Data for Ethene at 773 K
t_exp_C2 = [1.742, 6.722, 16.688, 28.963, 47.505]; % Contact time (s)
y_exp_C2_C2H4 = [39.912, 39.398, 38.875, 35.047, 31.407]/100; % Ethene
y_exp_C2_C3H6 = [0.007, 0.108, 0.660, 1.856, 4.407]/100; % Propene
y_exp_C2_C4H8 = [0.162, 0.497, 1.135, 1.611, 2.765]/100; % Butene
y_exp_C2_C5H10 = [0.001, 0.061, 0.345, 0.747, 1.620]/100; % Pentene
y_exp_C2_C6H12 = [0.011, 0.049, 0.141, 0.232, 0.404]/100; % Hexene
y_exp_C2_C7_plus = [0.000, 0.000, 0.000, 0.000, 0.067]/100; % C7+

% Plots
figure
hold on
plot(t_C2, y_C2_C2H4, 'b')
plot(t_exp_C2, y_exp_C2_C2H4, 'bx')
title('Ethene Mass Fraction with Ethene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C2 - Mod', 'C2 - Exp')

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
title('Olefin Mass Fractions with Ethene Feed at 773 K')
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
    + t1_mod_C2(27) + t1_mod_C2(29) + t1_mod_C2(30))]; % Model mass...
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
    'C2 Mod 773 K')
writetable(data_C2_exp,'Oligomerisation-Cracking_Data.xlsx', 'Sheet',...
    'C2 Exp 773 K')

%% Propene

% Duration (s)
ti_C3 = 0; % Intial time
tf_C3 = 6; % Final time
tspan_C3 = [ti_C3, tf_C3]; % Time span

% Feed Mass Fractions
yi_C3_CH3OH = 0; % Methanol
yi_C3_C2H4 = 0; % Ethene
yi_C3_C3H6 = 0.5; % Propene
yi_C3_N2 = 1 - sum(yi_C3_CH3OH + yi_C3_C2H4 + yi_C3_C3H6); % Nitrogen
yi_C3 = [yi_C3_CH3OH, 0, 0, 0, 0, yi_C3_C2H4, 0, yi_C3_C3H6, 0, 0, 0, 0,...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, yi_C3_N2];...
    % All species

% ODE Solver
[t_C3, y_C3] = ode15s(@Kinetic_Model, tspan_C3, yi_C3);

% Modelled Data for Propene at 773 K
y_C3_C2H4 = y_C3(:,6); % Ethene
y_C3_C3H6 = y_C3(:,8); % Propene
y_C3_C4H8 = y_C3(:,10); % Butene
y_C3_C5H10 = y_C3(:,12); % Pentene
y_C3_C6H12 = y_C3(:,15); % Hexene
y_C3_C7_plus = y_C3(:,18) + y_C3(:,21) + y_C3(:,24) + y_C3(:,27)...
    + y_C3(:,29) + y_C3(:,30); % C7+

% Experimental Data for Propene at 773 K
t_exp_C3 = [0.273, 0.469, 0.934, 2.039, 4.973]; % Contact time (s)
y_exp_C3_C2H4 = [0.025, 0.106, 0.423, 0.519, 1.718]/100; % Ethene
y_exp_C3_C3H6 = [46.237, 43.710, 40.040, 34.424, 27.009]/100; % Propene
y_exp_C3_C4H8 = [0.192, 0.766, 1.847, 5.910, 12.211]/100; % Butene
y_exp_C3_C5H10 = [0.165, 0.421, 0.882, 2.615, 4.743]/100; % Pentene
y_exp_C3_C6H12 = [3.374, 4.928, 6.549, 6.140, 3.834]/100; % Hexene
y_exp_C3_C7_plus = [0.043, 0.155, 0.287, 0.491, 0.590]/100; % C7+

% Plots
figure
hold on
plot(t_C3, y_C3_C3H6, 'r')
plot(t_exp_C3, y_exp_C3_C3H6, 'rx')
title('Propene Mass Fraction with Propene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C3 - Mod', 'C3 - Exp')

figure
hold on
plot(t_C3, y_C3_C2H4, 'b')
plot(t_exp_C3, y_exp_C3_C2H4, 'bx')
plot(t_C3, y_C3_C4H8, 'k')
plot(t_exp_C3, y_exp_C3_C4H8, 'kx')
plot(t_C3, y_C3_C5H10, 'c')
plot(t_exp_C3, y_exp_C3_C5H10, 'cx')
plot(t_C3, y_C3_C6H12, 'y')
plot(t_exp_C3, y_exp_C3_C6H12, 'yx')
plot(t_C3, y_C3_C7_plus, 'm')
plot(t_exp_C3, y_exp_C3_C7_plus, 'mx')
title('Olefin Mass Fractions with Propene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C2 - Mod', 'C2 - Exp', 'C4 - Mod', 'C4 - Exp', 'C5 - Mod',...
    'C5 - Exp', 'C6 - Mod', 'C6 - Exp', 'C7+ - Mod', 'C7+ - Exp')

% Error Calculation
sol = ode15s(@Kinetic_Model, tspan_C3, yi_C3); % ODE system solution

t0_mod_C3 = deval(sol, t_exp_C3(1))'; % Model solution at contact time 1
t1_mod_C3 = deval(sol, t_exp_C3(2))'; % Model solution at contact time 2
t2_mod_C3 = deval(sol, t_exp_C3(3))'; % Model solution at contact time 3
t3_mod_C3 = deval(sol, t_exp_C3(4))'; % Model solution at contact time 4
t4_mod_C3 = deval(sol, t_exp_C3(5))'; % Model solution at contact time 5

y_t0_mod_C3 = [t0_mod_C3(6), t0_mod_C3(8), t0_mod_C3(10), t0_mod_C3(12),...
    t0_mod_C3(15), (t0_mod_C3(18) + t0_mod_C3(21) + t0_mod_C3(24)...
    + t0_mod_C3(27) + t0_mod_C3(29) + t0_mod_C3(30))]; % Model mass...
    % fractions at contact time 1
y_t1_mod_C3 = [t1_mod_C3(6), t1_mod_C3(8), t1_mod_C3(10), t1_mod_C3(12),...
    t1_mod_C3(15), (t1_mod_C3(18) + t1_mod_C3(21) + t1_mod_C3(24)...
    + t1_mod_C3(27) + t1_mod_C3(29) + t1_mod_C3(30))]; % Model mass...
    % fractions at contact time 2
y_t2_mod_C3 = [t2_mod_C3(6), t2_mod_C3(8), t2_mod_C3(10), t2_mod_C3(12),...
    t2_mod_C3(15), (t2_mod_C3(18) + t2_mod_C3(21) + t2_mod_C3(24)...
    + t2_mod_C3(27) + t2_mod_C3(29) + t2_mod_C3(30))]; % Model mass...
    % fractions at contact time 3
y_t3_mod_C3 = [t3_mod_C3(6), t3_mod_C3(8), t3_mod_C3(10), t3_mod_C3(12),...
    t3_mod_C3(15), (t3_mod_C3(18) + t3_mod_C3(21) + t3_mod_C3(24)...
    + t3_mod_C3(27) + t3_mod_C3(29) + t3_mod_C3(30))]; % Model mass...
    % fractions at contact time 4
y_t4_mod_C3 = [t4_mod_C3(6), t4_mod_C3(8), t4_mod_C3(10), t4_mod_C3(12),...
    t4_mod_C3(15), (t4_mod_C3(18) + t4_mod_C3(21) + t4_mod_C3(24)...
    + t4_mod_C3(27) + t4_mod_C3(29) + t4_mod_C3(30))]; % Model mass...
    % fractions at contact time 5

mat_mod_C3 = [y_t0_mod_C3; y_t1_mod_C3; y_t2_mod_C3;...
    y_t3_mod_C3; y_t4_mod_C3]; % Modelled mass fraction values matrix
mat_exp_C3 = [y_exp_C3_C2H4', y_exp_C3_C3H6', y_exp_C3_C4H8',...
    y_exp_C3_C5H10', y_exp_C3_C6H12', y_exp_C3_C7_plus'];...
    % Experimental mass fraction values matrix

SSE_C3 = (mat_mod_C3 - mat_exp_C3).^2; % Squared errors

WSSE_C3_C2H4 = (SSE_C3(:,1))/(sum(mat_exp_C3(:,1))); % Weight for ethene
WSSE_C3_C3H6 = (SSE_C3(:,2))/(sum(mat_exp_C3(:,2))); % Weight for propene
WSSE_C3_C4H8 = (SSE_C3(:,3))/(sum(mat_exp_C3(:,3))); % Weight for butene
WSSE_C3_C5H10 = (SSE_C3(:,4))/(sum(mat_exp_C3(:,4))); % Weight for pentene
WSSE_C3_C6H12 = (SSE_C3(:,5))/(sum(mat_exp_C3(:,5))); % Weight for hexene
WSSE_C3_C7_plus = (SSE_C3(:,6))/(sum(mat_exp_C3(:,6))); % Weight for C7+

mat_WSSE_C3 = [WSSE_C3_C2H4, WSSE_C3_C3H6, WSSE_C3_C4H8, WSSE_C3_C5H10,...
    WSSE_C3_C6H12, WSSE_C3_C7_plus]; % Weighted sum of squared errors...
    % matrix

WSSE_C3 = sum(sum(mat_WSSE_C3))*100; % Weighted sum of sqaured errors...
    % (multiplied by 100 to avoid rounding errors)

fprintf('Propene WSSE = %.3f\n', WSSE_C3) % Error outputted

% Writing to Excel File
data_C3_mod = table(t_C3, y_C3_C2H4, y_C3_C3H6, y_C3_C4H8, y_C3_C5H10,...
    y_C3_C6H12, y_C3_C7_plus);
data_C3_exp = table(t_exp_C3', y_exp_C3_C2H4', y_exp_C3_C3H6',...
    y_exp_C3_C4H8', y_exp_C3_C5H10', y_exp_C3_C6H12',...
    y_exp_C3_C7_plus');

writetable(data_C3_mod, 'Oligomerisation-Cracking_Data.xlsx', 'Sheet',...
    'C3 Mod 773 K')
writetable(data_C3_exp,'Oligomerisation-Cracking_Data.xlsx', 'Sheet',...
    'C3 Exp 773 K')