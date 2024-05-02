%% CE40254
% Methanol-to-aromatics kinetic model at 643 K to determine rate constants 
% for olefin methylation
% Written by Ben Mayoh
% Last modified: 24/04/24

clc
clear all
close all

%% Ethene

% Duration (s)
ti = 0; % Intial time
tf = 2.5; % Final time
tspan = [ti, tf]; % Time span

% Feed Mass Fractions
yi_CH3OH = 0.25; % Methanol
yi_C2H4 = 0; % Ethene
yi_C3H6 = 0; % Propene
yi_N2 = 1 - sum(yi_CH3OH + yi_C2H4 + yi_C3H6); % Nitrogen
yi = [yi_CH3OH, 0, 0, 0, 0, yi_C2H4, 0, yi_C3H6, 0, 0, 0, 0, 0, 0, 0, 0,...
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, yi_N2]; % All species

% ODE Solver
[t, y] = ode15s(@Kinetic_Model, tspan, yi);

% Modelled Data for Ethene at 773 K
y_CH3OH = y(:,1); % Methanol
y_CH3OCH3 = y(:,2); % DME
y_H2O = y(:,3); % Water
y_CO = y(:,4); % Carbon Monoxide
y_CH4 = y(:,5); % Methane
y_C2H4 = y(:,6); % Ethene
y_C2H6 = y(:,7); % Ethane
y_C3H6 = y(:,8); % Propene
y_C3H8 = y(:,9); % Propane
y_C4H8 = y(:,10); % Butene
y_C4H10 = y(:,11); % Butane
y_C5H10 = y(:,12); % Pentene
y_C5H12 = y(:,13); % Pentane
y_C6H6 = y(:,14); % Benzene
y_C6H12 = y(:,15); % Hexene
y_C6H14 = y(:,16); % Hexane
y_C7H8 = y(:,17); % Toluene
y_C7_plus = y(:,18) + y(:,19) + y(:,21) + y(:,22) + y(:,24) + y(:,25)...
    + y(:,27) + y(:,28) + y(:,29) + y(:,30); % C7+
y_C8H10 = y(:,20); % Ethylbenzene
y_C9H12 = y(:,23); % Propylbenzene
y_C10H14 = y(:,26); % Butylbenzene
y_H2 = y(:,31); % Hydrogen
y_N2 = y(:,32); % Nitrogen

% Experimental Data for MtH at 643 K
t_exp = [0.46, 0.63, 0.85, 1.17, 2.12]; % Contact time (s)
y_exp_CH3OH = [15.892, 14.338, 12.689, 12.704, 10.873]/100; % Methanol
y_exp_CH3OCH3 = [6.719, 7.694, 8.936, 8.638, 7.668]/100; % DME
y_exp_H2O = [2.670, 3.143, 3.803, 3.833, 5.004]/100; % Water
y_exp_CO = [0.036, 0.041, 0.062, 0.093, 0.128]/100; % Carbon Monoxide
y_exp_CH4 = [0.026, 0.028, 0.042, 0.067, 0.087]/100; % Methane
y_exp_C2H4 = [0.020, 0.058, 0.097, 0.148, 0.402]/100; % Ethene
y_exp_C2H6 = [0.000, 0.000, 0.000, 0.000, 0.000]/100; % Ethane
y_exp_C3H6 = [0.017, 0.042, 0.097, 0.218, 0.634]/100; % Propene
y_exp_C3H8 = [0.000, 0.000, 0.000, 0.000, 0.000]/100; % Propane
y_exp_C4H8 = [0.000, 0.006, 0.028, 0.066, 0.211]/100; % Butene
y_exp_C4H10 = [0.000, 0.000, 0.000, 0.021, 0.062]/100; % Butane
y_exp_C5H10 = [0.000, 0.000, 0.005, 0.030, 0.064]/100; % Pentene
y_exp_C5H12 = [0.000, 0.000, 0.000, 0.018, 0.048]/100; % Pentane
y_exp_C6H6 = [0.000, 0.000, 0.000, 0.000, 0.000]/100; % Benzene
y_exp_C6H12 = [0.000, 0.000, 0.000, 0.015, 0.041]/100; % Hexene
y_exp_C6H14 = [0.000, 0.000, 0.000, 0.007, 0.029]/100; % Hexane
y_exp_C7H8 = [0.000, 0.000, 0.000, 0.000, 0.000]/100; % Toluene
y_exp_C7_plus = [0.000, 0.000, 0.000, 0.000, 0.000]/100; % C7+
y_exp_C8H10 = [0.000, 0.000, 0.000, 0.007, 0.057]/100; % Ethylbenzene
y_exp_C9H12 = [0.000, 0.000, 0.000, 0.000, 0.016]/100; % Propylbenzene
y_exp_C10H14 = [0.000, 0.000, 0.000, 0.000, 0.000]/100; % Butylbenzene

% Plots
figure
hold on
plot(t, y_CH3OH, 'k')
plot(t_exp, y_exp_CH3OH, 'kx')
plot(t, y_CH3OCH3, 'b')
plot(t_exp, y_exp_CH3OCH3, 'bx')
plot(t, y_H2O, 'r')
plot(t_exp, y_exp_H2O, 'rx')
plot(t, y_CO, 'g')
plot(t_exp, y_exp_CO, 'gx')
plot(t, y_CH4, 'm')
plot(t_exp, y_exp_CH4, 'mx')
title({'Methanol, DME, Water, Carbon Monoxide, and Methane',...
    'Mass Fractions with Methanol Feed at 643 K'})
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('CH3OH - Mod', 'CH3OH - Exp', 'CH3OCH3 - Mod', 'CH3OCH3 - Exp',...
    'H2O - Mod', 'H2O - Exp', 'CO - Mod', 'CO - Exp', 'CH4 - Mod',...
    'CH4 - Exp')

figure
hold on
plot(t, y_C2H4, 'b')
plot(t_exp, y_exp_C2H4, 'bx')
plot(t, y_C3H6, 'r')
plot(t_exp, y_exp_C3H6, 'rx')
plot(t, y_C4H8, 'k')
plot(t_exp, y_exp_C4H8, 'kx')
plot(t, y_C5H10, 'c')
plot(t_exp, y_exp_C5H10, 'cx')
plot(t, y_C6H12, 'm')
plot(t_exp, y_exp_C6H12, 'mx')
plot(t, y_C7_plus, 'y')
plot(t_exp, y_exp_C7_plus, 'yx')
title('Olefin Mass Fractions with Methanol Feed at 643 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C2H4 - Mod', 'C2H4 - Exp', 'C3H6 - Mod', 'C3H6 - Exp',...
    'C4H8 - Mod', 'C4H8 - Exp', 'C5H10 - Mod', 'C5H10 - Exp',...
    'C6H12 - Mod', 'C6H12 - Exp', 'C7+ - Mod', 'C7+ - Exp')

% Error Calculation
sol = ode15s(@Kinetic_Model, tspan, yi); % ODE system solution

t0_mod = deval(sol, t_exp(1))'; % Model solution at contact time 1
t1_mod = deval(sol, t_exp(2))'; % Model solution at contact time 2
t2_mod = deval(sol, t_exp(3))'; % Model solution at contact time 3
t3_mod = deval(sol, t_exp(4))'; % Model solution at contact time 4
t4_mod = deval(sol, t_exp(5))'; % Model solution at contact time 5

y_t0_mod = [t0_mod(1:6), t0_mod(8), t0_mod(10), t0_mod(12), t0_mod(15)];...
    % Model mass fractions at contact time 1
y_t1_mod = [t1_mod(1:6), t1_mod(8), t1_mod(10), t1_mod(12), t1_mod(15)];...
    % Model mass fractions at contact time 2
y_t2_mod = [t2_mod(1:6), t2_mod(8), t2_mod(10), t2_mod(12), t2_mod(15)];...
    % Model mass fractions at contact time 3
y_t3_mod = [t3_mod(1:6), t3_mod(8), t3_mod(10), t3_mod(12), t3_mod(15)];...
    % Model mass fractions at contact time 4
y_t4_mod = [t4_mod(1:6), t4_mod(8), t4_mod(10), t4_mod(12), t4_mod(15)];...
    % Model mass fractions at contact time 5

mat_mod = [y_t0_mod; y_t1_mod; y_t2_mod; y_t3_mod; y_t4_mod]; % Modelled...
    % mass fraction values matrix
mat_exp = [y_exp_CH3OH', y_exp_CH3OCH3', y_exp_H2O', y_exp_CO',...
    y_exp_CH4', y_exp_C2H4', y_exp_C3H6', y_exp_C4H8', y_exp_C5H10',...
    y_exp_C6H12']; % Experimental mass fraction values matrix

SSE = (mat_mod - mat_exp).^2; % Squared errors

WSSE_CH3OH = (SSE(:,1))/(sum(mat_exp(:,1))); % Weight for Methanol
WSSE_CH3OCH3 = (SSE(:,2))/(sum(mat_exp(:,2))); % Weight for DME
WSSE_H2O = (SSE(:,3))/(sum(mat_exp(:,3))); % Weight for water
WSSE_CO = (SSE(:,4))/(sum(mat_exp(:,4))); % Weight for carbon monoxide
WSSE_CH4 = (SSE(:,5))/(sum(mat_exp(:,5))); % Weight for methane
WSSE_C2H4 = (SSE(:,6))/(sum(mat_exp(:,6))); % Weight for ethene
WSSE_C3H6 = (SSE(:,7))/(sum(mat_exp(:,7))); % Weight for propene
WSSE_C4H8 = (SSE(:,8))/(sum(mat_exp(:,8))); % Weight for butene
WSSE_C5H10 = (SSE(:,9))/(sum(mat_exp(:,9))); % Weight for pentene
WSSE_C6H12 = (SSE(:,10))/(sum(mat_exp(:,10))); % Weight for hexene

mat_WSSE = [WSSE_CH3OH, WSSE_CH3OCH3, WSSE_H2O, WSSE_CO, WSSE_CH4,...
    WSSE_C2H4, WSSE_C3H6, WSSE_C4H8, WSSE_C5H10, WSSE_C6H12]; % Weighted...
    % sum of squared errors matrix

WSSE = sum(sum(mat_WSSE))*100; % Weighted sum of sqaured errors...
    % (multiplied by 100 to avoid rounding errors)

fprintf('WSSE = %.4f\n', WSSE) % Error outputted

% Writing to Excel File
data_mod = table(t, y_CH3OH, y_CH3OCH3, y_H2O, y_CO, y_CH4, y_C2H4,...
    y_C2H6, y_C3H6, y_C3H8, y_C4H8, y_C4H10, y_C5H10, y_C5H12, y_C6H6,...
    y_C6H12, y_C6H14, y_C7H8, y_C7_plus, y_C8H10, y_C9H12, y_C10H14);
data_exp = table(t_exp', y_exp_CH3OH', y_exp_CH3OCH3', y_exp_H2O',...
    y_exp_CO', y_exp_CH4', y_exp_C2H4', y_exp_C2H6', y_exp_C3H6',...
    y_exp_C3H8', y_exp_C4H8', y_exp_C4H10', y_exp_C5H10', y_exp_C5H12',...
    y_exp_C6H6', y_exp_C6H12', y_exp_C6H14', y_exp_C7H8',...
    y_exp_C7_plus', y_exp_C8H10', y_exp_C9H12', y_exp_C10H14');

writetable(data_mod, 'Methanol-to-Hydrocarbons_Data.xlsx', 'Sheet',...
    'Mod (3) 643 K')
writetable(data_exp,'Methanol-to-Hydrocarbons_Data.xlsx', 'Sheet',...
    'Exp (3) 643 K')