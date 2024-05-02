%% CE40254
% Methanol-to-aromatics kinetic model at 643 K to determine rate constants 
% for methanol-to-hydrocarbon reactions
% Written by Ben Mayoh
% Last modified: 24/04/24

clc
clear all
close all

%% Ethene

% Duration (s)
ti = 0; % Intial time
tf = 24; % Final time
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
t_exp = [0.46, 0.63, 0.85, 1.17, 2.12, 4.35, 8.80, 11.53, 20.02];...
    % Contact time (s)
y_exp_CH3OH = [15.892, 14.338, 12.689, 12.704, 10.873, 6.579, 2.944,...
    1.772, 0.000]/100; % Methanol
y_exp_CH3OCH3 = [6.719, 7.694, 8.936, 8.638, 7.668, 6.325, 1.445, 0.696,...
    0.000]/100; % DME
y_exp_H2O = [2.670, 3.143, 3.803, 3.833, 5.004, 7.710, 10.864, 11.578,...
    13.866]/100; % Water
y_exp_CO = [0.036, 0.041, 0.062, 0.093, 0.128, 0.229, 0.373, 0.435,...
    0.447]/100; % Carbon Monoxide
y_exp_CH4 = [0.026, 0.028, 0.042, 0.067, 0.087, 0.156, 0.248, 0.288,...
    0.296]/100; % Methane
y_exp_C2H4 = [0.020, 0.058, 0.097, 0.148, 0.402, 1.029, 1.511, 1.800,...
    1.423]/100; % Ethene
y_exp_C2H6 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.012,...
    0.024]/100; % Ethane
y_exp_C3H6 = [0.017, 0.042, 0.097, 0.218, 0.634, 1.555, 2.510, 2.726,...
    2.235]/100; % Propene
y_exp_C3H8 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.033, 0.123, 0.203,...
    0.551]/100; % Propane
y_exp_C4H8 = [0.000, 0.006, 0.028, 0.066, 0.211, 0.517, 1.255, 1.417,...
    1.104]/100; % Butene
y_exp_C4H10 = [0.000, 0.000, 0.000, 0.021, 0.062, 0.199, 0.485, 0.637,...
    1.358]/100; % Butane
y_exp_C5H10 = [0.000, 0.000, 0.005, 0.030, 0.064, 0.226, 0.561, 0.688,...
    0.519]/100; % Pentene
y_exp_C5H12 = [0.000, 0.000, 0.000, 0.018, 0.048, 0.140, 0.353, 0.403,...
    0.889]/100; % Pentane
y_exp_C6H6 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.005, 0.021,...
    0.041]/100; % Benzene
y_exp_C6H12 = [0.000, 0.000, 0.000, 0.015, 0.041, 0.105, 0.404, 0.256,...
    0.187]/100; % Hexene
y_exp_C6H14 = [0.000, 0.000, 0.000, 0.007, 0.029, 0.081, 0.252, 0.184,...
    0.342]/100; % Hexane
y_exp_C7H8 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.003, 0.034, 0.041,...
    0.157]/100; % Toluene
y_exp_C7_plus = [0.000, 0.000, 0.000, 0.000, 0.000, 0.015, 0.102, 0.134,...
    0.118]/100; % C7+
y_exp_C8H10 = [0.000, 0.000, 0.000, 0.007, 0.057, 0.160, 0.393, 0.461,...
    1.012]/100; % Ethylbenzene
y_exp_C9H12 = [0.000, 0.000, 0.000, 0.000, 0.016, 0.035, 0.131, 0.167,...
    0.338]/100; % Propylbenzene
y_exp_C10H14 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.011, 0.024,...
    0.039]/100; % Butylbenzene

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
title('Ethene and Propene Mass Fractions with Methanol Feed at 643 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C2H4 - Mod', 'C2H4 - Exp', 'C3H6 - Mod', 'C3H6 - Exp')

figure
hold on
plot(t, y_C4H8, 'k')
plot(t_exp, y_exp_C4H8, 'kx')
plot(t, y_C5H10, 'c')
plot(t_exp, y_exp_C5H10, 'cx')
plot(t, y_C6H12, 'm')
plot(t_exp, y_exp_C6H12, 'mx')
plot(t, y_C7_plus, 'y')
plot(t_exp, y_exp_C7_plus, 'yx')
title('C4+ Olefin Mass Fractions with Methanol Feed at 643 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C4H8 - Mod', 'C4H8 - Exp', 'C5H10 - Mod', 'C5H10 - Exp',...
    'C6H12 - Mod', 'C6H12 - Exp', 'C7+ - Mod', 'C7+ - Exp')

figure
hold on
plot(t, y_C2H6, 'b')
plot(t_exp, y_exp_C2H6, 'bx')
plot(t, y_C3H8, 'r')
plot(t_exp, y_exp_C3H8, 'rx')
plot(t, y_C4H10, 'k')
plot(t_exp, y_exp_C4H10, 'kx')
plot(t, y_C5H12, 'c')
plot(t_exp, y_exp_C5H12, 'cx')
plot(t, y_C6H14, 'm')
plot(t_exp, y_exp_C6H14, 'mx')
title('C2-6 Paraffin Mass Fractions with Methanol Feed at 643 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C2H6 - Mod', 'C2H6 - Exp', 'C3H8 - Mod', 'C3H8 - Exp',...
    'C4H10 - Mod', 'C4H10 - Exp', 'C5H12 - Mod', 'C5H12 - Exp',...
    'C6H14 - Mod', 'C6H14 - Exp')

figure
hold on
plot(t, y_C6H6, 'k')
plot(t_exp, y_exp_C6H6, 'kx')
plot(t, y_C7H8, 'b')
plot(t_exp, y_exp_C7H8, 'bx')
plot(t, y_C8H10, 'r')
plot(t_exp, y_exp_C8H10, 'rx')
plot(t, y_C9H12, 'g')
plot(t_exp, y_exp_C9H12, 'gx')
plot(t, y_C10H14, 'm')
plot(t_exp, y_exp_C10H14, 'mx')
title('Aromatic Mass Fractions with Methanol Feed at 643 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C6 - Mod', 'C6 - Exp', 'C7 - Mod', 'C7 - Exp', 'C8 - Mod',...
    'C8 - Exp', 'C9 - Mod', 'C9 - Exp', 'C10 - Mod', 'C10 - Exp')

% Error Calculation
sol = ode15s(@Kinetic_Model, tspan, yi); % ODE system solution

t0_mod = deval(sol, t_exp(1))'; % Model solution at contact time 1
t1_mod = deval(sol, t_exp(2))'; % Model solution at contact time 2
t2_mod = deval(sol, t_exp(3))'; % Model solution at contact time 3
t3_mod = deval(sol, t_exp(4))'; % Model solution at contact time 4
t4_mod = deval(sol, t_exp(5))'; % Model solution at contact time 5
t5_mod = deval(sol, t_exp(6))'; % Model solution at contact time 6
t6_mod = deval(sol, t_exp(7))'; % Model solution at contact time 7
t7_mod = deval(sol, t_exp(8))'; % Model solution at contact time 8
t8_mod = deval(sol, t_exp(9))'; % Model solution at contact time 9

y_t0_mod = [t0_mod(1:17), sum((t0_mod(18:19)) + sum(t0_mod(21:22))...
    + sum(t0_mod(24:25)) + sum(t0_mod(27:30))), t0_mod(20), t0_mod(23),...
    t0_mod(26)]; % Model mass fractions at contact time 1
y_t1_mod = [t1_mod(1:17), sum((t1_mod(18:19)) + sum(t1_mod(21:22))...
    + sum(t1_mod(24:25)) + sum(t1_mod(27:30))), t1_mod(20), t1_mod(23),...
    t1_mod(26)]; % Model mass fractions at contact time 2
y_t2_mod = [t2_mod(1:17), sum((t2_mod(18:19)) + sum(t2_mod(21:22))...
    + sum(t2_mod(24:25)) + sum(t2_mod(27:30))), t2_mod(20), t2_mod(23),...
    t2_mod(26)]; % Model mass fractions at contact time 3
y_t3_mod = [t3_mod(1:17), sum((t3_mod(18:19)) + sum(t3_mod(21:22))...
    + sum(t3_mod(24:25)) + sum(t3_mod(27:30))), t3_mod(20), t3_mod(23),...
    t3_mod(26)]; % Model mass fractions at contact time 4
y_t4_mod = [t4_mod(1:17), sum((t4_mod(18:19)) + sum(t4_mod(21:22))...
    + sum(t4_mod(24:25)) + sum(t4_mod(27:30))), t4_mod(20), t4_mod(23),...
    t4_mod(26)]; % Model mass fractions at contact time 5
y_t5_mod = [t5_mod(1:17), sum((t5_mod(18:19)) + sum(t5_mod(21:22))...
    + sum(t5_mod(24:25)) + sum(t5_mod(27:30))), t5_mod(20), t5_mod(23),...
    t5_mod(26)]; % Model mass fractions at contact time 6
y_t6_mod = [t6_mod(1:17), sum((t6_mod(18:19)) + sum(t6_mod(21:22))...
    + sum(t6_mod(24:25)) + sum(t6_mod(27:30))), t6_mod(20), t6_mod(23),...
    t6_mod(26)]; % Model mass fractions at contact time 7
y_t7_mod = [t7_mod(1:17), sum((t7_mod(18:19)) + sum(t7_mod(21:22))...
    + sum(t7_mod(24:25)) + sum(t7_mod(27:30))), t7_mod(20), t7_mod(23),...
    t7_mod(26)]; % Model mass fractions at contact time 8
y_t8_mod = [t8_mod(1:17), sum((t8_mod(18:19)) + sum(t8_mod(21:22))...
    + sum(t8_mod(24:25)) + sum(t8_mod(27:30))), t8_mod(20), t8_mod(23),...
    t8_mod(26)]; % Model mass fractions at contact time 9

mat_mod = [y_t0_mod; y_t1_mod; y_t2_mod; y_t3_mod; y_t4_mod; y_t5_mod;...
    y_t6_mod; y_t7_mod; y_t8_mod]; % Modelled mass fraction values matrix
mat_exp = [y_exp_CH3OH', y_exp_CH3OCH3', y_exp_H2O', y_exp_CO',...
    y_exp_CH4', y_exp_C2H4', y_exp_C2H6', y_exp_C3H6', y_exp_C3H8',...
    y_exp_C4H8', y_exp_C4H10', y_exp_C5H10', y_exp_C5H12', y_exp_C6H6',...
    y_exp_C6H12', y_exp_C6H14', y_exp_C7H8', y_exp_C7_plus',...
    y_exp_C8H10', y_exp_C9H12', y_exp_C10H14']; % Experimental mass...
    % fraction values matrix

SSE = (mat_mod - mat_exp).^2; % Squared errors

WSSE_CH3OH = (SSE(:,1))/(sum(mat_exp(:,1))); % Weight for Methanol
WSSE_CH3OCH3 = (SSE(:,2))/(sum(mat_exp(:,2))); % Weight for DME
WSSE_H2O = (SSE(:,3))/(sum(mat_exp(:,3))); % Weight for water
WSSE_CO = (SSE(:,4))/(sum(mat_exp(:,4))); % Weight for carbon monoxide
WSSE_CH4 = (SSE(:,5))/(sum(mat_exp(:,5))); % Weight for methane
WSSE_C2H4 = (SSE(:,6))/(sum(mat_exp(:,6))); % Weight for ethene
WSSE_C2H6 = (SSE(:,7))/(sum(mat_exp(:,7))); % Weight for ethane
WSSE_C3H6 = (SSE(:,8))/(sum(mat_exp(:,8))); % Weight for propene
WSSE_C3H8 = (SSE(:,9))/(sum(mat_exp(:,9))); % Weight for propane
WSSE_C4H8 = (SSE(:,10))/(sum(mat_exp(:,10))); % Weight for butene
WSSE_C4H10 = (SSE(:,11))/(sum(mat_exp(:,11))); % Weight for butane
WSSE_C5H10 = (SSE(:,12))/(sum(mat_exp(:,12))); % Weight for pentene
WSSE_C5H12 = (SSE(:,13))/(sum(mat_exp(:,13))); % Weight for pentane
WSSE_C6H6 = (SSE(:,14))/(sum(mat_exp(:,14))); % Weight for benzene
WSSE_C6H12 = (SSE(:,15))/(sum(mat_exp(:,15))); % Weight for hexene
WSSE_C6H14 = (SSE(:,16))/(sum(mat_exp(:,16))); % Weight for hexane
WSSE_C7H8 = (SSE(:,17))/(sum(mat_exp(:,17))); % Weight for toluene
WSSE_C7_plus = (SSE(:,18))/(sum(mat_exp(:,18))); % Weight for C7+
WSSE_C8H10 = (SSE(:,19))/(sum(mat_exp(:,19))); % Weight for ethylbenzene
WSSE_C9H12 = (SSE(:,20))/(sum(mat_exp(:,20))); % Weight for propylbenzene
WSSE_C10H14 = (SSE(:,21))/(sum(mat_exp(:,21))); % Weight for butylbenzene

mat_WSSE = [WSSE_CH3OH, WSSE_CH3OCH3, WSSE_H2O, WSSE_CO, WSSE_CH4,...
    WSSE_C2H4, WSSE_C2H6, WSSE_C3H6, WSSE_C3H8, WSSE_C4H8, WSSE_C4H10,...
    WSSE_C5H10, WSSE_C5H12, WSSE_C6H6, WSSE_C6H12, WSSE_C6H14,...
    WSSE_C7H8, WSSE_C7_plus, WSSE_C8H10, WSSE_C9H12, WSSE_C10H14];...
    % Weighted sum of squared errors matrix

WSSE = sum(sum(mat_WSSE))*100; % Weighted sum of sqaured errors...
    % (multiplied by 100 to avoid rounding errors)

fprintf('WSSE = %.5f\n', WSSE) % Error outputted

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
    'Mod 643 K')
writetable(data_exp,'Methanol-to-Hydrocarbons_Data.xlsx', 'Sheet',...
    'Exp 643 K')