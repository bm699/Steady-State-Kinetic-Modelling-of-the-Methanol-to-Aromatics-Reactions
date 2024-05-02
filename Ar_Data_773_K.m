%% CE40254
% Methanol-to-aromatics kinetic model at 773 K to determine rate constants 
% for olefin aromatisation
% Written by Ben Mayoh
% Last modified: 15/04/24

clc
clear all
close all

%% Ethene

% Duration (s)
ti_C2 = 0; % Intial time
tf_C2 = 1000; % Final time
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

% Modelled Data for Ethene at 773 K
y_C2_CH4 = y_C2(:,5); % Methane
y_C2_C2H4 = y_C2(:,6); % Ethene
y_C2_C2H6 = y_C2(:,7); % Ethane
y_C2_C3H6 = y_C2(:,8); % Propene
y_C2_C3H8 = y_C2(:,9); % Propane
y_C2_C4H8 = y_C2(:,10); % Butene
y_C2_C4H10 = y_C2(:,11); % Butane
y_C2_C5H10 = y_C2(:,12); % Pentene
y_C2_C5H12 = y_C2(:,13); % Pentane
y_C2_C6H6 = y_C2(:,14); % Benzene
y_C2_C6H12 = y_C2(:,15); % Hexene
y_C2_C6H14 = y_C2(:,16); % Hexane
y_C2_C7H8 = y_C2(:,17); % Toluene
y_C2_C7_plus = y_C2(:,18) + y_C2(:,19) + y_C2(:,21) + y_C2(:,22)...
    + y_C2(:,24) + y_C2(:,25) + y_C2(:,27) + y_C2(:,28) + y_C2(:,29)...
    + y_C2(:,30); % C7+
y_C2_C8H10 = y_C2(:,20); % Ethylbenzene
y_C2_C9H12 = y_C2(:,23); % Propylbenzene
y_C2_C10H14 = y_C2(:,26); % Butylbenzene

% Experimental Data for Ethene at 773 K
t_exp_C2 = [1.742, 6.722, 16.688, 28.963, 47.505, 121.29, 395.35,...
    833.71]; % Contact time (s)
y_exp_C2_CH4 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.053, 0.268, 0.489]...
    /100; % Methane
y_exp_C2_C2H4 = [39.912, 39.398, 38.875, 35.047, 31.407, 13.093, 3.169,...
    2.256]/100; % Ethene
y_exp_C2_C2H6 = [0.000, 0.000, 0.000, 0.000, 0.047, 0.156, 0.241,...
    0.675]/100; % Ethane
y_exp_C2_C3H6 = [0.007, 0.108, 0.660, 1.856, 4.407, 7.924, 1.976, 0.976]...
    /100; % Propene
y_exp_C2_C3H8 = [0.000, 0.000, 0.000, 0.000, 0.111, 2.586, 8.695,...
    9.357]/100; % Propane
y_exp_C2_C4H8 = [0.162, 0.497, 1.135, 1.611, 2.765, 6.281, 1.901, 0.698]...
    /100; % Butene
y_exp_C2_C4H10 = [0.000, 0.000, 0.000, 0.000, 0.107, 1.771, 7.986,...
    8.045]/100; % Butane
y_exp_C2_C5H10 = [0.001, 0.061, 0.345, 0.747, 1.620, 2.010, 0.449,...
    0.190]/100; % Pentene
y_exp_C2_C5H12 = [0.000, 0.000, 0.000, 0.000, 0.098, 0.642, 2.203,...
    1.605]/100; % Pentane
y_exp_C2_C6H6 = [0.000, 0.000, 0.000, 0.000, 0.019, 0.390, 1.460, 2.245]...
    /100; % Benzene
y_exp_C2_C6H12 = [0.011, 0.049, 0.141, 0.232, 0.404, 0.611, 0.201,...
    0.182]/100; % Hexene
y_exp_C2_C6H14 = [0.000, 0.000, 0.000, 0.000, 0.095, 0.203, 1.208,...
    1.798]/100; % Hexane
y_exp_C2_C7H8 = [0.000, 0.000, 0.000, 0.000, 0.085, 1.789, 4.641, 6.697]...
    /100; % Toluene
y_exp_C2_C7_plus = [0.000, 0.000, 0.000, 0.000, 0.067, 0.305, 0.115,...
    0.076]/100; % C7+
y_exp_C2_C8H10 = [0.000, 0.000, 0.000, 0.000, 0.099, 1.513, 4.482,...
    5.785]/100; % Ethylbenzene
y_exp_C2_C9H12 = [0.000, 0.000, 0.000, 0.000, 0.042, 0.189, 0.774,...
    0.855]/100; % Propylbenzene
y_exp_C2_C10H14 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.032, 0.495,...
    0.368]/100; % Butylbenzene

% Plots
figure
hold on
plot(t_C2, y_C2_C6H6, 'k')
plot(t_exp_C2, y_exp_C2_C6H6, 'kx')
plot(t_C2, y_C2_C7H8, 'b')
plot(t_exp_C2, y_exp_C2_C7H8, 'bx')
plot(t_C2, y_C2_C8H10, 'r')
plot(t_exp_C2, y_exp_C2_C8H10, 'rx')
plot(t_C2, y_C2_C9H12, 'g')
plot(t_exp_C2, y_exp_C2_C9H12, 'gx')
plot(t_C2, y_C2_C10H14, 'm')
plot(t_exp_C2, y_exp_C2_C10H14, 'mx')
title('Aromatic Mass Fractions with Ethene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C6 - Mod', 'C6 - Exp', 'C7 - Mod', 'C7 - Exp', 'C8 - Mod',...
    'C8 - Exp', 'C9 - Mod', 'C9 - Exp', 'C10 - Mod', 'C10 - Exp')

figure
hold on
plot(t_C2, y_C2_CH4, 'g')
plot(t_exp_C2, y_exp_C2_CH4, 'gx')
plot(t_C2, y_C2_C2H6, 'b')
plot(t_exp_C2, y_exp_C2_C2H6, 'bx')
plot(t_C2, y_C2_C3H8, 'r')
plot(t_exp_C2, y_exp_C2_C3H8, 'rx')
plot(t_C2, y_C2_C4H10, 'k')
plot(t_exp_C2, y_exp_C2_C4H10, 'kx')
plot(t_C2, y_C2_C5H12, 'c')
plot(t_exp_C2, y_exp_C2_C5H12, 'cx')
plot(t_C2, y_C2_C6H14, 'y')
plot(t_exp_C2, y_exp_C2_C6H14, 'yx')
plot(t_C2, y_C2_C7_plus, '--m')
plot(t_exp_C2, y_exp_C2_C7_plus, 'm+')
title('Paraffin Mass Fractions with Ethene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('CH4 - Mod', 'CH4 - Exp', 'C2H6 - Mod', 'C2H6 - Exp',...
    'C3H8 - Mod', 'C3H8 - Exp', 'C4H10 - Mod', 'C4H10 - Exp',...
    'C5H12 - Mod', 'C5H12 - Exp', 'C6H14 - Mod', 'C6H14 - Exp',...
    'C7+ - Mod', 'C7+ - Exp')

figure
hold on
plot(t_C2, y_C2_C2H4, 'b')
plot(t_exp_C2, y_exp_C2_C2H4, 'bx')
plot(t_C2, y_C2_C3H6, 'r')
plot(t_exp_C2, y_exp_C2_C3H6, 'rx')
plot(t_C2, y_C2_C4H8, 'k')
plot(t_exp_C2, y_exp_C2_C4H8, 'kx')
plot(t_C2, y_C2_C5H10, 'c')
plot(t_exp_C2, y_exp_C2_C5H10, 'cx')
plot(t_C2, y_C2_C6H12, 'y')
plot(t_exp_C2, y_exp_C2_C6H12, 'yx')
plot(t_C2, y_C2_C7_plus, '--m')
plot(t_exp_C2, y_exp_C2_C7_plus, 'm+')
title('Olefin Mass Fractions with Ethene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C2H4 - Mod', 'C2H4 - Exp', 'C3H6 - Mod', 'C3H6 - Exp',...
    'C4H8 - Mod', 'C4H8 - Exp', 'C5H10 - Mod', 'C5H10 - Exp',...
    'C6H12 - Mod', 'C6H12 - Exp', 'C7+ - Mod', 'C7+ - Exp')

% Error Calculation
sol = ode15s(@Kinetic_Model, tspan_C2, yi_C2); % ODE system solution

t0_mod_C2 = deval(sol, t_exp_C2(1))'; % Model solution at contact time 1
t1_mod_C2 = deval(sol, t_exp_C2(2))'; % Model solution at contact time 2
t2_mod_C2 = deval(sol, t_exp_C2(3))'; % Model solution at contact time 3
t3_mod_C2 = deval(sol, t_exp_C2(4))'; % Model solution at contact time 4
t4_mod_C2 = deval(sol, t_exp_C2(5))'; % Model solution at contact time 5
t5_mod_C2 = deval(sol, t_exp_C2(6))'; % Model solution at contact time 6
t6_mod_C2 = deval(sol, t_exp_C2(7))'; % Model solution at contact time 7
t7_mod_C2 = deval(sol, t_exp_C2(8))'; % Model solution at contact time 8

y_t0_mod_C2 = [t0_mod_C2(5:17), sum((t0_mod_C2(18:19))...
    + sum(t0_mod_C2(21:22)) + sum(t0_mod_C2(24:25))...
    + sum(t0_mod_C2(27:30))), t0_mod_C2(20), t0_mod_C2(23),...
    t0_mod_C2(26)]; % Model mass fractions at contact time 1
y_t1_mod_C2 = [t1_mod_C2(5:17), sum((t1_mod_C2(18:19))...
    + sum(t1_mod_C2(21:22)) + sum(t1_mod_C2(24:25))...
    + sum(t1_mod_C2(27:30))), t1_mod_C2(20), t1_mod_C2(23),...
    t1_mod_C2(26)]; % Model mass fractions at contact time 2
y_t2_mod_C2 = [t2_mod_C2(5:17), sum((t2_mod_C2(18:19))...
    + sum(t2_mod_C2(21:22)) + sum(t2_mod_C2(24:25))...
    + sum(t2_mod_C2(27:30))), t2_mod_C2(20), t2_mod_C2(23),...
    t2_mod_C2(26)]; % Model mass fractions at contact time 3
y_t3_mod_C2 = [t3_mod_C2(5:17), sum((t3_mod_C2(18:19))...
    + sum(t3_mod_C2(21:22)) + sum(t3_mod_C2(24:25))...
    + sum(t3_mod_C2(27:30))), t3_mod_C2(20), t3_mod_C2(23),...
    t3_mod_C2(26)]; % Model mass fractions at contact time 4
y_t4_mod_C2 = [t4_mod_C2(5:17), sum((t4_mod_C2(18:19))...
    + sum(t4_mod_C2(21:22)) + sum(t4_mod_C2(24:25))...
    + sum(t4_mod_C2(27:30))), t4_mod_C2(20), t4_mod_C2(23),...
    t4_mod_C2(26)]; % Model mass fractions at contact time 5
y_t5_mod_C2 = [t5_mod_C2(5:17), sum((t5_mod_C2(18:19))...
    + sum(t5_mod_C2(21:22)) + sum(t5_mod_C2(24:25))...
    + sum(t5_mod_C2(27:30))), t5_mod_C2(20), t5_mod_C2(23),...
    t5_mod_C2(26)]; % Model mass fractions at contact time 6
y_t6_mod_C2 = [t6_mod_C2(5:17), sum((t6_mod_C2(18:19))...
    + sum(t6_mod_C2(21:22)) + sum(t6_mod_C2(24:25))...
    + sum(t6_mod_C2(27:30))), t6_mod_C2(20), t6_mod_C2(23),...
    t6_mod_C2(26)]; % Model mass fractions at contact time 7
y_t7_mod_C2 = [t7_mod_C2(5:17), sum((t7_mod_C2(18:19))...
    + sum(t7_mod_C2(21:22)) + sum(t7_mod_C2(24:25))...
    + sum(t7_mod_C2(27:30))), t7_mod_C2(20), t7_mod_C2(23),...
    t7_mod_C2(26)]; % Model mass fractions at contact time 8

mat_mod_C2 = [y_t0_mod_C2; y_t1_mod_C2; y_t2_mod_C2; y_t3_mod_C2;...
    y_t4_mod_C2; y_t5_mod_C2; y_t6_mod_C2; y_t7_mod_C2]; % Modelled mass...
    % fraction values matrix
mat_exp_C2 = [y_exp_C2_CH4', y_exp_C2_C2H4', y_exp_C2_C2H6',...
    y_exp_C2_C3H6', y_exp_C2_C3H8', y_exp_C2_C4H8', y_exp_C2_C4H10',...
    y_exp_C2_C5H10', y_exp_C2_C5H12', y_exp_C2_C6H6', y_exp_C2_C6H12',...
    y_exp_C2_C6H14', y_exp_C2_C7H8', y_exp_C2_C7_plus', y_exp_C2_C8H10'...
    y_exp_C2_C9H12', y_exp_C2_C10H14']; % Experimental mass fraction...
    % values matrix

SSE_C2 = (mat_mod_C2 - mat_exp_C2).^2; % Squared errors

WSSE_C2_CH4 = (SSE_C2(:,1))/(sum(mat_exp_C2(:,1))); % Weight for methane
WSSE_C2_C2H4 = (SSE_C2(:,2))/(sum(mat_exp_C2(:,2))); % Weight for ethene
WSSE_C2_C2H6 = (SSE_C2(:,3))/(sum(mat_exp_C2(:,3))); % Weight for ethane
WSSE_C2_C3H6 = (SSE_C2(:,4))/(sum(mat_exp_C2(:,4))); % Weight for propene
WSSE_C2_C3H8 = (SSE_C2(:,5))/(sum(mat_exp_C2(:,5))); % Weight for propane
WSSE_C2_C4H8 = (SSE_C2(:,6))/(sum(mat_exp_C2(:,6))); % Weight for butene
WSSE_C2_C4H10 = (SSE_C2(:,7))/(sum(mat_exp_C2(:,7))); % Weight for butane
WSSE_C2_C5H10 = (SSE_C2(:,8))/(sum(mat_exp_C2(:,8))); % Weight for pentene
WSSE_C2_C5H12 = (SSE_C2(:,9))/(sum(mat_exp_C2(:,9))); % Weight for pentane
WSSE_C2_C6H6 = (SSE_C2(:,10))/(sum(mat_exp_C2(:,10))); % Weight for benzene
WSSE_C2_C6H12 = (SSE_C2(:,11))/(sum(mat_exp_C2(:,11))); % Weight for hexene
WSSE_C2_C6H14 = (SSE_C2(:,12))/(sum(mat_exp_C2(:,12))); % Weight for hexane
WSSE_C2_C7H8 = (SSE_C2(:,13))/(sum(mat_exp_C2(:,13))); % Weight for toluene
WSSE_C2_C7_plus = (SSE_C2(:,14))/(sum(mat_exp_C2(:,14))); % Weight for C7+
WSSE_C2_C8H10 = (SSE_C2(:,15))/(sum(mat_exp_C2(:,15))); % Weight for...
    % ethylbenzene
WSSE_C2_C9H12 = (SSE_C2(:,16))/(sum(mat_exp_C2(:,16))); % Weight for...
    % propylbenzene
WSSE_C2_C10H14 = (SSE_C2(:,17))/(sum(mat_exp_C2(:,17))); % Weight for...
    % butylbenzene

mat_WSSE_C2 = [WSSE_C2_CH4, WSSE_C2_C2H4, WSSE_C2_C2H6, WSSE_C2_C3H6,...
    WSSE_C2_C3H8, WSSE_C2_C4H8, WSSE_C2_C4H10, WSSE_C2_C5H10,...
    WSSE_C2_C5H12, WSSE_C2_C6H6, WSSE_C2_C6H12, WSSE_C2_C6H14,...
    WSSE_C2_C7H8, WSSE_C2_C7_plus, WSSE_C2_C8H10, WSSE_C2_C9H12,...
    WSSE_C2_C10H14]; % Weighted sum of squared errors matrix

WSSE_C2 = sum(sum(mat_WSSE_C2))*100; % Weighted sum of sqaured errors...
    % (multiplied by 100 to avoid rounding errors)

fprintf('Ethene WSSE = %.3f\n', WSSE_C2) % Error outputted

% Writing to Excel File
data_C2_mod = table(t_C2, y_C2_CH4, y_C2_C2H4, y_C2_C2H6, y_C2_C3H6,...
    y_C2_C3H8, y_C2_C4H8, y_C2_C4H10, y_C2_C5H10, y_C2_C5H12,...
    y_C2_C6H6, y_C2_C6H12, y_C2_C6H14, y_C2_C7H8, y_C2_C7_plus,...
    y_C2_C8H10, y_C2_C9H12, y_C2_C10H14);
data_C2_exp = table(t_exp_C2', y_exp_C2_CH4', y_exp_C2_C2H4',...
    y_exp_C2_C2H6', y_exp_C2_C3H6', y_exp_C2_C3H8', y_exp_C2_C4H8',...
    y_exp_C2_C4H10', y_exp_C2_C5H10', y_exp_C2_C5H12', y_exp_C2_C6H6',...
    y_exp_C2_C6H12', y_exp_C2_C6H14', y_exp_C2_C7H8', y_exp_C2_C7_plus',...
    y_exp_C2_C8H10', y_exp_C2_C9H12', y_exp_C2_C10H14');

writetable(data_C2_mod, 'Aromatisation_Data.xlsx', 'Sheet', 'C2 Mod 773 K')
writetable(data_C2_exp,'Aromatisation_Data.xlsx', 'Sheet', 'C2 Exp 773 K')

%% Propene

% Duration (s)
ti_C3 = 0; % Intial time
tf_C3 = 220; % Final time
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
y_C3_CH4 = y_C3(:,5); % Methane
y_C3_C2H4 = y_C3(:,6); % Ethene
y_C3_C2H6 = y_C3(:,7); % Ethane
y_C3_C3H6 = y_C3(:,8); % Propene
y_C3_C3H8 = y_C3(:,9); % Propane
y_C3_C4H8 = y_C3(:,10); % Butene
y_C3_C4H10 = y_C3(:,11); % Butane
y_C3_C5H10 = y_C3(:,12); % Pentene
y_C3_C5H12 = y_C3(:,13); % Pentane
y_C3_C6H6 = y_C3(:,14); % Benzene
y_C3_C6H12 = y_C3(:,15); % Hexene
y_C3_C6H14 = y_C3(:,16); % Hexane
y_C3_C7H8 = y_C3(:,17); % Toluene
y_C3_C7_plus = y_C3(:,18) + y_C3(:,19) + y_C3(:,21) + y_C3(:,22)...
    + y_C3(:,24) + y_C3(:,25) + y_C3(:,27) + y_C3(:,28) + y_C3(:,29)...
    + y_C3(:,30); % C7+
y_C3_C8H10 = y_C3(:,20); % Ethylbenzene
y_C3_C9H12 = y_C3(:,23); % Propylbenzene
y_C3_C10H14 = y_C3(:,26); % Butylbenzene

% Experimental Data for Propene at 773 K
t_exp_C3 = [0.273, 0.469, 0.934, 2.039, 4.973, 7.620, 74.950, 183.160];...
    % Contact time (s)
y_exp_C3_CH4 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.207, 0.227]...
    /100; % Methane
y_exp_C3_C2H4 = [0.025, 0.106, 0.423, 0.519, 1.718, 2.103, 3.881, 2.894]...
    /100; % Ethene
y_exp_C3_C2H6 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.011, 0.270, 0.284]...
    /100; % Ethane
y_exp_C3_C3H6 = [46.237, 43.710, 40.040, 34.424, 27.009, 23.559, 11.921,...
    6.422]/100; % Propene
y_exp_C3_C3H8 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.022, 5.515, 7.146]...
    /100; % Propane
y_exp_C3_C4H8 = [0.192, 0.766, 1.847, 5.910, 12.211, 15.187, 5.352,...
    2.729]/100; % Butene
y_exp_C3_C4H10 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.223, 6.030,...
    8.770]/100; % Butane
y_exp_C3_C5H10 = [0.165, 0.421, 0.882, 2.615, 4.743, 4.355, 2.008,...
    1.007]/100; % Pentene
y_exp_C3_C5H12 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.097, 1.800,...
    2.180]/100; % Pentane
y_exp_C3_C6H6 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.038, 1.039, 1.462]...
    /100; % Benzene
y_exp_C3_C6H12 = [3.374, 4.928, 6.549, 6.140, 3.834, 3.133, 1.315,...
    0.610]/100; % Hexene
y_exp_C3_C6H14 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.092, 1.120,...
    2.811]/100; % Hexane
y_exp_C3_C7H8 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.061, 3.254, 4.446]...
    /100; % Toluene
y_exp_C3_C7_plus = [0.043, 0.155, 0.287, 0.491, 0.590, 0.720, 0.430,...
    0.217]/100; % C7+
y_exp_C3_C8H10 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.118, 4.198,...
    6.151]/100; % Ethylbenzene
y_exp_C3_C9H12 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.017, 0.546,...
    0.974]/100; % Propylbenzene
y_exp_C3_C10H14 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.105,...
    0.195]/100; % Butylbenzene

% Plots
figure
hold on
plot(t_C3, y_C3_C6H6, 'k')
plot(t_exp_C3, y_exp_C3_C6H6, 'kx')
plot(t_C3, y_C3_C7H8, 'b')
plot(t_exp_C3, y_exp_C3_C7H8, 'bx')
plot(t_C3, y_C3_C8H10, 'r')
plot(t_exp_C3, y_exp_C3_C8H10, 'rx')
plot(t_C3, y_C3_C9H12, 'g')
plot(t_exp_C3, y_exp_C3_C9H12, 'gx')
plot(t_C3, y_C3_C10H14, 'm')
plot(t_exp_C3, y_exp_C3_C10H14, 'mx')
title('Aromatic Mass Fractions with Propene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C6 - Mod', 'C6 - Exp', 'C7 - Mod', 'C7 - Exp', 'C8 - Mod',...
    'C8 - Exp', 'C9 - Mod', 'C9 - Exp', 'C10 - Mod', 'C10 - Exp')

figure
hold on
plot(t_C3, y_C3_CH4, 'g')
plot(t_exp_C3, y_exp_C3_CH4, 'gx')
plot(t_C3, y_C3_C2H6, 'b')
plot(t_exp_C3, y_exp_C3_C2H6, 'bx')
plot(t_C3, y_C3_C3H8, 'r')
plot(t_exp_C3, y_exp_C3_C3H8, 'rx')
plot(t_C3, y_C3_C4H10, 'k')
plot(t_exp_C3, y_exp_C3_C4H10, 'kx')
plot(t_C3, y_C3_C5H12, 'c')
plot(t_exp_C3, y_exp_C3_C5H12, 'cx')
plot(t_C3, y_C3_C6H14, 'y')
plot(t_exp_C3, y_exp_C3_C6H14, 'yx')
plot(t_C3, y_C3_C7_plus, '--m')
plot(t_exp_C3, y_exp_C3_C7_plus, 'm+')
title('Paraffin Mass Fractions with Propene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('CH4 - Mod', 'CH4 - Exp', 'C2H6 - Mod', 'C2H6 - Exp',...
    'C3H8 - Mod', 'C3H8 - Exp', 'C4H10 - Mod', 'C4H10 - Exp',...
    'C5H12 - Mod', 'C5H12 - Exp', 'C6H14 - Mod', 'C6H14 - Exp',...
    'C7+ - Mod', 'C7+ - Exp')

figure
hold on
plot(t_C3, y_C3_C2H4, 'b')
plot(t_exp_C3, y_exp_C3_C2H4, 'bx')
plot(t_C3, y_C3_C3H6, 'r')
plot(t_exp_C3, y_exp_C3_C3H6, 'rx')
plot(t_C3, y_C3_C4H8, 'k')
plot(t_exp_C3, y_exp_C3_C4H8, 'kx')
plot(t_C3, y_C3_C5H10, 'c')
plot(t_exp_C3, y_exp_C3_C5H10, 'cx')
plot(t_C3, y_C3_C6H12, 'y')
plot(t_exp_C3, y_exp_C3_C6H12, 'yx')
plot(t_C3, y_C3_C7_plus, '--m')
plot(t_exp_C3, y_exp_C3_C7_plus, 'm+')
title('Olefin Mass Fractions with Propene Feed at 773 K')
xlabel('Contact Time (s)')
ylabel('Mass Fraction')
legend('C2H4 - Mod', 'C2H4 - Exp', 'C3H6 - Mod', 'C3H6 - Exp',...
    'C4H8 - Mod', 'C4H8 - Exp', 'C5H10 - Mod', 'C5H10 - Exp',...
    'C6H12 - Mod', 'C6H12 - Exp', 'C7+ - Mod', 'C7+ - Exp')

% Error Calculation
sol = ode15s(@Kinetic_Model, tspan_C3, yi_C3); % ODE system solution

t0_mod_C3 = deval(sol, t_exp_C3(1))'; % Model solution at contact time 1
t1_mod_C3 = deval(sol, t_exp_C3(2))'; % Model solution at contact time 2
t2_mod_C3 = deval(sol, t_exp_C3(3))'; % Model solution at contact time 3
t3_mod_C3 = deval(sol, t_exp_C3(4))'; % Model solution at contact time 4
t4_mod_C3 = deval(sol, t_exp_C3(5))'; % Model solution at contact time 5
t5_mod_C3 = deval(sol, t_exp_C3(6))'; % Model solution at contact time 6
t6_mod_C3 = deval(sol, t_exp_C3(7))'; % Model solution at contact time 7
t7_mod_C3 = deval(sol, t_exp_C3(8))'; % Model solution at contact time 8

y_t0_mod_C3 = [t0_mod_C3(5:17), sum((t0_mod_C3(18:19))...
    + sum(t0_mod_C3(21:22)) + sum(t0_mod_C3(24:25))...
    + sum(t0_mod_C3(27:30))), t0_mod_C3(20), t0_mod_C3(23),...
    t0_mod_C3(26)]; % Model mass fractions at contact time 1
y_t1_mod_C3 = [t1_mod_C3(5:17), sum((t1_mod_C3(18:19))...
    + sum(t1_mod_C3(21:22)) + sum(t1_mod_C3(24:25))...
    + sum(t1_mod_C3(27:30))), t1_mod_C3(20), t1_mod_C3(23),...
    t1_mod_C3(26)]; % Model mass fractions at contact time 2
y_t2_mod_C3 = [t2_mod_C3(5:17), sum((t2_mod_C3(18:19))...
    + sum(t2_mod_C3(21:22)) + sum(t2_mod_C3(24:25))...
    + sum(t2_mod_C3(27:30))), t2_mod_C3(20), t2_mod_C3(23),...
    t2_mod_C3(26)]; % Model mass fractions at contact time 3
y_t3_mod_C3 = [t3_mod_C3(5:17), sum((t3_mod_C3(18:19))...
    + sum(t3_mod_C3(21:22)) + sum(t3_mod_C3(24:25))...
    + sum(t3_mod_C3(27:30))), t3_mod_C3(20), t3_mod_C3(23),...
    t3_mod_C3(26)]; % Model mass fractions at contact time 4
y_t4_mod_C3 = [t4_mod_C3(5:17), sum((t4_mod_C3(18:19))...
    + sum(t4_mod_C3(21:22)) + sum(t4_mod_C3(24:25))...
    + sum(t4_mod_C3(27:30))), t4_mod_C3(20), t4_mod_C3(23),...
    t4_mod_C3(26)]; % Model mass fractions at contact time 5
y_t5_mod_C3 = [t5_mod_C3(5:17), sum((t5_mod_C3(18:19))...
    + sum(t5_mod_C3(21:22)) + sum(t5_mod_C3(24:25))...
    + sum(t5_mod_C3(27:30))), t5_mod_C3(20), t5_mod_C3(23),...
    t5_mod_C3(26)]; % Model mass fractions at contact time 6
y_t6_mod_C3 = [t6_mod_C3(5:17), sum((t6_mod_C3(18:19))...
    + sum(t6_mod_C3(21:22)) + sum(t6_mod_C3(24:25))...
    + sum(t6_mod_C3(27:30))), t6_mod_C3(20), t6_mod_C3(23),...
    t6_mod_C3(26)]; % Model mass fractions at contact time 7
y_t7_mod_C3 = [t7_mod_C3(5:17), sum((t7_mod_C3(18:19))...
    + sum(t7_mod_C3(21:22)) + sum(t7_mod_C3(24:25))...
    + sum(t7_mod_C3(27:30))), t7_mod_C3(20), t7_mod_C3(23),...
    t7_mod_C3(26)]; % Model mass fractions at contact time 8

mat_mod_C3 = [y_t0_mod_C3; y_t1_mod_C3; y_t2_mod_C3; y_t3_mod_C3;...
    y_t4_mod_C3; y_t5_mod_C3; y_t6_mod_C3; y_t7_mod_C3]; % Modelled mass...
    % fraction values matrix
mat_exp_C3 = [y_exp_C3_CH4', y_exp_C3_C2H4', y_exp_C3_C2H6',...
    y_exp_C3_C3H6', y_exp_C3_C3H8', y_exp_C3_C4H8', y_exp_C3_C4H10',...
    y_exp_C3_C5H10', y_exp_C3_C5H12', y_exp_C3_C6H6', y_exp_C3_C6H12',...
    y_exp_C3_C6H14', y_exp_C3_C7H8', y_exp_C3_C7_plus', y_exp_C3_C8H10'...
    y_exp_C3_C9H12', y_exp_C3_C10H14']; % Experimental mass fraction...
    % values matrix

SSE_C3 = (mat_mod_C3 - mat_exp_C3).^2; % Squared errors

WSSE_C3_CH4 = (SSE_C3(:,1))/(sum(mat_exp_C3(:,1))); % Weight for methane
WSSE_C3_C2H4 = (SSE_C3(:,2))/(sum(mat_exp_C3(:,2))); % Weight for ethene
WSSE_C3_C2H6 = (SSE_C3(:,3))/(sum(mat_exp_C3(:,3))); % Weight for ethane
WSSE_C3_C3H6 = (SSE_C3(:,4))/(sum(mat_exp_C3(:,4))); % Weight for propene
WSSE_C3_C3H8 = (SSE_C3(:,5))/(sum(mat_exp_C3(:,5))); % Weight for propane
WSSE_C3_C4H8 = (SSE_C3(:,6))/(sum(mat_exp_C3(:,6))); % Weight for butene
WSSE_C3_C4H10 = (SSE_C3(:,7))/(sum(mat_exp_C3(:,7))); % Weight for butane
WSSE_C3_C5H10 = (SSE_C3(:,8))/(sum(mat_exp_C3(:,8))); % Weight for pentene
WSSE_C3_C5H12 = (SSE_C3(:,9))/(sum(mat_exp_C3(:,9))); % Weight for pentane
WSSE_C3_C6H6 = (SSE_C3(:,10))/(sum(mat_exp_C3(:,10))); % Weight for benzene
WSSE_C3_C6H12 = (SSE_C3(:,11))/(sum(mat_exp_C3(:,11))); % Weight for hexene
WSSE_C3_C6H14 = (SSE_C3(:,12))/(sum(mat_exp_C3(:,12))); % Weight for hexane
WSSE_C3_C7H8 = (SSE_C3(:,13))/(sum(mat_exp_C3(:,13))); % Weight for toluene
WSSE_C3_C7_plus = (SSE_C3(:,14))/(sum(mat_exp_C3(:,14))); % Weight for C7+
WSSE_C3_C8H10 = (SSE_C3(:,15))/(sum(mat_exp_C3(:,15))); % Weight for...
    % ethylbenzene
WSSE_C3_C9H12 = (SSE_C3(:,16))/(sum(mat_exp_C3(:,16))); % Weight for...
    % propylbenzene
WSSE_C3_C10H14 = (SSE_C3(:,17))/(sum(mat_exp_C3(:,17))); % Weight for...
    % butylbenzene

mat_WSSE_C3 = [WSSE_C3_CH4, WSSE_C3_C2H4, WSSE_C3_C2H6, WSSE_C3_C3H6,...
    WSSE_C3_C3H8, WSSE_C3_C4H8, WSSE_C3_C4H10, WSSE_C3_C5H10,...
    WSSE_C3_C5H12, WSSE_C3_C6H6, WSSE_C3_C6H12, WSSE_C3_C6H14,...
    WSSE_C3_C7H8, WSSE_C3_C7_plus, WSSE_C3_C8H10, WSSE_C3_C9H12,...
    WSSE_C3_C10H14]; % Weighted sum of squared errors matrix

WSSE_C3 = sum(sum(mat_WSSE_C3))*100; % Weighted sum of sqaured errors...
    % (multiplied by 100 to avoid rounding errors)

fprintf('Propene WSSE = %.3f\n', WSSE_C3) % Error outputted

% Writing to Excel File
data_C3_mod = table(t_C3, y_C3_CH4, y_C3_C2H4, y_C3_C2H6, y_C3_C3H6,...
    y_C3_C3H8, y_C3_C4H8, y_C3_C4H10, y_C3_C5H10, y_C3_C5H12,...
    y_C3_C6H6, y_C3_C6H12, y_C3_C6H14, y_C3_C7H8, y_C3_C7_plus,...
    y_C3_C8H10, y_C3_C9H12, y_C3_C10H14);
data_C3_exp = table(t_exp_C3', y_exp_C3_CH4', y_exp_C3_C2H4',...
    y_exp_C3_C2H6', y_exp_C3_C3H6', y_exp_C3_C3H8', y_exp_C3_C4H8',...
    y_exp_C3_C4H10', y_exp_C3_C5H10', y_exp_C3_C5H12', y_exp_C3_C6H6',...
    y_exp_C3_C6H12', y_exp_C3_C6H14', y_exp_C3_C7H8', y_exp_C3_C7_plus',...
    y_exp_C3_C8H10', y_exp_C3_C9H12', y_exp_C3_C10H14');

writetable(data_C3_mod, 'Aromatisation_Data.xlsx', 'Sheet', 'C3 Mod 773 K')
writetable(data_C3_exp,'Aromatisation_Data.xlsx', 'Sheet', 'C3 Exp 773 K')