%% CE40254
% Methanol-to-aromatics kinetic model at 673 K to determine rate constants 
% for olefin aromatisation
% Written by Ben Mayoh
% Last modified: 15/04/24

clc
clear all
close all

%% Ethene

% Duration (s)
ti_C2 = 0; % Intial time
tf_C2 = 529; % Final time
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

% Experimental Data for Ethene at 673 K
t_exp_C2 = [8.03, 11.92, 20.64, 31.12, 67.75, 212.00, 441.00]; % Contact...
    % time (s)
y_exp_C2_CH4 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.020, 0.049]/100;...
    % Methane
y_exp_C2_C2H4 = [39.341, 38.441, 38.195, 36.620, 26.424, 1.515, 0.809]...
    /100; % Ethene
y_exp_C2_C2H6 = [0.000, 0.000, 0.000, 0.000, 0.017, 0.153, 0.248]/100;...
    % Ethane
y_exp_C2_C3H6 = [0.054, 0.170, 0.357, 0.814, 3.812, 1.503, 0.947]/100;...
    % Propene
y_exp_C2_C3H8 = [0.000, 0.000, 0.000, 0.000, 0.144, 4.491, 5.577]/100;...
    % Propane
y_exp_C2_C4H8 = [0.275, 0.341, 0.672, 0.940, 4.554, 2.593, 1.410]/100;...
    % Butene
y_exp_C2_C4H10 = [0.000, 0.000, 0.000, 0.000, 0.195, 9.298, 9.816]/100;...
    % Butane
y_exp_C2_C5H10 = [0.040, 0.089, 0.273, 0.665, 2.721, 1.610, 0.820]/100;...
    % Pentene
y_exp_C2_C5H12 = [0.000, 0.000, 0.000, 0.000, 0.125, 3.521, 3.823]/100;...
    % Pentane
y_exp_C2_C6H6 = [0.000, 0.000, 0.000, 0.000, 0.005, 0.588, 0.954]/100;...
    % Benzene
y_exp_C2_C6H12 = [0.112, 0.096, 0.232, 0.218, 0.906, 0.494, 0.203]/100;...
    % Hexene
y_exp_C2_C6H14 = [0.000, 0.000, 0.000, 0.000, 0.082, 1.903, 2.115]/100;...
    % Hexane
y_exp_C2_C7H8 = [0.000, 0.000, 0.000, 0.000, 0.047, 2.917, 4.771]/100;...
    % Toluene
y_exp_C2_C7_plus = [0.000, 0.000, 0.000, 0.000, 0.198, 0.868, 0.568]...
    /100; % C7+
y_exp_C2_C8H10 = [0.000, 0.000, 0.000, 0.000, 0.204, 5.201, 6.705]/100;...
    % Ethylbenzene
y_exp_C2_C9H12 = [0.000, 0.000, 0.000, 0.000, 0.038, 2.288, 2.435]/100;...
    % Propylbenzene
y_exp_C2_C10H14 = [0.000, 0.000, 0.000, 0.000, 0.000, 0.390, 0.969]/100;...
    % Butylbenzene

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
title('Aromatic Mass Fractions with Ethene Feed at 673 K')
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
title('Paraffin Mass Fractions with Ethene Feed at 673 K')
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
title('Olefin Mass Fractions with Ethene Feed at 673 K')
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

mat_mod_C2 = [y_t0_mod_C2; y_t1_mod_C2; y_t2_mod_C2; y_t3_mod_C2;...
    y_t4_mod_C2; y_t5_mod_C2; y_t6_mod_C2]; % Modelled mass fraction...
    % values matrix
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

writetable(data_C2_mod, 'Aromatisation_Data.xlsx', 'Sheet', 'C2 Mod 673 K')
writetable(data_C2_exp,'Aromatisation_Data.xlsx', 'Sheet', 'C2 Exp 673 K')