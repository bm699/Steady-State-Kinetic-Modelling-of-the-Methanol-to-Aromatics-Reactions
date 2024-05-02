%% CE40254
% Methanol-to-aromatics kinetic model function
% Written by Ben Mayoh
% Last modified: 26/03/24

function dmdt = Kinetic_Model(~, m)

% Pressure (bar)
P = 1;

% Molecular Weights (g/mol)
MW_CH3OH = 32; % Methanol
MW_CH3OCH3 = 46; % DME
MW_H2O = 18; % Water
MW_CO = 28; % Carbon Monoxide
MW_CH4 = 16; % Methane
MW_C2H4 = 28; % Ethene
MW_C2H6 = 30; % Ethane
MW_C3H6 = 42; % Propene
MW_C3H8 = 44; % Propane
MW_C4H8 = 56; % Butene
MW_C4H10 = 58; % Butane
MW_C5H10 = 70; % Pentene
MW_C5H12 = 72; % Pentane
MW_C6H6 = 78; % Benzene
MW_C6H12 = 84; % Hexene
MW_C6H14 = 86; % Hexane
MW_C7H8 = 92; % Toluene
MW_C7H14 = 98; % Heptene
MW_C7H16 = 100; % Heptane
MW_C8H10 = 104; % Ethylbenzene
MW_C8H16 = 112; % Octene
MW_C8H18 = 114; % Octane
MW_C9H12 = 117; % Propylbenzene
MW_C9H18 = 126; % Nonene
MW_C9H20 = 128; % Nonane
MW_C10H14 = 130; % Butylbenzene
MW_C10H20 = 140; % Decene
MW_C10H22 = 142; % Decane
MW_C11H22 = 154; % Undecene
MW_C12H24 = 168; % Dodecene
MW_H2 = 2; % Hydrogen
MW_N2 = 28; % Nitrogen

% Adsorption Equilibrium Constants (1/bar)
Ka1 = 10.44; % Methanol
Ka2 = 10.44; % DME
Ka3 = 4.30; % Water
Ka4 = 0.26; % Carbon Monoxide
Ka5 = 0.26; % Paraffins (C1-C10)
Ka6 = 4.30; % Olefins (C2-C12)
Ka7 = 70.63; % Aromatics (C6-C10)
Ka8 = 0; % Hydrogen

% Total Number of Moles (mol)
n_total = m(1)/MW_CH3OH + m(2)/MW_CH3OCH3 + m(3)/MW_H2O + m(4)/MW_CO...
    + m(5)/MW_CH4 + m(6)/MW_C2H4 + m(7)/MW_C2H6 + m(8)/MW_C3H6...
    + m(9)/MW_C3H8 + m(10)/MW_C4H8 + m(11)/MW_C4H10 + m(12)/MW_C5H10...
    + m(13)/MW_C5H12 + m(14)/MW_C6H6 + m(15)/MW_C6H12...
    + m(16)/MW_C6H14 + m(17)/MW_C7H8 + m(18)/MW_C7H14...
    + m(19)/MW_C7H16 + m(20)/MW_C8H10 + m(21)/MW_C8H16...
    + m(22)/MW_C8H18 + m(23)/MW_C9H12 + m(24)/MW_C9H18...
    + m(25)/MW_C9H20 + m(26)/MW_C10H14 + m(27)/MW_C10H20...
    + m(28)/MW_C10H22 + m(29)/MW_C11H22 + m(30)/MW_C12H24...
    + m(31)/MW_H2 + m(32)/MW_N2;

% Mole Fractions
y_CH3OH = (m(1)/MW_CH3OH)/n_total; % Methanol
y_CH3OCH3 = (m(2)/MW_CH3OCH3)/n_total; % DME
y_H2O = (m(3)/MW_H2O)/n_total; % Water
y_CO = (m(4)/MW_CO)/n_total; % Carbon Monoxide
y_CH4 = (m(5)/MW_CH4)/n_total; % Methane
y_C2H4 = (m(6)/MW_C2H4)/n_total; % Ethene
y_C2H6 = (m(7)/MW_C2H6)/n_total; % Ethane
y_C3H6 = (m(8)/MW_C3H6)/n_total; % Propene
y_C3H8 = (m(9)/MW_C3H8)/n_total; % Propane
y_C4H8 = (m(10)/MW_C4H8)/n_total; % Butene
y_C4H10 = (m(11)/MW_C4H10)/n_total; % Butane
y_C5H10 = (m(12)/MW_C5H10)/n_total; % Pentene
y_C5H12 = (m(13)/MW_C5H12)/n_total; % Pentane
y_C6H6 = (m(14)/MW_C6H6)/n_total; % Benzene
y_C6H12 = (m(15)/MW_C6H12)/n_total; % Hexene
y_C6H14 = (m(16)/MW_C6H14)/n_total; % Hexane
y_C7H8 = (m(17)/MW_C7H8)/n_total; % Toluene
y_C7H14 = (m(18)/MW_C7H14)/n_total; % Heptene
y_C7H16 = (m(19)/MW_C7H16)/n_total; % Heptane
y_C8H10 = (m(20)/MW_C8H10)/n_total; % Ethylbenzene
y_C8H16 = (m(21)/MW_C8H16)/n_total; % Octene
y_C8H18 = (m(22)/MW_C8H18)/n_total; % Octane
y_C9H12 = (m(23)/MW_C9H12)/n_total; % Propylbenzene
y_C9H18 = (m(24)/MW_C9H18)/n_total; % Nonene
y_C9H20 = (m(25)/MW_C9H20)/n_total; % Nonane
y_C10H14 = (m(26)/MW_C10H14)/n_total; % Butylbenzene
y_C10H20 = (m(27)/MW_C10H20)/n_total; % Decene
y_C10H22 = (m(28)/MW_C10H22)/n_total; % Decane
y_C11H22 = (m(29)/MW_C11H22)/n_total; % Undecene
y_C12H24 = (m(30)/MW_C12H24)/n_total; % Dodecene
y_H2 = (m(31)/MW_H2)/n_total; % Hydrogen
y_N2 = (m(32)/MW_N2)/n_total; % Nitrogen

% Partial Pressures (bar)
P_CH3OH = y_CH3OH*P; % Methanol
P_CH3OCH3 = y_CH3OCH3*P; % DME
P_H2O = y_H2O*P; % Water
P_CO = y_CO*P; % Carbon Monoxide
P_CH4 = y_CH4*P; % Methane
P_C2H4 = y_C2H4*P; % Ethene
P_C2H6 = y_C2H6*P; % Ethane
P_C3H6 = y_C3H6*P; % Propene
P_C3H8 = y_C3H8*P; % Propane
P_C4H8 = y_C4H8*P; % Butene
P_C4H10 = y_C4H10*P; % Butane
P_C5H10 = y_C5H10*P; % Pentene
P_C5H12 = y_C5H12*P; % Pentane
P_C6H6 = y_C6H6*P; % Benzene
P_C6H12 = y_C6H12*P; % Hexene
P_C6H14 = y_C6H14*P; % Hexane
P_C7H8 = y_C7H8*P; % Toluene
P_C7H14 = y_C7H14*P; % Heptene
P_C7H16 = y_C7H16*P; % Heptane
P_C8H10 = y_C8H10*P; % Ethylbenzene
P_C8H16 = y_C8H16*P; % Octene
P_C8H18 = y_C8H18*P; % Octane
P_C9H12 = y_C9H12*P; % Propylbenzene
P_C9H18 = y_C9H18*P; % Nonene
P_C9H20 = y_C9H20*P; % Nonane
P_C10H14 = y_C10H14*P; % Butylbenzene
P_C10H20 = y_C10H20*P; % Decene
P_C10H22 = y_C10H22*P; % Decane
P_C11H22 = y_C11H22*P; % Undecene
P_C12H24 = y_C12H24*P; % Dodecene
P_H2 = y_H2*P; % Hydrogen
P_N2 = y_N2*P; % Nitrogen

% Rate Constants (mol/g/s/bar) or (mol/g/s) (for olefin cracking)
k_df = 3.43e-2; % Methanol dehydration
k_db = 1.06e-1; % DME hydration

k_CO = 2.05e-3; % Carbon monoxide formation

k_c = 3.79e-1; % Surface acetyl group formation

k_ef = 2.00e-2; % Ethene formation

k_pf = 1.00e-2; % Propene formation

k_m2o = 4.81e-2; % Ethene methylation
k_m3o = 7.48e-2; % Propene methylation
k_m4o = 8.87e-2; % Butene methylation
k_m5o = 1.01e-1; % Pentene methylation
k_m6o = k_m5o; % Hexene methylation
k_m7o = k_m5o; % Heptene methylation
k_m8o = k_m5o; % Octene methylation
k_m9o = k_m5o; % Nonene methylation
k_m10o = k_m5o; % Decene methylation
k_m11o = k_m5o; % Undecene methylation

k_o1 = 2.24e-5; % Ethene + ethene oligomerisation
k_c1 = 2.05e-6; % Ethene + ethene cracking
k_o2 = 3.77e-3; % Ethene + propene oligomerisation
k_c2 = 7.81e-5; % Ethene + propene cracking
k_o3 = 3.14e-3; % Ethene + butene oligomerisation
k_c3 = 1.36e-4; % Ethene + butene cracking
k_o4 = k_o3; % Ethene + pentene oligomerisation
k_c4 = k_c3; % Ethene + pentene cracking
k_o5 = k_o3; % Ethene + hexene oligomerisation
k_c5 = k_c3; % Ethene + hexene cracking
k_o6 = k_o3; % Ethene + heptene oligomerisation
k_c6 = k_c3; % Ethene + heptene cracking
k_o7 = k_o3; % Ethene + octene oligomerisation
k_c7 = k_c3; % Ethene + octene cracking
k_o8 = k_o3; % Ethene + nonene oligomerisation
k_c8 = k_c3; % Ethene + nonene cracking
k_o9 = k_o3; % Ethene + decene oligomerisation
k_c9 = k_c3; % Ethene + decene cracking
k_o10 = 3.20e-2; % Propene + propene oligomerisation
k_c10 = 4.94e-3; % Propene + propene cracking
k_o11 = 3.50e-2; % Propene + butene oligomerisation
k_c11 = 3.11e-2; % Propene + butene cracking
k_o12 = k_o11; % Propene + pentene oligomerisation
k_c12 = k_c11; % Propene + pentene cracking
k_o13 = k_o11; % Propene + hexene oligomerisation
k_c13 = k_c11; % Propene + hexene cracking
k_o14 = k_o11; % Propene + heptene oligomerisation
k_c14 = k_c11; % Propene + heptene cracking
k_o15 = k_o11; % Propene + octene oligomerisation
k_c15 = k_c11; % Propene + octene cracking
k_o16 = k_o11; % Propene + nonene oligomerisation
k_c16 = k_c11; % Propene + nonene cracking
k_o17 = 1.69e-1; % Butene + butene oligomerisation
k_c17 = 6.01e-1; % Butene + butene cracking
k_o18 = k_o17; % Butene + pentene oligomerisation
k_c18 = k_c17; % Butene + pentene cracking
k_o19 = k_o17; % Butene + hexene oligomerisation
k_c19 = k_c17; % Butene + hexene cracking
k_o20 = k_o17; % Butene + heptene oligomerisation
k_c20 = k_c17; % Butene + heptene cracking
k_o21 = k_o17; % Butene + octene oligomerisation
k_c21 = k_c17; % Butene + octene cracking
k_o22 = k_o17; % Pentene + pentene oligomerisation
k_c22 = k_c17; % Pentene + pentene cracking
k_o23 = k_o22; % Pentene + hexene oligomerisation
k_c23 = k_c22; % Pentene + hexene cracking
k_o24 = k_o22; % Pentene + heptene oligomerisation
k_c24 = k_c22; % Pentene + heptene cracking
k_o25 = k_o22; % Hexene + hexene oligomerisation
k_c25 = k_c22; % Hexene + hexene cracking

k_ar1 = 1.16e-5; % Hexene + ethene aromatisation
k_ar2 = 9.31e-4; % Hexene + propene aromatisation
k_ar3 = 1.31e-3; % Hexene + butene aromatisation
k_ar4 = k_ar3; % Hexene + pentene aromatisation
k_ar5 = k_ar3; % Hexene + hexene aromatisation
k_ar6 = k_ar3; % Hexene + heptene aromatisation
k_ar7 = k_ar3; % Hexene + octene aromatisation
k_ar8 = k_ar3; % Hexene + nonene aromatisation
k_ar9 = k_ar3; % Hexene + decene aromatisation
k_ar10 = 3.18e-4; % Heptene + ethene aromatisation
k_ar11 = 2.54e-2; % Heptene + propene aromatisation
k_ar12 = 3.60e-2; % Heptene + butene aromatisation
k_ar13 = k_ar12; % Heptene + pentene aromatisation
k_ar14 = k_ar12; % Heptene + hexene aromatisation
k_ar15 = k_ar12; % Heptene + heptene aromatisation
k_ar16 = k_ar12; % Heptene + octene aromatisation
k_ar17 = k_ar12; % Heptene + nonene aromatisation
k_ar18 = k_ar12; % Heptene + decene aromatisation
k_ar19 = 2.38e-3; % Octene + ethene aromatisation
k_ar20 = 1.90e-1; % Octene + propene aromatisation
k_ar21 = 2.69e-1; % Octene + butene aromatisation
k_ar22 = k_ar21; % Octene + pentene aromatisation
k_ar23 = k_ar21; % Octene + hexene aromatisation
k_ar24 = k_ar21; % Octene + heptene aromatisation
k_ar25 = k_ar21; % Octene + octene aromatisation
k_ar26 = k_ar21; % Octene + nonene aromatisation
k_ar27 = k_ar21; % Octene + decene aromatisation
k_ar28 = 2.76e-3; % Nonene + ethene aromatisation
k_ar29 = 2.21e-1; % Nonene + propene aromatisation
k_ar30 = 3.11e-1; % Nonene + butene aromatisation
k_ar31 = k_ar30; % Nonene + pentene aromatisation
k_ar32 = k_ar30; % Nonene + hexene aromatisation
k_ar33 = k_ar30; % Nonene + heptene aromatisation
k_ar34 = k_ar30; % Nonene + octene aromatisation
k_ar35 = k_ar30; % Nonene + nonene aromatisation
k_ar36 = k_ar30; % Nonene + decene aromatisation
k_ar37 = 2.08e-3; % Decene + ethene aromatisation
k_ar38 = 1.66e-1; % Decene + propene aromatisation
k_ar39 = 2.34e-1; % Decene + butene aromatisation
k_ar40 = k_ar39; % Decene + pentene aromatisation
k_ar41 = k_ar39; % Decene + hexene aromatisation
k_ar42 = k_ar39; % Decene + heptene aromatisation
k_ar43 = k_ar39; % Decene + octene aromatisation
k_ar44 = k_ar39; % Decene + nonene aromatisation
k_ar45 = k_ar39; % Decene + decene aromatisation

k_m6a = 0; % Bezene methylation
k_m7a = 0; % Toluene methylation
k_m8a = 0; % Ethylbenzene methylation
k_m9a = 0; % Propylbenzene methylation

k_pcH3 = 0; % Propane protolytic cracking forming hydrogen
k_pcH4 = 0; % Butene protolytic cracking forming hydrogen
k_pcH5 = 0; % Pentane protolytic cracking forming hydrogen
k_pcH6 = 0; % Hexane protolytic cracking forming hydrogen
k_pcH7 = 0; % Heptane protolytic cracking forming hydrogen
k_pcH8 = 0; % Octane protolytic cracking forming hydrogen
k_pcC1 = 0; % Propane protolytic cracking forming methane
k_pcC2 = 0; % Butane protolytic cracking forming ethane
k_pcC3 = 0; % Butane protolytic cracking forming methane
k_pcC4 = 0; % Pentane protolytic cracking forming propane
k_pcC5 = 0; % Pentane protolytic cracking forming ethane
k_pcC6 = 0; % Pentane protolytic cracking forming methane
k_pcC7 = 0; % Hexane protolytic cracking forming butane
k_pcC8 = 0; % Hexane protolytic cracking forming propane
k_pcC9 = 0; % Hexane protolytic cracking forming ethane
k_pcC10 = 0; % Hexane protolytic cracking forming methane
k_pcC11 = 0; % Heptane protolytic cracking forming pentane
k_pcC12 = 0; % Heptane protolytic cracking forming butane
k_pcC13 = 0; % Heptane protolytic cracking forming propane
k_pcC14 = 0; % Heptane protolytic cracking forming ethane
k_pcC15 = 0; % Heptane protolytic cracking forming methane
k_pcC16 = 0; % Octane protolytic cracking forming hexane
k_pcC17 = 0; % Octane protolytic cracking forming pentane
k_pcC18 = 0; % Octane protolytic cracking forming butane
k_pcC19 = 0; % Octane protolytic cracking forming propane
k_pcC20 = 0; % Octane protolytic cracking forming ethane
k_pcC21 = 0; % Octane protolytic cracking forming methane

% Expression for [HZ]
HZ = 1/(1 + Ka1*P_CH3OH + Ka2*P_CH3OCH3 + Ka3*P_H2O + Ka4*P_CO...
    + (Ka5*P_CH4 + Ka5*P_C2H6 + Ka5*P_C3H8 + Ka5*P_C4H10...
    + Ka5*P_C5H12 + Ka5*P_C6H14 + Ka5*P_C7H16 + Ka5*P_C8H18...
    + Ka5*P_C9H20 + Ka5*P_C10H22) + (Ka6*P_C2H4 + Ka6*P_C3H6...
    + Ka6*P_C4H8 + Ka6*P_C5H10 + Ka6*P_C6H12 + Ka6*P_C7H14...
    + Ka6*P_C8H16 + Ka6*P_C9H18 + Ka6*P_C10H20 + Ka6*P_C11H22...
    + Ka6*P_C12H24) + (Ka7*P_C6H6 + Ka7*P_C7H8 + Ka7*P_C8H10...
    + Ka7*P_C9H12 + Ka7*P_C10H14) + Ka8*P_H2); % Fraction of empty active 
    %  sites

% Expression for [CH3COZ]
CH3COZ = (k_c*Ka1*P_CO*HZ)/(k_ef + k_pf); % Surface acetyl group 
    % concentration
%CH3COZ = 1e-100;

% Rate Equations (mol/g/s)
r_df = k_df*Ka1*(P_CH3OH^2)*HZ; % Methanol dehydration
r_db = k_db*Ka2*P_H2O*P_CH3OCH3*HZ; % DME hydration

r_CO = k_CO*Ka1*P_CH3OCH3*P_CH3OH*HZ; % Carbon monoxide formation

r_c = k_c*Ka1*P_CO*P_CH3OH*HZ; % Surface acetyl group formation

r_ef = k_ef*P_CH3OH*CH3COZ; % Ethene formation

r_pf = k_pf*P_CH3OH*CH3COZ; % Propene formation

r_m2o = k_m2o*Ka1*P_C2H4*P_CH3OH*HZ; % Ethene methylation
r_m3o = k_m3o*Ka1*P_C3H6*P_CH3OH*HZ; % Propene methylation
r_m4o = k_m4o*Ka1*P_C4H8*P_CH3OH*HZ; % Butene methylation
r_m5o = k_m5o*Ka1*P_C5H10*P_CH3OH*HZ; % Pentene methylation
r_m6o = k_m6o*Ka1*P_C6H12*P_CH3OH*HZ; % Hexene methylation
r_m7o = k_m7o*Ka1*P_C7H14*P_CH3OH*HZ; % Heptene methylation
r_m8o = k_m8o*Ka1*P_C8H16*P_CH3OH*HZ; % Octene methylation
r_m9o = k_m9o*Ka1*P_C9H18*P_CH3OH*HZ; % Nonene methylation
r_m10o = k_m10o*Ka1*P_C10H20*P_CH3OH*HZ; % Decene methylation
r_m11o = k_m11o*Ka1*P_C11H22*P_CH3OH*HZ; % Undecene methylation

r_o1 = k_o1*Ka6*(P_C2H4^2)*HZ; % Ethene + ethene oligomerisation
r_c1 = k_c1*Ka6*P_C4H8*HZ; % Ethene + ethene cracking
r_o2 = k_o2*Ka6*P_C2H4*P_C3H6*HZ; % Ethene + propene oligomerisation
r_c2 = k_c2*Ka6*P_C5H10*HZ; % Ethene + propene cracking
r_o3 = k_o3*Ka6*P_C2H4*P_C4H8*HZ; % Ethene + butene oligomerisation
r_c3 = k_c3*Ka6*P_C6H12*HZ; % Ethene + butene cracking
r_o4 = k_o4*Ka6*P_C2H4*P_C5H10*HZ; % Ethene + pentene oligomerisation
r_c4 = k_c4*Ka6*P_C7H14*HZ; % Ethene + pentene cracking
r_o5 = k_o5*Ka6*P_C2H4*P_C6H12*HZ; % Ethene + hexene oligomerisation
r_c5 = k_c5*Ka6*P_C8H16*HZ; % Ethene + hexene cracking
r_o6 = k_o6*Ka6*P_C2H4*P_C7H14*HZ; % Ethene + heptene oligomerisation
r_c6 = k_c6*Ka6*P_C9H18*HZ; % Ethene + heptene cracking
r_o7 = k_o7*Ka6*P_C2H4*P_C8H16*HZ; % Ethene + octene oligomerisation
r_c7 = k_c7*Ka6*P_C10H20*HZ; % Ethene + octene cracking
r_o8 = k_o8*Ka6*P_C2H4*P_C9H18*HZ; % Ethene + nonene oligomerisation
r_c8 = k_c8*Ka6*P_C11H22*HZ; % Ethene + nonene cracking
r_o9 = k_o9*Ka6*P_C2H4*P_C10H20*HZ; % Ethene + decene oligomerisation
r_c9 = k_c9*Ka6*P_C12H24*HZ; % Ethene + decene cracking
r_o10 = k_o10*Ka6*(P_C3H6^2)*HZ; % Propene + propene oligomerisation
r_c10 = k_c10*Ka6*P_C6H12*HZ; % Propene + propene cracking
r_o11 = k_o11*Ka6*P_C3H6*P_C4H8*HZ; % Propene + butene oligomerisation
r_c11 = k_c11*Ka6*P_C7H14*HZ; % Propene + butene cracking
r_o12 = k_o12*Ka6*P_C3H6*P_C5H10*HZ; % Propene + pentene oligomerisation
r_c12 = k_c12*Ka6*P_C8H16*HZ; % Propene + pentene cracking
r_o13 = k_o13*Ka6*P_C3H6*P_C6H12*HZ; % Propene + hexene oligomerisation
r_c13 = k_c13*Ka6*P_C9H18*HZ; % Propene + hexene cracking
r_o14 = k_o14*Ka6*P_C3H6*P_C7H14*HZ; % Propene + heptene oligomerisation
r_c14 = k_c14*Ka6*P_C10H20*HZ; % Propene + heptene cracking
r_o15 = k_o15*Ka6*P_C3H6*P_C8H16*HZ; % Propene + octene oligomerisation
r_c15 = k_c15*Ka6*P_C11H22*HZ; % Propene + octene cracking
r_o16 = k_o16*Ka6*P_C3H6*P_C9H18*HZ; % Propene + nonene oligomerisation
r_c16 = k_c16*Ka6*P_C12H24*HZ; % Propene + nonene cracking
r_o17 = k_o17*Ka6*(P_C4H8^2)*HZ; % Butene + butene oligomerisation
r_c17 = k_c17*Ka6*P_C8H16*HZ; % Butene + butene cracking
r_o18 = k_o18*Ka6*P_C4H8*P_C5H10*HZ; % Butene + pentene oligomerisation
r_c18 = k_c18*Ka6*P_C9H18*HZ; % Butene + pentene cracking
r_o19 = k_o19*Ka6*P_C4H8*P_C6H12*HZ; % Butene + hexene oligomerisation
r_c19 = k_c19*Ka6*P_C10H20*HZ; % Butene + hexene cracking
r_o20 = k_o20*Ka6*P_C4H8*P_C7H14*HZ; % Butene + heptene oligomerisation
r_c20 = k_c20*Ka6*P_C11H22*HZ; % Butene + heptene cracking
r_o21 = k_o21*Ka6*P_C4H8*P_C8H16*HZ; % Butene + octene oligomerisation
r_c21 = k_c21*Ka6*P_C12H24*HZ; % Butene + octene cracking
r_o22 = k_o22*Ka6*(P_C5H10^2)*HZ; % Pentene + pentene oligomerisation
r_c22 = k_c22*Ka6*P_C10H20*HZ; % Pentene + pentene cracking
r_o23 = k_o23*Ka6*P_C5H10*P_C6H12*HZ; % Pentene + hexene oligomerisation
r_c23 = k_c23*Ka6*P_C11H22*HZ; % Pentene + hexene cracking
r_o24 = k_o24*Ka6*P_C5H10*P_C7H14*HZ; % Pentene + heptene oligomerisation
r_c24 = k_c24*Ka6*P_C12H24*HZ; % Pentene + heptene cracking
r_o25 = k_o25*Ka6*(P_C6H12^2)*HZ; % Hexene + hexene oligomerisation
r_c25 = k_c25*Ka6*P_C12H24*HZ; % Hexene + hexene cracking

r_ar1 = k_ar1*Ka6*P_C6H12*P_C2H4*HZ; % Hexene + ethene aromatisation
r_ar2 = k_ar2*Ka6*P_C6H12*P_C3H6*HZ; % Hexene + propene aromatisation
r_ar3 = k_ar3*Ka6*P_C6H12*P_C4H8*HZ; % Hexene + butene aromatisation
r_ar4 = k_ar4*Ka6*P_C6H12*P_C5H10*HZ; % Hexene + pentene aromatisation
r_ar5 = k_ar5*Ka6*(P_C6H12^2)*HZ; % Hexene + hexene aromatisation
r_ar6 = k_ar6*Ka6*P_C6H12*P_C7H14*HZ; % Hexene + heptene aromatisation
r_ar7 = k_ar7*Ka6*P_C6H12*P_C8H16*HZ; % Hexene + octene aromatisation
r_ar8 = k_ar8*Ka6*P_C6H12*P_C9H18*HZ; % Hexene + nonene aromatisation
r_ar9 = k_ar9*Ka6*P_C6H12*P_C10H20*HZ; % Hexene + decene aromatisation
r_ar10 = k_ar10*Ka6*P_C7H14*P_C2H4*HZ; % Heptene + ethene aromatisation
r_ar11 = k_ar11*Ka6*P_C7H14*P_C3H6*HZ; % Heptene + propene aromatisation
r_ar12 = k_ar12*Ka6*P_C7H14*P_C4H8*HZ; % Heptene + butene aromatisation
r_ar13 = k_ar13*Ka6*P_C7H14*P_C5H10*HZ; % Heptene + pentene aromatisation
r_ar14 = k_ar14*Ka6*P_C7H14*P_C6H12*HZ; % Heptene + hexene aromatisation
r_ar15 = k_ar15*Ka6*(P_C7H14^2)*HZ; % Heptene + heptene aromatisation
r_ar16 = k_ar16*Ka6*P_C7H14*P_C8H16*HZ; % Heptene + octene aromatisation
r_ar17 = k_ar17*Ka6*P_C7H14*P_C9H18*HZ; % Heptene + nonene aromatisation
r_ar18 = k_ar18*Ka6*P_C7H14*P_C10H20*HZ; % Heptene + decene aromatisation
r_ar19 = k_ar19*Ka6*P_C8H16*P_C2H4*HZ; % Octene + ethene aromatisation
r_ar20 = k_ar20*Ka6*P_C8H16*P_C3H6*HZ; % Octene + propene aromatisation
r_ar21 = k_ar21*Ka6*P_C8H16*P_C4H8*HZ; % Octene + butene aromatisation
r_ar22 = k_ar22*Ka6*P_C8H16*P_C5H10*HZ; % Octene + pentene aromatisation
r_ar23 = k_ar23*Ka6*P_C8H16*P_C6H12*HZ; % Octene + hexene aromatisation
r_ar24 = k_ar24*Ka6*P_C8H16*P_C7H14*HZ; % Octene + heptene aromatisation
r_ar25 = k_ar25*Ka6*(P_C8H16^2)*HZ; % Octene + octene aromatisation
r_ar26 = k_ar26*Ka6*P_C8H16*P_C9H18*HZ; % Octene + nonene aromatisation
r_ar27 = k_ar27*Ka6*P_C8H16*P_C10H20*HZ; % Octene + decene aromatisation
r_ar28 = k_ar28*Ka6*P_C9H18*P_C2H4*HZ; % Nonene + ethene aromatisation
r_ar29 = k_ar29*Ka6*P_C9H18*P_C3H6*HZ; % Nonene + propene aromatisation
r_ar30 = k_ar30*Ka6*P_C9H18*P_C4H8*HZ; % Nonene + butene aromatisation
r_ar31 = k_ar31*Ka6*P_C9H18*P_C5H10*HZ; % Nonene + pentene aromatisation
r_ar32 = k_ar32*Ka6*P_C9H18*P_C6H12*HZ; % Nonene + hexene aromatisation
r_ar33 = k_ar33*Ka6*P_C9H18*P_C7H14*HZ; % Nonene + heptene aromatisation
r_ar34 = k_ar34*Ka6*P_C9H18*P_C8H16*HZ; % Nonene + octene aromatisation
r_ar35 = k_ar35*Ka6*(P_C9H18^2)*HZ; % Nonene + nonene aromatisation
r_ar36 = k_ar36*Ka6*P_C9H18*P_C10H20*HZ; % Nonene + decene aromatisation
r_ar37 = k_ar37*Ka6*P_C10H20*P_C2H4*HZ; % Decene + ethene aromatisation
r_ar38 = k_ar38*Ka6*P_C10H20*P_C3H6*HZ; % Decene + propene aromatisation
r_ar39 = k_ar39*Ka6*P_C10H20*P_C4H8*HZ; % Decene + butene aromatisation
r_ar40 = k_ar40*Ka6*P_C10H20*P_C5H10*HZ; % Decene + pentene aromatisation
r_ar41 = k_ar41*Ka6*P_C10H20*P_C6H12*HZ; % Decene + hexene aromatisation
r_ar42 = k_ar42*Ka6*P_C10H20*P_C7H14*HZ; % Decene + heptene aromatisation
r_ar43 = k_ar43*Ka6*P_C10H20*P_C8H16*HZ; % Decene + octene aromatisation
r_ar44 = k_ar44*Ka6*P_C10H20*P_C9H18*HZ; % Decene + nonene aromatisation
r_ar45 = k_ar45*Ka6*(P_C10H20^2)*HZ; % Decene + decene aromatisation

r_m6a = k_m6a*Ka7*P_CH3OH*P_C6H6*HZ; % Bezene methylation
r_m7a = k_m7a*Ka7*P_CH3OH*P_C7H8*HZ; % Toluene methylation
r_m8a = k_m8a*Ka7*P_CH3OH*P_C8H10*HZ; % Ethylbenzene methylation
r_m9a = k_m9a*Ka7*P_CH3OH*P_C9H12*HZ; % Propylbenzene methylation

r_pcH3 = k_pcH3*P_C3H8*HZ; % Propane protolytic cracking forming hydrogen
r_pcH4 = k_pcH4*P_C4H10*HZ; % Butane protolytic cracking forming hydrogen
r_pcH5 = k_pcH5*P_C5H12*HZ; % Pentane protolytic cracking forming hydrogen
r_pcH6 = k_pcH6*P_C6H14*HZ; % Hexane protolytic cracking forming hydrogen
r_pcH7 = k_pcH7*P_C7H16*HZ; % Heptane protolytic cracking forming hydrogen
r_pcH8 = k_pcH8*P_C8H18*HZ; % Octane protolytic cracking forming hydrogen
r_pcC1 = k_pcC1*P_C3H8*HZ; % Propane protolytic cracking forming methane
r_pcC2 = k_pcC2*P_C4H10*HZ; % Butane protolytic cracking forming ethane
r_pcC3 = k_pcC3*P_C4H10*HZ; % Butane protolytic cracking forming methane
r_pcC4 = k_pcC4*P_C5H12*HZ; % Pentane protolytic cracking forming propane
r_pcC5 = k_pcC5*P_C5H12*HZ; % Pentane protolytic cracking forming ethane
r_pcC6 = k_pcC6*P_C5H12*HZ; % Pentane protolytic cracking forming methane
r_pcC7 = k_pcC7*P_C6H14*HZ; % Hexane protolytic cracking forming Butane
r_pcC8 = k_pcC8*P_C6H14*HZ; % Hexane protolytic cracking forming Propane
r_pcC9 = k_pcC9*P_C6H14*HZ; % Hexane protolytic cracking forming ethane
r_pcC10 = k_pcC10*P_C6H14*HZ; % Hexane protolytic cracking forming methane
r_pcC11 = k_pcC11*P_C7H16*HZ; % Heptane protolytic cracking forming pentane
r_pcC12 = k_pcC12*P_C7H16*HZ; % Heptane protolytic cracking forming butane
r_pcC13 = k_pcC13*P_C7H16*HZ; % Heptane protolytic cracking forming propane
r_pcC14 = k_pcC14*P_C7H16*HZ; % Heptane protolytic cracking forming ethane
r_pcC15 = k_pcC15*P_C7H16*HZ; % Heptane protolytic cracking forming methane
r_pcC16 = k_pcC16*P_C8H18*HZ; % Octane protolytic cracking forming hexane
r_pcC17 = k_pcC17*P_C8H18*HZ; % Octane protolytic cracking forming pentane
r_pcC18 = k_pcC18*P_C8H18*HZ; % Octane protolytic cracking forming butane
r_pcC19 = k_pcC19*P_C8H18*HZ; % Octane protolytic cracking forming propane
r_pcC20 = k_pcC20*P_C8H18*HZ; % Octane protolytic cracking forming ethane
r_pcC21 = k_pcC21*P_C8H18*HZ; % Octane protolytic cracking forming methane

% Species Mass Balances
R_CH3OH = -2*r_df + 2*r_db - r_CO - r_c - r_ef - 2*r_pf - r_m2o...
    - r_m3o - r_m4o - r_m5o - r_m6o - r_m7o - r_m8o - r_m9o - r_m10o...
    - r_m11o - r_m6a - r_m7a - r_m8a - r_m9a; % Methanol
R_CH3OCH3 = r_df - r_db - r_CO; % DME
R_H2O = r_df - r_db + r_CO + r_c + r_ef + 2*r_pf + r_m2o + r_m3o...
    + r_m4o + r_m5o + r_m6o + r_m7o + r_m8o + r_m9o + r_m10o + r_m11o...
    + r_m6a + r_m7a + r_m8a + r_m9a; % Water
R_CO = r_CO - r_c + r_ef + r_pf; % Carbon Monoxide
R_CH4 = 2*r_CO + r_pcC1 + r_pcC3 + r_pcC6 + r_pcC10 + r_pcC15...
    + r_pcC21; % Methane
R_C2H4 = r_ef - r_m2o - 2*r_o1 + 2*r_c1 - r_o2 + r_c2 - r_o3 + r_c3...
    - r_o4 + r_c4 - r_o5 + r_c5 - r_o6 + r_c6 - r_o7 + r_c7 - r_o8...
    + r_c8 - r_o9 + r_c9 - 3*r_ar1 - 3*r_ar10 - 3*r_ar19 - 3*r_ar28...
    - 3*r_ar37 + r_pcC1 + r_pcC2 + r_pcC4 + r_pcC7 + r_pcC11...
    + r_pcC16; % Ethene
R_C2H6 = 3*r_ar1 + 3*r_ar10 + 3*r_ar19 + 3*r_ar28 + 3*r_ar37 + r_pcC2...
    + r_pcC5 + r_pcC9 + r_pcC14 + r_pcC20; % Ethane
R_C3H6 = r_pf + r_m2o - r_m3o - r_o2 + r_c2 - 2*r_o10 + 2*r_c10...
    - r_o11 + r_c11 - r_o12 + r_c12 - r_o13 + r_c13 - r_o14 + r_c14...
    - r_o15 + r_c15 - r_o16 + r_c16 - 3*r_ar2 - 3*r_ar11 - 3*r_ar20...
    - 3*r_ar29 - 3*r_ar38 + r_pcC3 + r_pcC5 + r_pcC8 + r_pcC12...
    + r_pcC17; % Propene
R_C3H8 = 3*r_ar2 + 3*r_ar11 + 3*r_ar20 + 3*r_ar29 + 3*r_ar38 - r_pcH3...
    - r_pcC1 + r_pcC4 + r_pcC8 + r_pcC13 + r_pcC19; % Propane
R_C4H8 = r_m3o - r_m4o + r_o1 - r_c1 - r_o3 + r_c3 - r_o11 + r_c11...
    - 2*r_o17 + 2*r_c17 - r_o18 + r_c18 - r_o19 + r_c19 - r_o20...
    + r_c20 - r_o21 + r_c21 - 3*r_ar3 - 3*r_ar12 - 3*r_ar21...
    - 3*r_ar30 - 3*r_ar39 + r_pcC6 + r_pcC9 + r_pcC13 + r_pcC18;...
    % Butene
R_C4H10 = 3*r_ar3 + 3*r_ar12 + 3*r_ar21 + 3*r_ar30 + 3*r_ar39...
    - r_pcH4 - r_pcC2 - r_pcC3 + r_pcC7 + r_pcC12 + r_pcC18; % Butane
R_C5H10 = r_m4o - r_m5o + r_o2 - r_c2 - r_o4 + r_c4 - r_o12 + r_c12...
    - r_o18 + r_c18 - 2*r_o22 + 2*r_c22 - r_o23 + r_c23 - r_o24...
    + r_c24 - 3*r_ar4 - 3*r_ar13 - 3*r_ar22 - 3*r_ar31 - 3*r_ar40...
    + r_pcC10 + r_pcC14 + r_pcC19; % Pentene
R_C5H12 = 3*r_ar4 + 3*r_ar13 + 3*r_ar22 + 3*r_ar31 + 3*r_ar40...
    - r_pcH5 - r_pcC4 - r_pcC5 - r_pcC6 + r_pcC11 + r_pcC17; % Pentane
R_C6H6 = r_ar1 + r_ar2 + r_ar3 + r_ar4 + r_ar5 + r_ar6 + r_ar7...
    + r_ar8 + r_ar9 - r_m6a; % Benzene
R_C6H12 = r_m5o - r_m6o + r_o3 - r_c3 - r_o5 + r_c5 + r_o10 - r_c10...
    - r_o13 + r_c13 - r_o19 + r_c19 - r_o23 + r_c23 - 2*r_o25...
    + 2*r_c25 - r_ar1 - r_ar2 - r_ar3 - r_ar4 - 4*r_ar5 - r_ar6...
    - r_ar7 - r_ar8 - r_ar9 - 3*r_ar14 - 3*r_ar23 - 3*r_ar32...
    - 3*r_ar41 + r_pcC15 + r_pcC20; % Hexene
R_C6H14 = 3*r_ar5 + 3*r_ar14 + 3*r_ar23 + 3*r_ar32 + 3*r_ar41...
    - r_pcH6 - r_pcC7 - r_pcC8 - r_pcC9 - r_pcC10 + r_pcC16; % Hexane
R_C7H8 = r_ar10 + r_ar11 + r_ar12 + r_ar13 + r_ar14 + r_ar15 + r_ar16...
    + r_ar17 + r_ar18 - r_m7a; % Toluene
R_C7H14 = r_m6o - r_m7o + r_o4 - r_c4 - r_o6 + r_c6 + r_o11 - r_c11...
    - r_o14 + r_c14 - r_o20 + r_c20 - r_o24 + r_c24 - 3*r_ar6...
    - r_ar10 - r_ar11 - r_ar12 - r_ar13 - r_ar14 - 4*r_ar15 - r_ar16...
    - r_ar17 - r_ar18 - 3*r_ar24 - 3*r_ar33 - 3*r_ar42 + r_pcC21;...
    % Heptene
R_C7H16 = 3*r_ar6 + 3*r_ar15 + 3*r_ar24 + 3*r_ar33 + 3*r_ar42...
    - r_pcH7 - r_pcC11 - r_pcC12 - r_pcC13 - r_pcC14 - r_pcC15;...
    % Heptane
R_C8H10 = r_ar19 + r_ar20 + r_ar21 + r_ar22 + r_ar23 + r_ar24...
    + r_ar25 + r_ar26 + r_ar27 - r_m8a; % Ethylbenzene
R_C8H16 = r_m7o - r_m8o + r_o5 - r_c5 - r_o7 + r_c7 + r_o12 - r_c12...
    - r_o15 + r_c15 + r_o17 - r_c17 - r_o21 + r_c21 - 3*r_ar7...
    - 3*r_ar16 - r_ar19 - r_ar20 - r_ar21 - r_ar22 - r_ar23 - r_ar24...
    - 4*r_ar25 - r_ar26 - r_ar27 - 3*r_ar34 - 3*r_ar43; % Octene
R_C8H18 = 3*r_ar7 + 3*r_ar16 + 3*r_ar25 + 3*r_ar34 + 3*r_ar43...
    - r_pcH8 - r_pcC16 - r_pcC17 - r_pcC18 - r_pcC19 - r_pcC20...
    - r_pcC21; % Octane
R_C9H12 = r_ar28 + r_ar29 + r_ar30 + r_ar31 + r_ar32 + r_ar33...
    + r_ar34 + r_ar35 + r_ar36 - r_m9a; % Propylbenzene
R_C9H18 = r_m8o - r_m9o + r_o6 - r_c6 - r_o8 + r_c8 + r_o13 - r_c13...
    - r_o16 + r_c16 + r_o18 - r_c18 - 3*r_ar8 - 3*r_ar17 - 3*r_ar26...
    - r_ar28 - r_ar29 - r_ar30 - r_ar31 - r_ar32 - r_ar33 - r_ar34...
    - 4*r_ar35 - r_ar36 - 3*r_ar44; % Nonene
R_C9H20 = 3*r_ar8 + 3*r_ar17 + 3*r_ar26 + 3*r_ar35 + 3*r_ar44; % Nonane
R_C10H14 = r_ar37 + r_ar38 + r_ar39 + r_ar40 + r_ar41 + r_ar42...
    + r_ar43 + r_ar44 + r_ar45; % Butylbenzene
R_C10H20 = r_m9o - r_m10o + r_o7 - r_c7 - r_o9 + r_c9 + r_o14 - r_c14...
    + r_o19 - r_c19 + r_o22 - r_c22 - 3*r_ar9 - 3*r_ar18 - 3*r_ar27...
    - 3*r_ar36 - r_ar37 - r_ar38 - r_ar39 - r_ar40 - r_ar41 - r_ar42...
    - r_ar43 - r_ar44 - 4*r_ar45; % Decene
R_C10H22 = 3*r_ar9 + 3*r_ar18 + 3*r_ar27 + 3*r_ar36 + 3*r_ar45; % Decane
R_C11H22 = r_m10o - r_m11o + r_o8 - r_c8 + r_o15 - r_c15 + r_o20...
    - r_c20 + r_o23 - r_c23; % Undecene
R_C12H24 = r_m11o + r_o9 - r_c9 + r_o16 - r_c16 + r_o21 - r_c21...
    + r_o24 - r_c24 + r_o25 - r_c25; % Dodecene
R_H2 = r_pcH3 + r_pcH4 + r_pcH5 + r_pcH6 + r_pcH7 + r_pcH8; % Hydrogen
R_N2 = 0; % Nitrogen

% Zero Column Vector
dmdt = zeros(32,1);

% ODEs
dmdt(1) = R_CH3OH*MW_CH3OH; % Methanol
dmdt(2) = R_CH3OCH3*MW_CH3OCH3; % DME
dmdt(3) = R_H2O*MW_H2O; % Water
dmdt(4) = R_CO*MW_CO; % Carbon Monoxide
dmdt(5) = R_CH4*MW_CH4; % Methane
dmdt(6) = R_C2H4*MW_C2H4; % Ethene
dmdt(7) = R_C2H6*MW_C2H6; % Ethane
dmdt(8) = R_C3H6*MW_C3H6; % Propene
dmdt(9) = R_C3H8*MW_C3H8; % Propane
dmdt(10) = R_C4H8*MW_C4H8; % Butene
dmdt(11) = R_C4H10*MW_C4H10; % Butane
dmdt(12) = R_C5H10*MW_C5H10; % Pentene
dmdt(13) = R_C5H12*MW_C5H12; % Pentane
dmdt(14) = R_C6H6*MW_C6H6; % Benzene
dmdt(15) = R_C6H12*MW_C6H12; % Hexene
dmdt(16) = R_C6H14*MW_C6H14; % Hexane
dmdt(17) = R_C7H8*MW_C7H8; % Toluene
dmdt(18) = R_C7H14*MW_C7H14; % Heptene
dmdt(19) = R_C7H16*MW_C7H16; % Heptane
dmdt(20) = R_C8H10*MW_C8H10; % Ethylbenzene
dmdt(21) = R_C8H16*MW_C8H16; % Octene
dmdt(22) = R_C8H18*MW_C8H18; % Octane
dmdt(23) = R_C9H12*MW_C9H12; % Propylbenzene
dmdt(24) = R_C9H18*MW_C9H18; % Nonene
dmdt(25) = R_C9H20*MW_C9H20; % Nonane
dmdt(26) = R_C10H14*MW_C10H14; % Butylbenzene
dmdt(27) = R_C10H20*MW_C10H20; % Decene
dmdt(28) = R_C10H22*MW_C10H22; % Decane
dmdt(29) = R_C11H22*MW_C11H22; % Undecene
dmdt(30) = R_C12H24*MW_C12H24; % Dodecene
dmdt(31) = R_H2*MW_H2; % Hydrogen
dmdt(32) = R_N2*MW_N2; % Nitrogen

end