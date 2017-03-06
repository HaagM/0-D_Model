function out = MartNel_ATHL( m, p, h, d, l, Zeta)
%m: mass flux [kg/s], p: pressure [bar], h: enthalpie [kJ/kj], d: diameter [m], l: tube length [m],
%Zeta: friction loss coeff. [-], r: wall roughness [m]
%Calculates the pressure drop for 1- and 2-phase flow using the
%correlations which are implemented in ATHLET. 

%%

global  r 

% Stoffdaten
p_crit = 220.6;
v_l = XSteam('vL_p',p);
v_v = XSteam('vV_p',p);
v = XSteam('v_ph',p,h);
eta_l = XSteam('my_ph',p,XSteam('hL_p',p));
eta_v = XSteam('my_ph',p,XSteam('hV_p',p));
eta = XSteam('my_ph',p,h);
x = XSteam('x_ph',p,h);
% Rechengrößen
u = m * v / (pi/4 * d^2); %Strömungsgschwindigkeit
Re_l = max(u * d / eta_l / v_l , 500);
Re_v = max(u * d / eta_v / v_v , 1);
A = pi/4 * d^2;

%% Darcy-Weissbach friction factor lambda
lambda_lam_l = 64 / Re_l;
lambda_turb_l = 1 / (-2 * log10(r / (3.7 * d) + (6.81 / Re_l)^0.9))^2;
lambda_l = max(lambda_lam_l , lambda_turb_l);
lambda_v = 1 / (-2 * log10(r / (3.7 * d) + (6.81 / Re_v)^0.9))^2;

%% Two-phase friction factor 
if x == 0 
    k = lambda_l * l / (d * A^2) *0.5 * v_l;
else
    %% Martinelli Nelson
    X_tt = (v_l / v_l)^0.571 * ((eta_l / eta_v)^0.143 * (1 - x) / x);
    a = 1 + 1.53 * (p_crit - p) / p_crit;
    PHI2_tt_L = ((1 + X_tt^(1 / a))/X_tt^(1 / a))^(1.75*a);
    C_PHI_L_only = PHI2_tt_L * (1 - x)^1.75;
    k_PHI_L_only = 0.5 * lambda_l * l / ((pi/4 * d^2)^2 * d) * C_PHI_L_only * v_l;
 
    PHI2_tt_V = (1 + X_tt^(1 / a))^(1.75 * a);
    C_PHI_V_only = PHI2_tt_V * x^1.75;
    k_PHI_V_only = 0.5 * lambda_v * l / ((pi/4 * d^2)^2 * d) * C_PHI_V_only * v_v;
    %
    k = (1 - x) * k_PHI_L_only + x * k_PHI_V_only;
end
out = (k + 0.5 * v * Zeta / A^2) * m^2;
end
