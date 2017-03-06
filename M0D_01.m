clear;
clc;

global T_ini p0 P_in dH r d_dc l_dc Zeta_dc T_dc p_dc d_ris l_ris Zeta_ris p_ris g

%% input
%
T_ini = 40;
p0 = 1;
P_in = 80;
dH = 7.348;
r = 4e-5; % Wand Rauigkeit
% downcomer
d_dc = 41.8e-3;
l_dc = 11.436;
Zeta_dc = 7.202;
T_dc = T_ini;
p_dc = p0;
% riser
d_ris = 53e-3;
l_ris = 7.375;
Zeta_ris = 1.5;

%% Konstanten
g = 9.81;

Sum_dp_eq = @press_eq;

%% Path
path(path, 'Druckverlust');

%% Stoffwerte
v_dc = XSteam('v_pT',p_dc,T_dc);
h_dc = XSteam('h_pT',p_dc,T_dc);

interv = 0.02;
n_sol = 0;
Sum_dp = zeros(1,2);
for m = 0 : interv : 2 %Sucht in welchen Interval sich eine 0-Stelle befindet für Sum_dp 
    Sum_dp = [Sum_dp_eq(m) Sum_dp];
    Sum_dp(3) = [];
    if m == 0
        continue
    end
    if (Sum_dp(1) * Sum_dp(2)) > 0 || isnan(Sum_dp(1)) || isnan(Sum_dp(2)) %Prüft ob sich das Vorzeichen über ein Interval ändert
        continue
    elseif Sum_dp(:) == 0;
        n_sol = n_sol +1;
        m_sol(n_sol) = m;
    else  
        n_sol = n_sol + 1;
        m_sol(n_sol) = fzero(Sum_dp_eq,[m-interv m]); % 0-Stellensuche für Druckgleichgewicht im durchsuchten Interval
    end
end
        

