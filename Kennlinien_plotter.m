clear;
clc;

%% 
path(path, 'Druckverlust');
%% input
global P_in dH d_dc l_dc Zeta_dc T_dc p_dc d_ris l_ris Zeta_ris p_ris g p0
%
T_ini = 40;
p = 1;
P_in = 200;
dH = 7.348;
r = 4e-5; % Wall roughness
% downcomer
d_dc = 41.8e-3;
l_dc = 11.436;
Zeta_dc = 7.202;
T_dc = T_ini;
p_dc = p;
% riser
d_ris = 53e-3;
l_ris = 7.375;
Zeta_ris = 1.5;
p_ris = p;
p0 = 1;

%% constants
g = 9.81;
r = 4e-5;
%% Stoffwerte
h_dc = XSteam('h_pT',p,T_dc);
v_dc = XSteam('v_pT',p,T_dc);

%% Druckdifferenzen
% Antreibende Druckdifferenz
i = 1;
for m = 0.05 : 0.0025 : 2
    h_ris = P_in / m + h_dc;
    p_ris_eq = @(dp_ris)((p0 - dp_ris)*1e5 + g * dH / XSteam('v_ph',dp_ris,h_ris)); %Berechnung des Druckes im SR 
    p_ris = fzero(p_ris_eq, [1 3]);
    sum_dp(i) = press_eq(m);
    
    v_ris = XSteam('v_ph',p_ris,h_ris);
    dp_acc(i) = dH * g * (1/v_dc - 1/v_ris);
    dp_fric_dc(i) = MartNel_ATHL(m,p_dc,h_dc,d_dc,l_dc,Zeta_dc);
    dp_fric_ris_MN(i) = MartNel_ATHL(m,p_ris,h_ris,d_ris,l_ris,Zeta_ris);
    dp_fric_ris_Ch(i) = Chis_ATHL(m,p_ris,h_ris,d_ris,l_ris,Zeta_ris,r);
    x(i) = XSteam('x_ph',p,h_ris) * 1e5;
    m_x(i) = m;
    i = i + 1;
end

dp_fric_MN = dp_fric_dc + dp_fric_ris_MN;
dp_fric_Ch = dp_fric_dc + dp_fric_ris_Ch;
%%Plot
%plot(m_x,dp_acc,m_x,dp_fric_MN,m_x,dp_fric_Ch,m_x,x)
%legend('dp acc','dp MN','dp Ch','x')
plot(m_x,dp_acc,m_x,dp_fric_MN,m_x,sum_dp)
legend('dp acc','dp MN','sum dp MN')
grid on
