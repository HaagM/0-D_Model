function Sum_dp = press_eq(m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global P_in dH d_dc l_dc Zeta_dc T_dc p_dc d_ris l_ris Zeta_ris p_ris g p0

h_dc = XSteam('h_pT',p_dc,T_dc);
v_dc = XSteam('v_pT',p_dc,T_dc);

h_ris = P_in / m + h_dc;
if h_ris > 7350 %Gültigkeitsbereich für Wasserdampfmakro prüfen: 7350 kJ/kg enstpricht h(p=99bar, T=1999 °C) 
    Sum_dp = NaN;
    return;
end

p_ris_eq = @(dp_ris)((p0 - dp_ris)*1e5 + g * dH / XSteam('v_ph',dp_ris,h_ris)); %Berechnung des Druckes im SR 
p_ris = fzero(p_ris_eq, [1 3]);

v_ris = XSteam('v_ph',p_ris,h_ris);
dp_acc = dH * g * (1/v_dc - 1/v_ris); % Antreibende Druckdifferenz
dp_fric_dc = MartNel_ATHL(m,p_dc,h_dc,d_dc,l_dc,Zeta_dc); %Reibungsdruckverlust Fallrohr
dp_fric_ris = MartNel_ATHL(m,p_ris,h_ris,d_ris,l_ris,Zeta_ris); %Reibungsdruckverlust Steigrohr
 
Sum_dp = dp_acc - (dp_fric_dc + dp_fric_ris);

end

