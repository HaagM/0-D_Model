function [Sum_dp, dp_acc, dp_fric, dp_fric_ris, H_VF] = press_eq(m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global P_in H d_dc l_dc Zeta_dc T_dc p_dc d_ris Zeta_ris g p0 h_ris

h_dc = XSteam('h_pT',p_dc,T_dc);
v_dc = XSteam('v_pT',p_dc,T_dc);

h_ris = P_in / m + h_dc;

if h_ris > 3200 %Gültigkeitsbereich für Wasserdampfmakro prüfen: 7350 kJ/kg enstpricht h(p=99bar, T=1999 °C) 
    Sum_dp = NaN;
    dp_acc = NaN;
    dp_fric = NaN;
    dp_fric_ris = NaN; 
    H_VF = NaN;
    return;
end

%Druck im Eintritt des SR (Ausgangspunkt: Wassersäulfe: gesät. Flüssigkeit)
p_ris_eq_b = @(p_ris_b)((p0 - p_ris_b) + (g * H / XSteam('vL_p',p_ris_b))*1e-5); 
p_ris_b = fzero(p_ris_eq_b, [1 3]);

%% Berechnung der Höhe der Verdampfungsfront
H_VF = H - ((XSteam('psat_T',XSteam('T_ph',p_ris_b,h_ris)) - p0) * 1e5 * XSteam('v_ph',p_ris_b,h_ris)) / g; %Berechnung der Höhe der Verdampfungsfront
if H_VF < 0 % Dampfbildung in Heizzone
    H_VF = 0;
elseif H_VF >= H % Einphasiger Umlauf
    H_VF = H;
end

%% 
p_ris_eq_v = @(p_ris_v)((p0 - p_ris_v)*1e5 + g * (H - H_VF) / XSteam('v_ph',p_ris_v,h_ris)); %Berechnung des Druckes im SR im oberen- (dampf-) sub-KV
p_ris_v = fzero(p_ris_eq_v, [1 2]);

%Druck der Flüssigphase auf p_ris_b festsetzten. Würde p_ris_l wie p_ris_v
%berechnet werden, ist der Druck zu niedrig und es kommt zu weiteren
%verdampfen der Flüssigphase (Flashing)

%so nicht: p_ris_eq_l = @(p_ris_l)((p_ris_v - p_ris_l)*1e5 + g * H_VF / XSteam('v_ph',p_ris_l,h_ris)); 
% p_ris_l = fzero(p_ris_eq_l, [1 3]);
p_ris_l = p_ris_b;

v_ris_v = XSteam('v_ph',p_ris_v,h_ris);
v_ris_l = XSteam('v_ph',p_ris_l,h_ris);
%Dampfgehalt nach Verdampfungsfront:
x_ris_v = XSteam('vx_ph',p_ris_v,h_ris);

dp_fric_ris_v = MartNel_ATHL(m,p_ris_v,h_ris,d_ris,(H-H_VF),(Zeta_ris / 2));
dp_fric_ris_l = MartNel_ATHL(m,p_ris_l,h_ris,d_ris,H_VF,(Zeta_ris / 2));

%% Druckdifferenzen    
dp_acc = g * (H/v_dc - (H_VF / v_ris_l + (H-H_VF) / v_ris_v));  %Antreibende Druckdifferenz
dp_fric_dc = MartNel_ATHL(m,p_dc,h_dc,d_dc,l_dc,Zeta_dc);	%Reibungsdruckverlust Fallrohr

dp_fric_ris = dp_fric_ris_v + dp_fric_ris_l;    
dp_fric = dp_fric_dc + dp_fric_ris; 
Sum_dp = dp_acc - dp_fric;

end

