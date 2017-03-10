function Sum_dp = press_eq(m)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

global P_in H d_dc l_dc Zeta_dc T_dc p_dc d_ris Zeta_ris g p0 h_ris

h_dc = XSteam('h_pT',p_dc,T_dc);
v_dc = XSteam('v_pT',p_dc,T_dc);

h_ris = P_in / m + h_dc;
if h_ris > 7350 %Gültigkeitsbereich für Wasserdampfmakro prüfen: 7350 kJ/kg enstpricht h(p=99bar, T=1999 °C) 
    Sum_dp = NaN;
    return;
end

p_ris_b = p0 + (g * H / XSteam('v_pT',1.8,80))*1e-5;

if h_ris < XSteam('hL_p', p0)
    H_VF = H;   %kein Flashing
    p_ris_eq_l = @(dp_ris_l)((p0 - dp_ris_l)*1e5 + g * H_VF / 2 / XSteam('v_ph',dp_ris_l,h_ris)); %Berechnung des Druckes im SR im oberen- (dampf-) unter KV
    p_ris_l = fzero(p_ris_eq_l, [1 3]);
    p_ris_v = 0;
    v_ris_v = inf;
    v_ris_l = XSteam('v_ph',p_ris_l,h_ris);
elseif h_ris > XSteam('hL_p', p_ris_b)
    H_VF = 0;   %Dampfbildung in Heizzone
    p_ris_eq_v = @(dp_ris_v)((p0 - dp_ris_v)*1e5 + g * H / 2 / XSteam('v_ph',dp_ris_v,h_ris)); %Berechnung des Druckes im SR im oberen- (dampf-) unter KV
    p_ris_v = fzero(p_ris_eq_v, [1 3]);
    p_ris_l = 0;
    v_ris_v = XSteam('v_ph',p_ris_v,h_ris);
    v_ris_l = inf;
else
    H_VF = H - ((XSteam('psat_T',XSteam('T_ph',p_ris_b,h_ris)) - p0) * 1e5 * XSteam('v_ph',p_ris_b,h_ris)) / g;
    
    p_ris_eq_v = @(dp_ris_v)((p0 - dp_ris_v)*1e5 + g * (H - H_VF) / 2 / XSteam('v_ph',dp_ris_v,h_ris)); %Berechnung des Druckes im SR im oberen- (dampf-) unter KV
    p_ris_v = fzero(p_ris_eq_v, [1 3]);
    
    p_ris_eq_l = @(dp_ris_l)(((p_ris_v - p0) * 2 + p0 - dp_ris_l)*1e5 + g * H_VF / 2 / XSteam('v_ph',dp_ris_l,h_ris)); %Berechnung des Druckes im SR im oberen- (dampf-) unter KV
    p_ris_l = fzero(p_ris_eq_l, [1 3]);
    
    v_ris_v = XSteam('v_ph',p_ris_v,h_ris);
    v_ris_l = XSteam('v_ph',p_ris_l,h_ris);
end
    

dp_acc = g * (H/v_dc - (H_VF / v_ris_l + (H-H_VF) / v_ris_v)); % Antreibende Druckdifferenz
dp_fric_dc = MartNel_ATHL(m,p_dc,h_dc,d_dc,l_dc,Zeta_dc); %Reibungsdruckverlust Fallrohr
dp_fric_ris_v = MartNel_ATHL(m,p_ris_v,h_ris,d_ris,(H-H_VF),(Zeta_ris / 2));
dp_fric_ris_l = MartNel_ATHL(m,p_ris_l,h_ris,d_ris,H_VF,(Zeta_ris / 2));
dp_fric_ris = dp_fric_ris_v + dp_fric_ris_l;
 
Sum_dp = dp_acc - (dp_fric_dc + dp_fric_ris);

end

