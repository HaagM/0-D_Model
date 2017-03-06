
function out = Chis_ATHL( m, p, h, d, l, Zeta, r )
%m: mass flux [kg/s], p: pressure [bar], h: enthalpie [kJ/kj], d: diameter [m], l: tube length [m],
%Zeta: friction loss coeff. [-], r: wall roughness [m]
%Calculates the pressure drop for 1- and 2-phase flow using the
%correlations which are implemented in ATHLET. For 2-phase flow, 

%%
% Stoffdaten
v_l = XSteam('vL_p',p);
v_v = XSteam('vV_p',p);
v = XSteam('v_ph',p,h);
eta_l = XSteam('my_ph',p,XSteam('hL_p',p));
x = XSteam('x_ph',p,h);
% Rechengrößen
u = m * v / (pi/4 * d^2); %Strömungsgschwindigkeit
Re_l = max(u * d / eta_l / v_l , 500);
A = pi/4 * d^2;
G_A = m / A;
%
%% Darcy-Weissbach friction factor lambda
lambda_lam_l = 64 / Re_l;
lambda_turb_l = 1 / (-2 * log10(r / (3.7 * d) + (6.81 / Re_l)^0.9))^2;
lambda_l = max(lambda_lam_l , lambda_turb_l);
%
%% Two-phase friction factor 
if x == 0 
    k = lambda_l * l / (d * A^2) *0.5 * v_l;
else
    %% Chisholm model
    R = sqrt(v_v / v_l);
    if R <= 9.5
        if G_A <= 500
            B = 4.8;
        elseif (G_A > 500) && (G_A < 1900)
            B = 2400 / G_A;
        else
            B = 55 / (G_A)^0.5;
        end
    elseif (R > 9.5) && (R < 28)
        if G_A <= 600
            B = 520 / (R * G_A^0.5);
        else
            B = 21 / R;
        end
    else 
        B = 1.5e4 / (R^2 * G_A^0.5);
    end
    C_PHI = 1 + (R^2 - 1) * (B * x * (1 - x) + x^2);
    k = lambda_l * l / (d * A^2) * v_l / 2 * C_PHI;
end
out = (k + 0.5 * v * Zeta / A^2) * m^2;
end


