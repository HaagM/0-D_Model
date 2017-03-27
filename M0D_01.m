clear;
clc;

global T_ini p0 P_in H r d_dc l_dc Zeta_dc T_dc p_dc d_ris l_ris Zeta_ris g 

%% input
%
T_ini = 20;
p0 = 1;
P_in = 300;
H = 7.348;
r = 4e-5; % Wand Rauigkeit
% downcomer
d_dc = 41.8e-3;
l_dc = 11.436;
Zeta_dc = 7.202;
T_dc = T_ini;
p_dc = p0;
% riser
d_ris = 53e-3;
l_ris = H;
Zeta_ris = 1.5;

%% Konstanten
g = 9.81;

Sum_dp_eq = @press_eq;

%% Path
path(path, 'Druckverlust');

m_max = 2;
interv = 0.04;
n_sol = 0;
Sum_dp = zeros(m_max / interv + 1 , 4);
M = zeros(m_max / interv + 1 , 1);
counter = 0;
for m = 0 : interv : m_max %Sucht in welchen Interval sich eine 0-Stelle befindet für Sum_dp 
    counter = counter + 1;
    M(counter) = m; 
    [Sum_dp(counter,1), Sum_dp(counter,2), Sum_dp(counter,3), Sum_dp(counter,4), Sum_dp(counter,5)] = Sum_dp_eq(m); %Berechnen der Druckbilanz an den Intervallgrenzen
    if m == 0
        continue
    end
    if (Sum_dp(counter , 1) * Sum_dp(counter - 1 , 1)) > 0 || isnan(Sum_dp(counter - 1 , 1)) || isnan(Sum_dp(counter , 1)) %Prüft ob sich das Vorzeichen über ein Interval ändert
        continue
    elseif Sum_dp(counter , 1) == 0;
        n_sol = n_sol +1;
        m_sol(n_sol) = m;
    else  
        n_sol = n_sol + 1;
        m_sol(n_sol) = fzero(Sum_dp_eq,[m-interv m]); % 0-Stellensuche für Druckgleichgewicht im durchsuchten Interval
    end
end
 
plot(M,Sum_dp(:,1)*1e-5 , M,Sum_dp(:,2)*1e-5 , M,Sum_dp(:,3)*1e-5 , M,Sum_dp(:,4)*1e-5 , M,Sum_dp(:,5)/10)
legend('Sum dp','dp acc','dp fric','dp fric ris','HVF')
grid on
