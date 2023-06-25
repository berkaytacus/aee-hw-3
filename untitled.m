%% Constants
clc
W_0 = 10000;
W_1 = exp(-1/4) * W_0;
S = 16;
P_A_max = 150000;
n = 1;
c = 4 * 10^(-6);
V_max = 100;
p = 1;

%% Q1
clc

V_stall = 25; %constants
thr = 1/6; 

P_A = P_A_max * thr * n;
Thrust_i = P_A / V_stall; %thrust when v = v_stall


C_L = (2 * W_0 / (p * V_stall^2 * S)) % finding CL
C_D = (2 * Thrust_i / (p * V_stall^2 * S)) % finding CD

Range = n*C_L/(c*C_D) * log(W_0 / W_1)
Endurance =  n*C_L^(3/2)/(c*C_D) * sqrt(2*p*S) * (W_1^(-1/2) - W_0^(-1/2))

V_end = sqrt(2*W_1 / (p*S*C_L));
Thrust_f = (V_end^2*S*p*C_D)/2  ;   %results
Power_f = V_end * Thrust_f;
thr_f = Power_f / P_A_max;

%% Q2

clc
CL_min = sqrt(2 * W_0 / (p * V_max^2 * S));
Thrust_1 = P_A_max / V_max;
CD_min = sqrt(2 * Thrust_1 / (p * V_max^2 * S));
Range_max = n/c * log(W_0 / W_1) * CL_1;

CL_1 = sqrt(CD0 / K);
CD_1 = 2*CD0;

V_initial2 = sqrt(2*W_1 / (p*S*CL_1))
V_final2 = sqrt(2*W_0 / (p*S*CL_1))


t_initial2 = (V_initial2^2*S*p*CD_1)/2;
t_final2 = (V_final2^2*S*p*CD_1)/2;

pa_initial2 = t_initial2 * V_initial2;
pa_final2 = V_final2 * t_final2;

thr_initial = pa_initial2 / P_A_max
thr_final = pa_final2 / P_A_max

%% Q3
clc

CL_2 = sqrt(3*CD0/K)
CD_2 = 4*CD0
Endurance_max = n*CL_2^(3/2)/(c*CD_2) * sqrt(2*p*S) * (W_1^(-1/2) - W_0^(-1/2));

V_initial2 = sqrt(2*W_1 / (p*S*CL_2))
V_final2 = sqrt(2*W_0 / (p*S*CL_2))


t_initial2 = (V_initial2^2*S*p*CD_2)/2;
t_final2 = (V_final2^2*S*p*CD_2)/2;

pa_initial2 = t_initial2 * V_initial2;
pa_final2 = V_final2 * t_final2;

thr_initial2 = pa_initial2 / P_A_max
thr_final2 = pa_final2 / P_A_max


%% Q4
clc

CD0 = 0.068;
K = 0.033;
CLmax = 2;
V = 25:100;
n_max_aero = (p * V.^2. * S * CLmax) / (2*W_0);
n_max_prop = sqrt(p.*V.^2./(2*K*W_0/S).*(P_A_max./V./W_0 - (p.*V.^2.*CD0./(W_0/S)/2)));
 
plot(V, n_max_prop)
hold on 
plot(V, n_max_aero)

title("V-n graph")
xlabel("V(m/s)")
ylabel("n")


n = n_max_aero(V==61);


%% Q5

V = (2 * P_A_max / (S*(K*CLmax^2 + CD0)))^(1/3);
n = V^2*S*CLmax/W_0/2;

R1 = 25^2 / (9.8 * sqrt(n^2-1));
Rm = V^2 / (9.8 * sqrt(n^2-1));
R3 = 90^2 / (9.8 * sqrt(n^2-1));
