%% Rotating Heated Tube with Pre-Heat, Spin + Icing Dip, and Blade Wall Channels
% clear; clc;

%% GEOMETRY
D       = 4.7e-3;
t_w     = 1.32e-3;
L_A     = 0.190;      % Section A length
L_B     = 0.4318;     % Section B length
Ltot    = L_A + L_B;  %#ok<NASGU>

r_head_in  = 0.190;
r_head_out = 0.6604;   % 26"

T_inf = 268.15;        % [K] ambient
p_amb = 101325;        % [Pa]

% Fraction of inner area actually heated by flow
f_hot  = 0.56;

% A-section loss factor: fraction of temp rise that reaches B
alpha_A = 0.40;  % tune so Start-B ~ desired at ~t_spin_on

% ===== Motor -> air heating as POWER (more linear recovery) =====
Q_air_motor_max = 2;     % [W] max motor->air heating into Section A (KNOB)
tau_air_power   = 1.0;    % [s] how fast that power ramps in after ramp end (KNOB)
t_air_ramp_lin  = 150;   % [s] duration of linear ramp after ramp end (KNOB)
% ===============================================================


%% ROTOR SPEED PROFILE: pre-heat then spin
rpm_end    = 2400;     % [RPM] final speed
t_spin_on  = 650;      % [s] rotor starts spinning
t_ramp     = 60;       % [s] ramp from 0 -> rpm_end

omega_end  = 2*pi*(rpm_end/60);
omega_of_t = @(t) omega_end * max(min((t - t_spin_on)/t_ramp, 1), 0);

%% DISCRETIZATION / TIME STEP
Nx_A = 60;
Nx_B = 120;
dt   = 5e-4;                 % [s]
tEnd_max = 1000;             % [s]
T_stop_K = 273.15;           %#ok<NASGU>

% Stabilization
Hi_cap      = 0.15;
mdotE_floor = 1e-4;          % [kg/s]

%% CONSTANT HEATED INLET (ALWAYS ON)
mDot_in_of_t = @(t) 0.003;           % [kg/s]

T_in_final = 40 + 273.15;   % [K]
% smooth heater warm-up over first 600 s
T_in_of_t  = @(t) T_inf + (T_in_final - T_inf) * min(max(t/600,0),1);

%% ===== MOTOR THERMAL MODEL (FAST RAMP + CONTINUED RISE AT SPEED) =====
Tmotor_init_C = 26;                  % [C]
Tmotor_min_C  = Tmotor_init_C;

% 1) Steady-state rise due to RPM-dependent losses
dTmotor_ss_rpm_C = 65;               % [C] rise due to RPM at full speed  (KNOB)
motor_rise_exp   = 2.0;              % omega exponent (KNOB)

% 2) Extra "soak" rise that ramps in AFTER reaching full speed
%    (captures motor/ESC/hub soaking even at constant omega)
dTmotor_ss_soak_C = 35;              % [C] additional rise after ramp (KNOB)
tau_soak          = 180;             % [s] soak timescale (KNOB)

% 3) Make motor respond faster during ramp than later (2 time constants)
tau_motor_ramp = 5;                 % [s] fast during ramp (KNOB)
tau_motor_run  = 40;                % [s] slower after ramp (KNOB)

% Helpful timing
t_ramp_end = t_spin_on + t_ramp;

% RPM fraction
omega_frac = @(t) min(max(omega_of_t(t)/omega_end, 0), 1);

% RPM-based steady state
Tss_rpm_K = @(t) (Tmotor_init_C + 273.15) + ...
                 dTmotor_ss_rpm_C * (omega_frac(t).^motor_rise_exp);

% Soak contribution: only after ramp ends, rises exponentially
soak_gain = @(t) (t >= t_ramp_end) .* (1 - exp(-(t - t_ramp_end)/tau_soak));
Tss_soak_K = @(t) dTmotor_ss_soak_C * soak_gain(t);

% Total steady-state target
Tmotor_ss_K = @(t) Tss_rpm_K(t) + Tss_soak_K(t);

% Effective tau: fast during ramp, slow after
tau_motor_eff = @(t) (t < t_ramp_end)*tau_motor_ramp + (t >= t_ramp_end)*tau_motor_run;
% ===================================================================

%% INLET HEATING NEAR MOTOR (Section A)
G_contact       = 100.0;     % [W/K] motor->inner wall (cell 1) (TUNING KNOB)
h_motor_air     = 18.0;      % [W/m^2-K] motor->air exchange (TUNING KNOB)
A_motor_air     = 700e-6;    % [m^2]
Nsrc_A          = 12;        % first N cells of A receive motor->air heating (TUNING KNOB)
dT_air_src_cap  = 2.0;       % [K/pass] cap for stability (leave)
% ===== EXTRA distributed heating in B (to linearize recovery) =====
Nsrc_B          = 25;      % [-] first N cells of B get extra heat (try 15–40)
Q_B_max         = 35;      % [W] max extra power into B air (try 10–60)
t_B_ramp_lin    = 260;     % [s] linear ramp after ramp end (try 120–400)
dT_B_src_cap    = 2.0;     % [K/pass] stability cap (leave)
% ================================================================
G_hub_B       = 0.0;     % [W/K] hub->B inner wall conduction strength (TUNING KNOB)
Nhub_B        = 20;       % [-] number of B cells near root that receive hub heat
t_hub_on      = t_spin_on + t_ramp;  % start after ramp completes (turnaround point)
tau_hub_rise  = 50;       % [s] how fast this conduction ramps in (TUNING KNOB)

Qhub_cap_W = 20;     % [W] HARD CAP total hub power into B (PHYSICAL LIMIT)
hub_allow_cooling = false;  % set true only if you want hub to cool blade

% hub_gain_of_t = @(t) (t>=t_hub_on) .* (1 - exp(-(t - t_hub_on)/tau_hub_rise));
t_hub_ramp_lin = 100;  % [s] linear ramp duration (KNOB)
hub_gain_of_t  = @(t) (t<=t_hub_on)*0 + ...
                      (t> t_hub_on)*min(max((t - t_hub_on)/max(t_hub_ramp_lin,1e-6),0),1);

%% CONSTANTS / HELPERS
R_air   = 287.058;
cp_air  = 1010;
Pr_air  = 0.70;
mu0     = 18.27e-6;
T0      = 291.15;
S       = 120;
mu_of_T = @(T) mu0.*((T0+S)./(T+S)).*(T./T0).^(3/2);

rho_cu = 8960;
cp_cu  = 385;

k_cu_radial = 500;  % [W/m-K] radial
k_cu_axial  = 800;  % [W/m-K] axial

haaland    = @(Re,epsD)(-1.8*log10((epsD/3.7).^1.11 + 6.9./max(Re,100))).^(-2);
k_air_amb  = 0.024;
nu_amb     = 1.4e-5;
Pr_amb     = 0.72;
hout_of    = @(om,rrot)(0.62*(om.*rrot.^2/nu_amb).^0.5 .* Pr_amb.^0.37).*k_air_amb./(2*rrot);

ri  = D/2;
ro  = ri + t_w;
Acs = pi*(D^2)/4;
P_in  = pi*D;
P_out = 2*pi*ro;
dx_A = L_A/Nx_A;
dx_B = L_B/Nx_B;

r_mid = sqrt(ri*ro);
Aax_in  = pi*(r_mid^2 - ri^2);
Aax_out = pi*(ro^2    - r_mid^2);

V_in_A  = Aax_in  * dx_A;
V_out_A = Aax_out * dx_A;
V_in_B  = Aax_in  * dx_B;
V_out_B = Aax_out * dx_B;

Ci_A = rho_cu*cp_cu*V_in_A;
Co_A = rho_cu*cp_cu*V_out_A;
Ci_B = rho_cu*cp_cu*V_in_B;
Co_B = rho_cu*cp_cu*V_out_B;

% Extra thermal mass for Section B
mass_scale_B = 10;
Ci_B = mass_scale_B * Ci_B;
Co_B = mass_scale_B * Co_B;

Gcond_A = (2*pi*k_cu_radial*dx_A)/log(ro/ri);
Gcond_B = (2*pi*k_cu_radial*dx_B)/log(ro/ri);

secDiff_over_dx = @(T,dx)[T(2)-T(1);
                          T(3:end)-2*T(2:end-1)+T(1:end-2);
                          T(end-1)-T(end)] ./ dx;

tiny      = 1e-12;
mDot_floor= 1e-6;
eps_rel   = 0.0;

%% ICING PARAMETERS (Section B outer wall)
icing_on   = true;

LWC        = 0.05e-3;     % [kg/m^3]
T_drop     = 263.15;      % [K]
cp_w       = 4180;
L_f        = 3.34e5;

eta_const = 0.10;

% --- NEW: taper collection efficiency vs radius (reduces tip cooling) ---
eta_tip_factor = 0.2;   % [-] tip collects 55% of root (KNOB)  (0.4–0.8)
eta_shape      = 1;    % [-] shape exponent (KNOB)           (1–3)

r_B_vec    = linspace(r_head_in, r_head_out, Nx_B)';

r_norm = (r_B_vec - r_head_in) / (r_head_out - r_head_in);  % 0 at root, 1 at tip
eta_col_B = eta_const * (1 - (1 - eta_tip_factor)*(r_norm.^eta_shape));
% -----------------------------------------------------------------------


r_B_vec    = linspace(r_head_in, r_head_out, Nx_B)';

t_ice_on   = t_spin_on + t_ramp;
t_ice_off  = 840;

dT_ice_cap = 0.12;

%% EXTERNAL CONVECTION
h_nat      = 8.0;   % [W/m^2-K]
hout_factor= 10;    % multiplier on rotating correlation (TUNING KNOB)

% ===== RAMP-ONLY forced convection boost (deepens dip) =====
hout_ramp_boost = 2.5;   % [–] >1 deepens dip only during ramp 
% ==========================================================

% ===== Motor->air pickup lag knobs (sharpens recovery) =====
tau_air_pickup_fast = 6.0;    % [s] fast “snap” right after ramp end 
tau_air_pickup_slow = 25.0;   % [s] slower tail 
w_air_fast          = 0.9;    % [–] fraction of fast component 
% ==========================================================

%% “CHANNEL” INDICES IN B
iB_start = round(0.1 * Nx_B);
iB_half  = round(0.5 * Nx_B);
iB_end   = Nx_B;

%% INITIAL CONDITIONS
Tf_A  = T_inf*ones(Nx_A,1);
Twi_A = T_inf*ones(Nx_A,1);
Two_A = T_inf*ones(Nx_A,1);
Tf_B  = T_inf*ones(Nx_B,1);
Twi_B = T_inf*ones(Nx_B,1);
Two_B = T_inf*ones(Nx_B,1);

% Motor temperature state [K]
Tmotor_state = Tmotor_init_C + 273.15;

t  = 0;
n  = 0;

t_hist          = [];
Twi_start_hist  = [];
Twi_half_hist   = [];
Twi_end_hist    = [];
ice_active_hist = [];
houtB_hist      = [];
Tmotor_hist     = []; 
TfB_start_hist = [];
TfB_half_hist  = [];
TfB_end_hist   = [];
TwoB_start_hist = [];
TwoB_half_hist  = [];
TwoB_end_hist   = [];

% debug

%% MAIN LOOP
while (t < tEnd_max) && isfinite(min(Two_B)) && isreal(min(Two_B))
    n = n+1; 
    t = t + dt;
    om = omega_of_t(t);

    % Icing active
    ice_active = icing_on && (t >= t_ice_on) && (t <= t_ice_off);
    ice_active_hist(end+1,1) = ice_active;

    %% Motor temperature state update 
    Tss = Tmotor_ss_K(t);
    Tmotor_state = Tmotor_state + dt*(Tss - Tmotor_state)/tau_motor_eff(t);
    Tmotor_state = max(Tmotor_state, Tmotor_min_C + 273.15);
    Tmotor_hist(end+1,1) = Tmotor_state;

    T_inlet_base = T_in_of_t(t);
    T_inlet = T_inlet_base;


    T_inlet_mix = T_inlet;

    % Mean gas properties
    Tmean_gas = mean([Tf_A; Tf_B]);
    Tm        = max(Tmean_gas, 200);
    rho_m     = p_amb/(R_air*Tm); 

    mDot_in = max(mDot_in_of_t(t), 0);
    mDot    = max(mDot_in, mDot_floor);
    mdotE   = max(mDot, mdotE_floor);

    % Inner h (A)
    muA  = mu_of_T(Tf_A);
    rhoA = p_amb./(R_air*Tf_A);
    uA   = mDot./(rhoA*Acs+tiny);
    ReA  = max(rhoA.*uA*D./(muA+tiny),1);
    fA   = haaland(ReA,eps_rel);
    kA   = muA.*cp_air./Pr_air;
    NuA  = (fA./8).*(ReA-1000).*Pr_air./(1+12.7*sqrt(fA./8).*(Pr_air.^(2/3)-1));
    NuA(ReA<3000) = 3.66;
    hA = 2.0 * NuA.*kA/max(D,tiny);

    % Inner h (B)
    muB  = mu_of_T(Tf_B);
    rhoB = p_amb./(R_air*Tf_B);
    uB   = mDot./(rhoB*Acs+tiny);
    ReB  = max(rhoB.*uB*D./(muB+tiny),1);
    fB   = haaland(ReB,eps_rel);
    kB   = muB.*cp_air./Pr_air;
    NuB  = (fB./8).*(ReB-1000).*Pr_air./(1+12.7*sqrt(fB./8).*(Pr_air.^(2/3)-1));
    NuB(ReB<3000) = 3.66;
    hB = 2.0 * NuB.*kB/max(D,tiny);

   %% Outer Convection (spin active after t_spin_on)
hout_A = h_nat;
hout_B = h_nat;

if t >= t_spin_on
    houtA_forced = hout_factor * hout_of(om, r_head_in);
    houtB_forced = hout_factor * hout_of(om, r_head_out);

    blend = max(0, min(1, (t - t_spin_on)/t_ramp));
    
    if (blend > 0) && (blend < 1)
        houtA_forced = houtA_forced * hout_ramp_boost;
        houtB_forced = houtB_forced * hout_ramp_boost;
    end

    hout_A = (1 - blend)*h_nat + blend*houtA_forced;
    hout_B = (1 - blend)*h_nat + blend*houtB_forced;
end

houtB_hist(end+1,1) = hout_B;


    %% Fluid march A
    TfA_new = Tf_A;
    Tf_up   = T_inlet_mix;

    % ===== Motor -> air pickup DELAY (physical mixing lag) =====
% ===== Motor -> air heating: linear available power + speed-dependent UA cap =====
% Goal: make recovery follow linear Q_avail (not exponential motor state),
% while still respecting physics via UA*(Tmotor - Tair).

% ---- KNOBS ----
Q_air_motor_max   = 40;    % [W] max motor->air power available (your original)
t_air_ramp_lin    = 100;   % [s] duration of linear ramp after ramp end

UA0   = h_motor_air * A_motor_air;  % [W/K] baseline UA from your h & area
UA_speed_mult = 10.0;       % [-] UA at full speed relative to baseline 
t_UA_rise      = 4.0;      % [s] how fast UA ramps up after ramp end 

Q_cap_hard = 35;           % [W] hard physical ceiling 
% --------------------------------------------

% (1) Linear "available" power ramp after ramp end 
if t <= t_ramp_end
    Q_avail = 0;
else
    dt_since = t - t_ramp_end;
    Q_avail  = Q_air_motor_max * min(max(dt_since/max(t_air_ramp_lin,1e-6),0),1); % linear 0->max
end

% (2) UA cap increases after ramp end (forced mixing / better coupling at speed)
if t <= t_ramp_end
    UA_gain = 1.0;
else
    dt_since = t - t_ramp_end;
    UA_gain  = 1.0 + (UA_speed_mult - 1.0) * (1 - exp(-dt_since/max(t_UA_rise,1e-6)));
end

UA_eff = UA0 * UA_gain;    % [W/K]

% (3) Cap from exchanger physics
DeltaT_ma = max(Tmotor_state - TfA_new(1), 0);   % [K] use local A air temp near motor (better than inlet)
Q_cap_UA  = UA_eff * DeltaT_ma;                  % [W]
Q_cap_UA  = min(Q_cap_UA, Q_cap_hard);           % [W] hard ceiling

% (4) Final applied motor->air power
Qdot_ma_total = min(Q_avail, Q_cap_UA);

% distribute to the first Nsrc_A cells
Qdot_ma_cell  = Qdot_ma_total / max(Nsrc_A,1);
% ========================================================================


    for i=1:Nx_A
        Hi = hA(i)*P_in*dx_A / (mdotE*cp_air + 1e-12);
        Hi = min(Hi, Hi_cap);

        Ssrc = 0;
        if i <= Nsrc_A
            S_uncapped = Qdot_ma_cell / (mdotE*cp_air + 1e-12);
            Ssrc = sign(S_uncapped) * min(abs(S_uncapped), dT_air_src_cap);
        end

        TfA_new(i) = Tf_up + Hi*(Twi_A(i) - Tf_up) + Ssrc;
        Tf_up = TfA_new(i);
    end

    Tin_B_raw = TfA_new(end);
    Tin_B = T_inf + alpha_A * (Tin_B_raw - T_inf);

    %% Wall A
    Qin_i_A  = f_hot*hA.*P_in.*(TfA_new - Twi_A)*dx_A;
    Qout_o_A =        hout_A.*P_out.*(Two_A   - T_inf)*dx_A;
    Qrad_A   = Gcond_A.*(Two_A - Twi_A);

    % motor -> wall conduction (use motor_state, physical)
    Q_contact_motor = G_contact*(Tmotor_state - Twi_A(1));
    Qin_i_A(1) = Qin_i_A(1) + Q_contact_motor;

    Qax_in_A  = k_cu_axial*Aax_in  * secDiff_over_dx(Twi_A,dx_A);
    Qax_out_A = k_cu_axial*Aax_out * secDiff_over_dx(Two_A,dx_A);

    Twi_A = Twi_A + dt*((Qin_i_A + Qrad_A + Qax_in_A) ./ Ci_A);
    Two_A = Two_A + dt*((-Qrad_A - Qout_o_A + Qax_out_A) ./ Co_A);
    Tf_A  = TfA_new;

    %% Fluid march B
    TfB_new = Tf_B;
    Tf_up   = Tin_B;
% ---- B distributed power schedule (linear after ramp end) ----
if t <= t_ramp_end
    Q_B_avail = 0;
else
    dt_since  = t - t_ramp_end;
    Q_B_avail = Q_B_max * min(max(dt_since/max(t_B_ramp_lin,1e-6),0),1);
end
Q_B_cell = Q_B_avail / max(Nsrc_B,1);
% -------------------------------------------------------------
   for i=1:Nx_B
    Hi = hB(i)*P_in*dx_B / (mdotE*cp_air + 1e-12);
    Hi = min(Hi, Hi_cap);

    % Extra distributed heating into B air (first Nsrc_B cells)
    SsrcB = 0;
    if i <= Nsrc_B
        S_uncapped = Q_B_cell / (mdotE*cp_air + 1e-12);
        SsrcB = sign(S_uncapped) * min(abs(S_uncapped), dT_B_src_cap);
    end

    TfB_new(i) = Tf_up + Hi*(Twi_B(i) - Tf_up) + SsrcB;
    Tf_up = TfB_new(i);
end


    %% Wall B base terms
    Qin_i_B  = f_hot*hB.*P_in.*(TfB_new - Twi_B)*dx_B;
    Qout_o_B =        hout_B.*P_out.*(Two_B   - T_inf)*dx_B;
    Qrad_B   = Gcond_B.*(Two_B - Twi_B);

    Qax_in_B  = k_cu_axial*Aax_in  * secDiff_over_dx(Twi_B,dx_B);
    Qax_out_B = k_cu_axial*Aax_out * secDiff_over_dx(Two_B,dx_B);

    %% Icing load on B
    Qice_B = zeros(Nx_B,1);
    if ice_active
        V_rel_B        = om .* r_B_vec;
        m_dot_wpp_B    = LWC .* V_rel_B .* eta_col_B;
        A_cell_B       = P_out * dx_B;
        m_dot_w_cell_B = m_dot_wpp_B * A_cell_B;
        dH_ice         = cp_w*(273.15 - T_drop) + L_f;

        mask_ice = (Two_B < 273.15);
        Qice_B(mask_ice) = m_dot_w_cell_B(mask_ice) * dH_ice;

        Qcap_cell = Co_B * dT_ice_cap / dt;
        Qice_B = min(Qice_B, Qcap_cell);
    end
%% ===== hub conduction heating to B inner wall near root (CAPPED, PHYSICAL) =====
Qhub_B = zeros(Nx_B,1);
g_hub  = hub_gain_of_t(t);

Nuse = min(Nhub_B, Nx_B);

dT_hub = (Tmotor_state - Twi_B(1));     % [K]
Qhub_total = g_hub * G_hub_B * dT_hub;  % [W]

% Physical behavior:
%  - cap magnitude (finite conduction path / contact area)
%  - optionally forbid cooling
if ~hub_allow_cooling
    Qhub_total = max(Qhub_total, 0);
end
Qhub_total = max(min(Qhub_total, Qhub_cap_W), -Qhub_cap_W);

Qhub_B(1:Nuse) = Qhub_total / Nuse;
% ============================================================================


    %% Update B walls
    Twi_B = Twi_B + dt*((Qin_i_B + Qrad_B + Qax_in_B + Qhub_B) ./ Ci_B);
    Two_B = Two_B + dt*((-Qrad_B - Qout_o_B - Qice_B + Qax_out_B) ./ Co_B);
    Tf_B  = TfB_new;

    %% Logs
    t_hist(end+1,1)         = t;
    Twi_start_hist(end+1,1) = Twi_B(iB_start);
    Twi_half_hist(end+1,1)  = Twi_B(iB_half);
    Twi_end_hist(end+1,1)   = Twi_B(iB_end);
    TfB_start_hist(end+1,1) = Tf_B(iB_start);
    TfB_half_hist(end+1,1)  = Tf_B(iB_half);
    TfB_end_hist(end+1,1)   = Tf_B(iB_end);
    TwoB_start_hist(end+1,1) = Two_B(iB_start);
    TwoB_half_hist(end+1,1)  = Two_B(iB_half);
    TwoB_end_hist(end+1,1)   = Two_B(iB_end);


end

%% PLOTS
figure;
plot(t_hist, Twi_start_hist - 273.15, 'LineWidth', 1.5); hold on;
plot(t_hist, Twi_half_hist  - 273.15, 'LineWidth', 1.5);
plot(t_hist, Twi_end_hist   - 273.15, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Inner Wall T [°C]');
legend('Start B','Half B','End B','Location','best');
title('Simulated Blade Inner-Wall Temperature vs Time');

figure;
plot(t_hist, TfB_start_hist - 273.15, 'LineWidth', 1.5); hold on;
plot(t_hist, TfB_half_hist  - 273.15, 'LineWidth', 1.5);
plot(t_hist, TfB_end_hist   - 273.15, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Temperature [°C]');
xlim([0,1000]);
ylim([-10,25]);
legend('Start B','Half B','End B','Location','best');
title('Simulated Section B Air Temperature vs Time');

figure;
plot(t_hist, TwoB_start_hist - 273.15, 'LineWidth', 1.5); hold on;
plot(t_hist, TwoB_half_hist  - 273.15, 'LineWidth', 1.5);
plot(t_hist, TwoB_end_hist   - 273.15, 'LineWidth', 1.5);
grid on;

xlabel('Time [s]');
ylabel('Outer Wall Temperature [°C]');
legend('Start B','Half B','End B','Location','best');
title('Section B External (Outer Wall) Temperature vs Time');

figure;
stairs(t_hist, (t_hist>=t_ice_on) & (t_hist<=t_ice_off), 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Icing Active (0/1)');
title('Icing Window');

figure;
plot(t_hist, houtB_hist, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('h_{out,B} [W/m^2-K]');
title('Outer Convection Coefficient vs Time');

figure;
plot(t_hist, Tmotor_hist - 273.15, 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Motor Temp [°C]');
title('Motor Temperature State (debug)');
