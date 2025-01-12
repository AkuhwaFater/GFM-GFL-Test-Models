%
% Parameters file for SPS model: HVDC_MMC.slx
%
%
GFMGFLSelect = 0; %1 for GFM 0 for GFL
Fnom= 50;                      % Nominal system frequency (Hz)
Pnom= 1000e6;                  % Converter 3-phase rated power (MVA)
Vnom_prim= 400e3;              % Nominal primary voltage (V)
Vnom_sec= 333e3;               % Nominal secondary voltage (V)
Nb_PM=36;                      % Number of power module per arm
Vnom_dc= 640e3;                % DC nominal voltage (V)

C_PM= 1.758e-3; % Power module capacitor (F)
% Energy in kJ/MVA
W_kJ_MVA= 0.5 * C_PM * (Vnom_dc/Nb_PM)^2 * Nb_PM * 6 / (Pnom/1e6)/1e3;
Vc0_PM=0;                     % Capacitors initial voltage (V)

%% GFM Converter Parameters
gridInverter.apparentPower   = Pnom;  % kVA, Apparent power
gridInverter.frequency       = Fnom;   % Hz, Grid frequency
gridInverter.DCVoltage       = Vnom_dc; % V, DC bus voltage
gridInverter.lineRMSVoltage  = Vnom_sec;  % V, Line RMS voltage at the point of interconnection
gridInverter.measurementSampleTime = 40e-6; % s, Power measurement time constant

% Estimating the base values
base.power = gridInverter.apparentPower; % kVA
base.frequency = gridInverter.frequency; % Hz
base.lineVoltage = gridInverter.lineRMSVoltage; % V
base.basePhasePower = base.power*1e3/3; % kVA
base.basePhaseVoltage = base.lineVoltage/sqrt(3); % V
base.voltage = base.basePhaseVoltage*sqrt(2); % V
base.basePhaseCurrent = base.basePhasePower/base.basePhaseVoltage; % A
base.current = base.basePhaseCurrent*sqrt(2); % A

base.impedance = base.basePhaseVoltage/base.basePhaseCurrent; % Ohm
base.inductance = base.impedance/(2*pi*base.frequency); % H
base.capacitance = 1/(base.impedance*2*pi*base.frequency); % F

% grid.gridVoltageLL  = Vnom_prim; % V, Grid line RMS voltage
% grid.gridResistance = 3;     % Ohm, Grid source resistance
% grid.gridInductance = 0.05;  % H, Grid source inductance

%Droop Active Power Control Parameters
gridInverter.droopControl.freqSlopeMp = 0.01; % pu, Hz/W, Power droop value

gridInverter.droopControl.lpfTimeConst = 0.015; % s, Low pass filter time constant

% Lead-lag parameter for three phase power measurement
gridInverter.droopControl.T2 = 0.006; % s, Denominator time constant
gridInverter.droopControl.T1 = 0.005; % s, Numerator time constant
gridInverter.droopControl.sampleTime = 40e-6; % s, Sampling time

gridInverter.freqMeasTimeConst = 150e-3; % s, Frequency measurement time constant

%Virtual Synchrnous Machine (VSM) Active Power Control Parameters 
gridInverter.vsm.inertiaConstant = 1; % s, Mechanical time constant
gridInverter.vsm.dampingCoefficent = 1.056; % pu/Hz, Damping coefficient
gridInverter.vsm.freqDroop  = 10; % pu, W(pu)/Hz(pu) VSM Frequency droop
gridInverter.vsm.PmeasTimeConst = 1e-3; % s, Power measurement filter time constant
gridInverter.vsm.maxDampingPower = 0.7; % pu
gridInverter.vsm.minDampingPower = -0.6; % pu
gridInverter.vsm.samplingTime = gridInverter.droopControl.sampleTime; % s
gridInverter.vsm.dampingPowerOption = 'Constant Frequency Reference'; % Selecting the damping frequency option

%Reactive Power Droop Control
gridInverter.Qcontrol.voltageDroop = 0.3; % pu, V/VAR
gridInverter.Qcontrol.QmeasTimeConst = 1e-3; % s, Power measurement time constant
gridInverter.Qcontrol.voltageReference = 1.0; % pu
gridInverter.Qcontrol.lowVoltageSupportGain = 1.5; % pu.A/V, it adds more reactive current (Iq), during the low voltage condition

%Virtual Impedance Parameters
gridInverter.currentLimit.virImpResistanceCoeff = 0.1875; % pu, Resistance
gridInverter.currentLimit.virImpXbyR = 13.2; % X/R ratio
gridInverter.currentLimit.viCurrentLimit = 1.2; % A, Maximum current
gridInverter.currentLimit.viFilterTimeConst = 1e-3; % s, Filter time constant

%Current Limiting Parameters
gridInverter.currentLimit.maxSaturationCurrent = 1.2; % pu
gridInverter.currentLimit.maxSaturationDelay = 1e-3; % s

%Virtual Impedance and Current Limiting Parameters
gridInverter.currentLimit.satCurrentRunTime = 1e-3; % s % Current saturation run time

%Current Controller Parameters
gridInverter.controller.CurrentControlSampleTime = 40e-6; % s, Controller sampling time
gridInverter.controller.ctControllerKp = 1.5; % Proportional gain
gridInverter.controller.ctControllerKi = 10; % Integral gain

%Voltage Controller Parameters
gridInverter.controller.VoltageControlSampleTime = 40e-6; % s, Controller sampling time
gridInverter.controller.voltControllerKp = 0.3; % Proportional gain
gridInverter.controller.voltControllerKi = 180; % Integral gain

gridInverter.controller.voltageMaxId = 1.2; % Id controller saturation maximum limit
gridInverter.controller.voltageMinId = -1.2; % Id controller saturation minimum limit
gridInverter.controller.voltageMaxIq = 1.2; % Iq controller saturation maximum limit
gridInverter.controller.voltageMinIq = -1.2; % Iq controller saturation minimum limit
gridInverter.controller.voltMeasTimeConst = 5e-3; % Voltage measurement low pass filter time constant

%GFM Filter Inductor Design
gridInverter.ratedrmsCurrent = gridInverter.apparentPower*1e3/(sqrt(3)*gridInverter.lineRMSVoltage); % A, Filter rated current
gridInverter.L = (0.1*gridInverter.lineRMSVoltage/(gridInverter.ratedrmsCurrent*gridInverter.frequency*2*pi*sqrt(3))); % H, Filter inductance
gridInverter.lineResistance = 0.1; % Ohm, Filter resistance
%%
% Sequencer timing:
Tbrk1_On=0.1;                 % Closing time of breaker 1 (converter energizing)
Tbrk2_On=1.0;                 % Closing time (s) of breaker 2 (across start-up resistor)
%
Tdeblock=1.5;                 % Converter de-block time (s)
Ton_VDCreg=1.5;               % VDC regulator turn-on time (s) - VDC Regulation
Tramping_Vdc_ref=2;           % Start time Vdc_ref ramping to nominal (s)
Slope_Vdc_ref=Vnom_dc/5;      % Sloge Vdc_ref ramping (V/s)
%
Ton_PQreg=4;                  % Preg & Qreg regulators turn-on time (s) - PQ regulation
Tramping_Pref=Ton_PQreg+0.2;  % Start time Pref ramping(s)
Slope_Pref=0.5;               % Sloge Pref ramping (V/s)
Tramping_Qref=Ton_PQreg+3.5;  % Start time Pref ramping(s)
Slope_Qref=0.5;               % Sloge Pref ramping (V/s)
%
Ton_Converter2=4;             % Converter 2 equivalent switched-on time (s)
%%
Tfault= 9999;             % DC Fault timing (s)
Rfault=1;                 % DC Fault resistance (Ohms)
%
%%
% PWM Output pulses selector
pp=0;
for p=1:2:72
    pp=pp+1;
    SelectPulses1(p)=pp;
    SelectPulses1(p+1)=pp+36;
end
%
Ts_Power= 20e-6;    % SPS Simulation time step(s)
Ts_Control=40e-6;   % Control system time step (s)
Ts=Ts_Control;
%
% Transformer impedance
Lxfo= 0.12;         % Total Leakage inductance (pu)
Rxfo= 0.003;        % Total winding resistance (pu)
%
Zbase= Vnom_sec^2/Pnom;
%
Larm_pu=0.15;
Rarm_pu=Larm_pu/100;
%
Zbase= Vnom_sec^2/Pnom;
Larm=Larm_pu*(Zbase/(2*pi*Fnom));
Rarm=Rarm_pu*Zbase;
w=2*pi*Fnom;
wc2=(2*w)^2;
Cfilter=1/(Larm*wc2);      % Capacitor value for 2th harmonic filter(F)
Rfilter=1/(Cfilter*w)*30;  % Resistance value for 2th harmonic filter (Ohm)
Topen_Filter=1e6;          % Breaker opening time for second-harmonic filters (s)
%
% *****************************************************************
%                         CONTROL PARAMETERS
% *****************************************************************
%
% Modulator Parameters
Fc=Fnom*3.37;        % Carriers frequency (Hz)
%
% dq and Vdc measurement filter cut-off frequency:
Fn_filter=1000;
Zeta_filter=1;
%
% Active power regulator (Preg)
Kp_Preg= 0.5/3;                % Proportional gain
Ki_Preg= 1.0;                  % Integral gain
Limits_Preg = [ 1.2, 0.8 ] ;   % Output (Vdc_ref) Upper/Lower limits (pu)

%
% Reactive power regulator (Qreg)
Kp_Qreg= 0.5/3;                % Proportional gain
Ki_Qreg= 1.0;                  % Integral gain
Limits_Qreg = [ 0.25, -0.25 ]; % Output (Iq_ref) Upper/Lower limit (pu)
%
% VDC regulator (VDCreg)
Kp_VDCreg=4;                   % Proportional gain
Ki_VDCreg=100;                 % Integral gain
Limits_VDCreg= [ 2.0  -2.0];   % Output Idref [Upper Lower] limits (pu)
%
% Current regulator (Ireg)
Kp_Ireg= 0.6;                  % Proportional gain
Ki_Ireg= 6;                    % Integral gain
Limits_Ireg= [ 2.0  -2.0];     % Output Vdq_conv [Upper Lower] limits (pu)
% Feedforward coefficients:
Lff=Larm_pu/2;
Rff= Rarm_pu/2;
%
% ******************************
% Power system parameters
% ******************************
%
Psc= Pnom*20;     % Short circuit power (MVA)
X_R= 7;           % X/R ratio
P_Ld1= Psc/30;   % load (primary bus) (MW)
R_startup= 400;   % Startup resistance (Ohm)
%
% Cable data
R_cable = 0.5;      % ohm
L_cable= 15e-3;   % (H)
%
% Grounding reference (series RC)
Rg= 100;              % (Ohms)
Cg= 50e-9;            % (F)



