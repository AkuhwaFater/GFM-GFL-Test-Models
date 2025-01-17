%
% Parameters file for SPS model: HVDC_MMC.slx
%
%load sound                     % Sound file used by the Mode of Operation Panel
%
%% Model Configurations
%Model = 'FaultAnalysis';
%FltProgBlock = Simulink.findBlocks(Model,'Name','Fault Programming');
Fnom= 50;

%Plant Simulation config
PlantSim.SamplesPerCycle = 512;
PlantSim.Freq = Fnom;
%PltSim.Ts = 1/PltSim.Freq/PltSim.samplesPerCycle;
PlantSim.Ts = 2e-5;
PlantSim.Source_S_Type = "IBR_DYg";

%Remote source parameters
Source_R.SCL    = 25460e6*500; %VA
Source_R.Vbase  = 735e3; %Vrmp ph-ph
Source_R.XoverR = 10;
Source_R.Freq   = Fnom; %Hz
Source_R.Vmag   = 735e3*0.99; %Vrms
Source_R.Vphase = 6.63-26.02; %degrees

Source_R.L1     = Source_R.Vbase^2/Source_R.SCL *1/(2*pi*Source_R.Freq); %H
Source_R.R1     = 2*pi*Source_R.Freq*Source_R.L1/Source_R.XoverR; %ohm
Source_R.Z1     = Source_R.R1 + 1i*Source_R.L1*2*pi*Source_R.Freq; %ohm
Source_R.Z2     = Source_R.Z1; %ohm
Source_R.Z0Z1Ratio = 3;
Source_R.L0     = Source_R.Z0Z1Ratio * Source_R.L1; %H
Source_R.R0     = Source_R.Z0Z1Ratio * Source_R.R1; %ohm
Source_R.Z0     = Source_R.R0 + 1i*Source_R.L0*2*pi*Source_R.Freq; %ohm

%Local source parameters

%Voltage Source Parameters
Source_S.SCL    = 31800e6*500; %VA
Source_S.Vbase  = 735e3; %Vrmp ph-ph
Source_S.XoverR = 10;
Source_S.Freq   = Fnom; %Hz
Source_S.Vmag   = 735e3*1.00701; %Vrms
Source_S.Vphase =  6.6-16; %degrees

Source_S.L1     = Source_S.Vbase^2/Source_S.SCL *1/(2*pi*Source_S.Freq); %H
Source_S.R1     = 2*pi*Source_S.Freq*Source_S.L1/Source_S.XoverR; %ohm
Source_S.Z1     = Source_S.R1 + 1i*Source_S.L1*2*pi*Source_S.Freq; %ohm
Source_S.Z2     = Source_S.Z1; %ohm
Source_S.Z0Z1Ratio = 3;
Source_S.L0     = Source_S.Z0Z1Ratio * Source_S.L1; %H
Source_S.R0     = Source_S.Z0Z1Ratio * Source_S.R1; %ohm
Source_S.Z0     = Source_S.R0 + 1i*Source_S.L0*2*pi*Source_S.Freq; %ohm

%IBR Parameters
IBR.Freq = Fnom;


%Local side transformer impedance
Trf.Pnom       = 1000e6; %VA
Trf.Freq       = Fnom; %Hz
Trf.Vprim      = 400e3; %V
Trf.Vsecondary = 333e3; %V
Trf.R1         = 0.003;        % Total winding resistance (pu)
Trf.L1         = 0.12;         % Total Leakage inductance (pu)
Trf.Rm         = 1000; %pu
Trf.Lm         = 1000; %pu
Trf.ZbasePrim  = Trf.Vprim^2/Trf.Pnom;
Trf.Z1         = (Trf.R1 + Trf.L1*1i)*Trf.ZbasePrim; %ohm
Trf.Z0         = Trf.Z1; %ohm
Trf.ZbaseSec  = Trf.Vsecondary^2/Trf.Pnom;

%Transmission line parameters

TL.R1 = 0.01273;%TransmissionLine.R; %ohm/km
TL.R0 = 0.3864;%3*TransmissionLine.R; %ohm/km
TL.L1 = 0.9337e-3;%TransmissionLine.L; %H/km
TL.L0 = 4.1264e-3;%3*TransmissionLine.L; %H/km

TL.Length = 100; %km
TL.Freq = Fnom; %Hz
TL.C1 = 12.74e-9; %F/km
TL.C0 = 7.751e-9; %F/km
TL.Z1 = TL.Length*(TL.R1 + 1i*(TL.L1*(2*pi*TL.Freq))); %ohm
TL.Z0 = TL.Length*(TL.R0 + 1i*(TL.L0*(2*pi*TL.Freq))); %ohm 
TL.K0 = (TL.Z0-TL.Z1)/TL.Z1; %residual compensation factor

%Fault
Fault.Type = 'ABC-G'; %Options are {'No fault'  'A-G'  'AB'  'AB-G'  'ABC'  'ABC-G'}
Fault.m = 0.5; %pu length
Fault.R = 1; %ohm
Fault.StartTime = 10;%2.505; %Fault injection simulation time in s
Fault.EndTime = 12.6; %Fault removal simulation time in s

%set_param(FltProgBlock,'FaultType',Fault.Type); %Update the fault programming block in the model to use the selected fault type

%Converter
ConvPlt.Fc=Fnom*3.37;        % Carriers frequency (Hz)
ConvPlt.Nb_PM =36;                      % Number of power module per arm
ConvPlt.Vnom_dc = 640000;                % DC nominal voltage (V)
ConvPlt.C_PM= 1.758e-3; % Power module capacitor (F)
ConvPlt.Larm_pu = 0.15; %pu
ConvPlt.Rarm_pu = ConvPlt.Larm_pu/100; %pu
ConvPlt.Larm = ConvPlt.Larm_pu*(Trf.ZbaseSec/(2*pi*Fnom)); %Henry
ConvPlt.Rarm = ConvPlt.Rarm_pu*Trf.ZbaseSec; %ohm

% Energy in kJ/MVA
ConvPlt.W_kJ_MVA= 0.5 * ConvPlt.C_PM * (ConvPlt.Vnom_dc/ConvPlt.Nb_PM)^2 * ConvPlt.Nb_PM * 6 / (Trf.Pnom/1e6)/1e3;
ConvPlt.Vc0_PM=0;                     % Capacitors initial voltage (V)

%% Converter Control Parameters
% Sequencer timing:
HMI.Tbrk1_On=0;                 % Closing time of breaker 1 (converter energizing)
HMI.Tbrk2_On=0;                 % Closing time (s) of breaker 2 (across start-up resistor)
%
HMI.Tdeblock=0;                 % Converter de-block time (s)
HMI.Ton_VDCreg=0;               % VDC regulator turn-on time (s) - VDC Regulation
HMI.Tramping_Vdc_ref=2;           % Start time Vdc_ref ramping to nominal (s)
HMI.Slope_Vdc_ref=ConvPlt.Vnom_dc/5;      % Sloge Vdc_ref ramping (V/s)
%
HMI.Ton_PQreg=0.25;                  % Preg & Qreg regulators turn-on time (s) - PQ regulation
HMI.Tramping_Pref = HMI.Ton_PQreg;  % Start time Pref ramping(s)
HMI.Slope_Pref=0.5;               % Sloge Pref ramping (V/s)
HMI.Tramping_Qref = 9999;  % Start time Pref ramping(s)
HMI.Slope_Qref=0.5;               % Sloge Pref ramping (V/s)
%
HMI.Ton_Converter2=4;             % Converter 2 equivalent switched-on time (s)

%%
% *****************************************************************
%                         CONTROL PARAMETERS
% *****************************************************************

% Parameters
ConvCtrl.Freq = Fnom;
ConvCtrl.Ts=40e-6;   % Control system time step (s)

% dq and Vdc measurement filter cut-off frequency:
ConvCtrl.Fn_filter=1000;
ConvCtrl.Zeta_filter=1;

%PLL parameters
ConvCtrl.PLL.MinFreq = 0.9*Fnom;
ConvCtrl.PLL.InitFreq = Fnom;
ConvCtrl.PLL.Kp = 180;
ConvCtrl.PLL.Ki = 3200;
ConvCtrl.PLL.Kd = 1;
ConvCtrl.PLL.Ts = 1e-4;
ConvCtrl.PLL.MaxRocof = 12;
ConvCtrl.PLL.Fc_filter = 25;

% Active power regulator (Preg)
ConvCtrl.Preg.Kp= 0.5/3;                % Proportional gain
ConvCtrl.Preg.Ki= 1.0;                  % Integral gain
ConvCtrl.Preg.Limits = [ 1.2, 0.8 ] ;   % Output (Vdc_ref) Upper/Lower limits (pu)

% Reactive power regulator (Qreg)
ConvCtrl.Qreg.Kp= 0.5/3;                % Proportional gain
ConvCtrl.Qreg.Ki= 1.0;                  % Integral gain
ConvCtrl.Qreg.Limits = [ 1.2, -1.2]; % Output (Iq_ref) Upper/Lower limit (pu)

% VDC regulator (VDCreg)
ConvCtrl.VDCreg.Kp =4;                   % Proportional gain
ConvCtrl.VDCreg.Ki=100;                 % Integral gain
ConvCtrl.VDCreg.Limits= [ 1.2  -1.2];   % Output Idref [Upper Lower] limits (pu)

% Current regulator (Ireg)
ConvCtrl.Ireg.Kp = 0.6/2;                  % Proportional gain
ConvCtrl.Ireg.Ki = 6/2;                    % Integral gain
ConvCtrl.Ireg.Limits = [ 1.2  -1.2];     % Output Vdq_conv [Upper Lower] limits (pu)

% Feedforward coefficients:
ConvCtrl.FF.L = 0.0750;
ConvCtrl.FF.R = 7.5000e-04;

%% Clean up
clear Model FltProgBlock Fnom