% Select the operating scenarios
mode = [1]; %1 for GFM 0 for GFL
% Test Cases Parameters
faultImpedanceArray = [0.001 1 5 10 15];%[0.05 0.1 0.25 0.4]; % Three phase fault impedance
faultDistanceArray = [0.2 0.5 0.7];
gridSCLArray = [5 20 40]*1e9; %GVA

%Run test cases
% Defining simulink simulation input
simIn = Simulink.SimulationInput('HVDCMMCInterconnection');
%set_param('HVDCMMCInterconnection', 'SignalLogging', 'off');
for GFMGFLSelect = mode
    disp(['Control Mode = ' num2str(GFMGFLSelect)])
    resIdx = 1;
    try
        clear TestCaseResults
    end

    for k = 1: length(gridSCLArray)
        tic
        Psc = gridSCLArray(k); % pu
        for j = 1:length(faultDistanceArray)
            Fault.m = faultDistanceArray(j);
        
            for i = 1:length(faultImpedanceArray)
                
                Fault.R = faultImpedanceArray(i); % pu
            
                % Running the simulation
                disp(['Running Test Case: Rf = ' num2str(Fault.R) ' ohm ; m = ' num2str(Fault.m) ' pu ; SCL = ' num2str(Psc/1e9) ' GVA']);
                warning('off','all');
                outData = sim(simIn);
    
                TestCaseResults(resIdx).mZ1L_Est = outData.mZ1L_Est;
                TestCaseResults(resIdx).V1_Local = outData.V1_Local;
                TestCaseResults(resIdx).I1_Local = outData.I1_Local;
                TestCaseResults(resIdx).V1_Remote = outData.V1_Remote;
                TestCaseResults(resIdx).I1_Remote = outData.I1_Remote;
                TestCaseResults(resIdx).faultResistance = Fault.R;
                TestCaseResults(resIdx).faultDistance = Fault.m;
                TestCaseResults(resIdx).SCR = Psc;
                resIdx = resIdx +1;
                
                SimTime = toc;
            end
        end
    end
    set_param('HVDCMMCInterconnection', 'SignalLogging', 'on');
    
    %Save results to .mat file
    FilePath = [pwd '\SimData\'];
    
    if GFMGFLSelect == 1
        save([FilePath 'HVDCMMC_TestCasesResults_GFM.mat'],"TestCaseResults");
    else
        save([FilePath 'HVDCMMC_TestCasesResults_GFL.mat'],"TestCaseResults");
    end
end

toc

%% Plot Results for GFL
TestCasesResults_GFL = load('HVDCMMC_TestCasesResults_GFL.mat');
TestCasesResults_GFL = TestCasesResults_GFL.TestCaseResults;
shape = 'x';

%plot figure for varying fault resistance
%Results_VaryingFaultResistance = TestCasesResults_GFL(1:4);
Rf = [];
m = 0.5;
SCL = 20e9;
Results_VaryingFaultResistance = FindTestCases(TestCasesResults_GFL, Rf, m, SCL);

FigNumber = 100;
for i = 1:length(Results_VaryingFaultResistance)
    PlotDistanceElement(TL,Results_VaryingFaultResistance(i),FigNumber,shape)
end
figure(FigNumber)
title('Distance Element With Varying Fault Resistance - GFL')

%plot figure for varying fault distance
%Results_VaryingFaultDistance = TestCasesResults_GFL(1:4:12);

Rf = 0.001;
m = [];
SCL = 20e9;
Results_VaryingFaultDistance = FindTestCases(TestCasesResults_GFL, Rf, m, SCL);

FigNumber = 101;
for i = 1:length(Results_VaryingFaultDistance)
    PlotDistanceElement(TL,Results_VaryingFaultDistance(i),FigNumber,shape)
end
figure(FigNumber)
title('Distance Element With Varying Fault Distance - GFL')

%plot figure for varying SCR
%Results_VaryingSCR = TestCasesResults_GFL(1:12:36);

Rf = 0.001;
m = 0.5;
SCL = [];
Results_VaryingSCR = FindTestCases(TestCasesResults_GFL, Rf, m, SCL);

FigNumber = 102;
for i = 1:length(Results_VaryingSCR)
    PlotDistanceElement(TL,Results_VaryingSCR(i),FigNumber,shape)
end
figure(FigNumber)
title('Distance Element With Varying Grid SRC - GFL')

%% Plot Results for GFM
TestCasesResults_GFM = load('HVDCMMC_TestCasesResults_GFM.mat');
TestCasesResults_GFM = TestCasesResults_GFM.TestCaseResults;
shape = 'o';

%plot figure for varying fault resistance
%Results_VaryingFaultResistance = TestCasesResults_GFM(1:4);
Rf = [];
m = 0.5;
SCL = 20e9;
Results_VaryingFaultResistance = FindTestCases(TestCasesResults_GFM, Rf, m, SCL);

FigNumber = 200;
for i = 1:length(Results_VaryingFaultResistance)
    PlotDistanceElement(TL,Results_VaryingFaultResistance(i),FigNumber,shape)
end
figure(FigNumber)
title('Distance Element With Varying Fault Resistance - GFM')

%plot figure for varying fault distance
%Results_VaryingFaultDistance = TestCasesResults_GFM(1:4:12);

Rf = 0.001;
m = [];
SCL = 20e9;
Results_VaryingFaultDistance = FindTestCases(TestCasesResults_GFM, Rf, m, SCL);

FigNumber = 201;
for i = 1:length(Results_VaryingFaultDistance)
    PlotDistanceElement(TL,Results_VaryingFaultDistance(i),FigNumber,shape)
end
figure(FigNumber)
title('Distance Element With Varying Fault Distance - GFM')

%plot figure for varying SCR
%Results_VaryingSCR = TestCasesResults_GFM(1:12:36);

Rf = 0.001;
m = 0.5;
SCL = [];
Results_VaryingSCR = FindTestCases(TestCasesResults_GFM, Rf, m, SCL);

FigNumber = 202;
for i = 1:length(Results_VaryingSCR)
    PlotDistanceElement(TL,Results_VaryingSCR(i),FigNumber,shape)
end
figure(FigNumber)
title('Distance Element With Varying Grid SRC - GFM')

%% Function
function Results = FindTestCases(TestCases, Rf, m, SCR)

    if isempty(Rf)
        Results = TestCases(cell2mat({TestCases.faultDistance}) == m);
        Results = Results(cell2mat({Results.SCR}) == SCR);
    elseif isempty(m)
        Results = TestCases(cell2mat({TestCases.faultResistance}) == Rf);
        Results = Results(cell2mat({Results.SCR}) == SCR);
    elseif isempty(SCR)
        Results = TestCases(cell2mat({TestCases.faultDistance}) == m);
        Results = Results(cell2mat({Results.faultResistance}) == Rf);
    else
        Results = TestCases(cell2mat({TestCases.faultDistance}) == m);
        Results = Results(cell2mat({Results.faultResistance}) == Rf);
        Results = Results(cell2mat({Results.SCR}) == SCR);
    end

end