% Select the operating scenarios
testCondition.testCondition = 'Permanent three-phase fault'; 
% Test Cases Parameters
faultImpedanceArray = [0.005 1 5 50];%[0.05 0.1 0.25 0.4]; % Three phase fault impedance
faultDistanceArray = [0.2 0.5 0.7];
gridSCRArray = [1 2.5 4];
%Run test cases
TestCaseResults = FaultAnalysisSimRun(faultImpedanceArray, faultDistanceArray, gridSCRArray, testCondition, Fault);

%% Plot Results
load('TestCasesResults_GFL.mat')

%plot figure for varying fault resistance
%Results_VaryingFaultResistance = TestCaseResults(1:4);
Rf = [];
m = 0.5;
SCR = 2.5;
Results_VaryingFaultResistance = FindTestCases(TestCaseResults, Rf, m, SCR);

FigNumber = 100;
for i = 1:length(Results_VaryingFaultResistance)
    PlotDistanceElement(TL,Results_VaryingFaultResistance(i),FigNumber)
end
figure(FigNumber)
title('Distance Element With Varying Fault Resistance')

%plot figure for varying fault distance
%Results_VaryingFaultDistance = TestCaseResults(1:4:12);

Rf = 0.005;
m = [];
SCR = 2.5;
Results_VaryingFaultDistance = FindTestCases(TestCaseResults, Rf, m, SCR);

FigNumber = 101;
for i = 1:length(Results_VaryingFaultDistance)
    PlotDistanceElement(TL,Results_VaryingFaultDistance(i),FigNumber)
end
figure(FigNumber)
title('Distance Element With Varying Fault Distance')

%plot figure for varying SCR
%Results_VaryingSCR = TestCaseResults(1:12:36);

Rf = 0.005;
m = 0.5;
SCR = [];
Results_VaryingSCR = FindTestCases(TestCaseResults, Rf, m, SCR);

FigNumber = 102;
for i = 1:length(Results_VaryingSCR)
    PlotDistanceElement(TL,Results_VaryingSCR(i),FigNumber)
end
figure(FigNumber)
title('Distance Element With Varying Grid SRC')


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