function TestCaseResults = FaultAnalysisSimRun(faultImpedanceArray, faultDistanceArray, gridSCRArray, testCondition, Fault)
    % Copyright 2023 The MathWorks, Inc.

    testConditionSelected = testCondition; % Storing the selected test condition
    
    % Setting up the grid forming converter parameters
    run("GridFormingConverterInputParameters.mlx");
    
    % Applying the original test condition
    testCondition.activePowerMethod = testConditionSelected.activePowerMethod;
    testCondition.currentLimitMethod = testConditionSelected.currentLimitMethod;
    testCondition.testCondition = testConditionSelected.testCondition;
    testCondition.XbyR = testConditionSelected.XbyR;
    testCondition.SCR = testConditionSelected.SCR;
    
    % Defining function output variable
    run('GridFormingConverterTestCondition.mlx');
    
    
    % Defining simulink simulation input
    simIn = Simulink.SimulationInput('GridFormingConverter');
    resIdx = 1;
    for k = 1: length(gridSCRArray)
        testCondition.SCR = gridSCRArray(k); % pu
        for j = 1:length(faultDistanceArray)
            Fault.m = faultDistanceArray(j);
    
            for i = 1:length(faultImpedanceArray)
            
                testCondition.faultResistance = faultImpedanceArray(i); % pu
        
                % Setting up the simulation test condition
                simIn = setVariable(simIn,'testCondition',testCondition);
                simIn = setVariable(simIn, 'Fault',Fault);
                % Running the simulation
                outData = sim(simIn);
            
                TestCaseResults(resIdx).mZ1L_Est = outData.mZ1L_Est;
                TestCaseResults(resIdx).faultResistance = testCondition.faultResistance;
                TestCaseResults(resIdx).faultDistance = Fault.m;
                TestCaseResults(resIdx).SCR = testCondition.SCR;
                resIdx = resIdx +1;
            end
        end
    end
end

