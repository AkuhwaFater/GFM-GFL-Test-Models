function PlotDistanceElement(TL,Data,FigNumber,shape,size,resetColor)
        
    if ~ishandle(FigNumber)
    figure(FigNumber);
    set(gcf,'Position',[10 10 900 600]);
    set(gcf,'Renderer','painters');
    
    plot([0+0*1i TL.Z1],'b','LineWidth',2,'DisplayName','Line Impedance')
    grid on
    hold on
    %Settings
    Reactive_Reach = 0+imag(TL.Z1*0.8)*1i; %Reactive Reach
    Resistive_Reach = 50; %Resistive reach
    alpha_1_Left_Blinder = 115; %Left Blinder angle
    alpha_2_Lower_Blinder = -2; %Lower Blinder angle
    
    Lower_Blinder = Resistive_Reach/cos(deg2rad(alpha_2_Lower_Blinder))*exp(1i*deg2rad(alpha_2_Lower_Blinder));
    Left_Blinder = abs(Reactive_Reach)/sin(deg2rad(alpha_1_Left_Blinder))*exp(1i*deg2rad(alpha_1_Left_Blinder));
    plot([0 Lower_Blinder],'r','LineWidth',2,'HandleVisibility','off'); %Lower Blinder
    plot([Resistive_Reach+imag(Lower_Blinder)*1i Resistive_Reach+Reactive_Reach],'r','LineWidth',2,'HandleVisibility','off'); %right Blinder
    plot([real(Left_Blinder)+Reactive_Reach Reactive_Reach+Resistive_Reach],'r','LineWidth',2,'HandleVisibility','off'); %Top Blinder
    plot([0 Left_Blinder],'r','LineWidth',2,'HandleVisibility','off'); %left Blinder
    
    grid on
    else 
    figure(FigNumber);
    end

    if resetColor == 1
        ax = gca;
        ax.ColorOrderIndex = 1;
    end
    
    plot(Data.mZ1L_Est.Data(end),shape,'LineWidth',2,'MarkerSize',size);
    legend('Interpreter','Latex');
    xlabel('Resistance (\Omega)','FontSize',13);
    ylabel('Reactance (\Omega)','FontSize',13);
    
end