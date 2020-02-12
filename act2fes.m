function [stim, e] = act2fes(t, act, modulation, Usat, Utr, FT, gamma, lF)

% Inputs:

    if nargin<8; lF = 0.144; end           % Fiber length [m]
    if nargin<7; gamma = 5; end            % Pennation angle [º]
    if nargin<6; FT = 1250; end            % Muscle force [N]
    if nargin<5; Utr = 20; end             % Treshold amplitude
    if nargin<4; Usat = 35; end            % Saturation value
    if nargin<3; modulation = 1; end       % Type of modulation: 1 Amplitude; 2 frequency
    
    S0 = 45;                               % Muscle tension [N/cm^2]  (see Gfohler 2002)
    dm = 1054;                             % Muscle density [kg/m^3]
    Cfast = 0.35;                          % Percentage of fast fibres    
    time = t;
    
    N = length(time);

    PCSA = FT /(S0*cosd(gamma));
    m = PCSA * lF * dm /10000;

    % Time constants (see Gfohler 2004)
    Te = (25+0.1*m*(1-Cfast))/1000;
    Trise = (5+0.2*m*(1-Cfast)^2)/1000*15;
    Tfall = (30+0.5*m*(1-Cfast)^2)/1000*1.9;
    
    a = act;
    a_p = gradient(act)./mean(diff(time));
    a_pp = gradient(a_p)./mean(diff(time));
    e = zeros(N,1);
    
    
    for i = 1:N
        if time<time(end)
            if a(i)<a(i+1)
                k1 = Te*Trise;
                k2 = Trise+Te;
            else
                k1 = Te*Tfall;
                k2 = Tfall+Te;
            end
        else
            k1 = Te*Tfall;
            k2 = Tfall+Te;
        end
        if k1*a_pp(i) + k2*a_p(i) + a(i) >= 0
            e(i) = k1*a_pp(i) + k2*a_p(i) + a(i);
        else
            e(i)=0;
        end
    end
    
   
    R = 15;         % See Gofhler et al. (2004) for values
    a2 = 2.5;       % See Gofhler et al. (2004) for values

    % Watanabe et al. (1999) IEEE Trans Rehab Eng 7 (1)
    rCF = (R-91.2)/(-1.03);             % -->  R = 91.2-1.03*rCF
    r0 = R*log((a2-1)*exp(rCF/R)-a2);
    a1 = -a2*exp(-r0/R);
    
    stim = zeros(N-2,1);
    
    switch modulation
        case 1 % Amplitude
            fstim = 50;
            Sv = (-a2*(1-exp(-r0/R)))/(1+exp(fstim-r0/R))+a2;
            
%             Sv = ((a1 - a2)/(1+exp((fstim-r0)/R))+a2);%*0.2772;
%             Utr = (5+((70-1)*(50-5))/((127-1)));
%             Usat = (5+((110-1)*(50-5))/((127-1)));
%             Utr = 20;
%             Usat = 35;
%             Su = e./Sv;
%             Su = e./(max(e)-min(e));
            stim = e.*(Usat-Utr)+Utr;
                      

        case 2 % Frequency
            for i = 1:length(e)
                Ustim = (5+((109-1)*(50-5))/((127-1)));
                Utr = (5+((70-1)*(50-5))/((127-1)));
                Usat = (5+((110-1)*(50-5))/((127-1)));
                Su = (Ustim-Utr)/(Usat-Utr);
                stim(i) = R*log((a1-a2)/(e(i)/Su-a2)-1)+r0;
            end
            
    end

%     plot(time, stim, 'linewidth', 2); hold on; 
%     plot(time, act, 'r','linewidth', 2);
%     legend('Inverted Signal', 'Activation Signal'); grid on;
    disp('>> end act2fes')
end

