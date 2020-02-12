function [a,e] = fes2act(t, stimFrec, stimAmp, plotFlag)

% -------------------------------------------------------------------------
% Author: Francisco Romero Sánchez
% Dept. of Mechanical, Energetics and Meterials Engineering
% University of Extremadura
% 
% fes2act generates a '.mat' file containing muscle activation signal from
% a FES profile modulated either in amplitude or frequency.
%
% Input: 
%   - stimFrec: array containing fes signal in a FREQUENCY modulation
%   profile. If stimAmp is selected as input, stimFrec must be an array of
%   ones of the same length as stimAmp.
%   - stimAmp: array containing fes signal in a AMPLITUDE modulation
%   profile. If stimFrec is selected as input, stimAmp must be an array of
%   ones of the same length as stimFrec.
%               
% Output: 
% 'Activ.mat': file containing 
%   - t: time vector
%   - e: ecitation signal prior to linear dynamics in Hammerstein model
%   - a: resulting activation profile
%   - ampFES & frecFES: input signals
% Function fes2act returns a and e values to workspace
%--------------------------------------------------------------------------

if nargin<4; plotFlag = 1; end
if nargin~=3
    disp('Time, amplitude and frequency profiles are required'); 
    return; 
end


    N = size(stimFrec,1);
    ampFES = stimAmp;
    frecFES = stimFrec;

    Utres = (5+((70-1)*(50-5))/((127-1)));  % Measured on test
    Usat = (5+((110-1)*(50-5))/((127-1)));  % Measured on test
    R = 40;                                 % See Gfohler et al. (2004)
    a2 = 2.5;                               % See Gfohler et al. (2004)
    Su = zeros(N,1);    Sv = zeros(N,1);    e  = zeros(N,1);
    
    for i = 1:N
        if ampFES(i) < Utres
            Su(i) = 0; end
        if ampFES(i)>= Utres && ampFES(i) < Usat
            Su(i) = (ampFES(i)-Utres)/(Usat-Utres); end
        if ampFES(i)>Usat
            Su(i) = 1;
        end
        Su(i) = stimAmp(i);
        % Watanabe et al. (1999) IEEE Trans Rehab Eng 7 (1) 12-18
        rCF = (R-91.2)/(-1.03);          % -->  R = 91.2-1.03*rCF
        r0 = R*log((a2-1)*exp(rCF/R)-a2);
        a1 = -a2*exp(-r0/R);
        Sv(i) = (a1 - a2)/(1 + exp((frecFES(i)-r0)/R)) + a2;
        e(i) = Su(i)*Sv(i);
    end
    save FESdata.mat t e

    x0 = [0.00001, 0];
    tspan=t;

    % Gfolher et al. J. of Mech Med Biol, 2004, 4, 77-92
    FT = 2549;              % [N]
    S0 = 450000;            % [N/cm^2]
    gamma = 10 * pi/180;    % [rad]
    dm = 1054;              % [kg/m^3]
    lF = 0.1;               % [m]
    Cfast = 0.52;

    PCSA = FT /(S0*cos(gamma));
    m = PCSA * lF * dm;

    global Te Trise Tfall
    Te = (25+0.1*m*(1-Cfast))/1000;
    Trise = (5+0.2*m*(1-Cfast)^2)/1000*15;
    Tfall = (30+0.5*m*(1-Cfast)^2)/1000*19;

    [time, a] = ode23(@FES, tspan, x0);
    [time2, Nact] = ode23(@FES2, tspan, x0);

    figure; 
    plot(time, e, 'linewidth', 2); hold on; plot(time, a(:,1),'r','linewidth',2);
    legend('e(t)', 'a_{FES}(t)');
    if plotFlag == 1
        hold on; plot(time2,Nact,'g','linewidth',2);
        legend('e(t)', 'a_{FES}(t)', 'a_{PHY}(t)');
    end
     xlabel('Time [s]'); ylabel('Signal'); grid on;

    save Activ.mat t e a ampFES frecFES
    disp('>> end fes2act')
end

% FES2 --> calculate activations for an ARTIFICIALLY activated muscle
function dx = FES(t,x)

    global Te Trise Tfall       
    A = load('FESdata.mat');
    time = A.t;
    e = A.e;
    e_int= interp1(time,e,t);
    e_int2=interp1(time,e,t+0.01);
    
    if t<(time(end))
        if e_int < e_int2
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
    
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = 1/k1*(-k2*x(2)-x(1)+e_int);
end

% FES2 --> calculate activations for a PHYSIOLOGICALLY activated muscle
function act = FES2(t,x)

global Trise Tfall       
    A = load('FESdata.mat');
    time = A.t;
    e = A.e;
    e_int= interp1(time,e,t);
    act=(e_int-x)*((e_int/Trise)+(1-e_int)/Tfall);
end


