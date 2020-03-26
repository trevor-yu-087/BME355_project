function [f, lm_n] = MuscleForce(a, Lm, Vm, Lceopt, Fmax, type)


    Br = 0.28;
    Ar = 0.17;
%     Lceopt = 0.116;
%     Lceopt = 0.128647692672029;
    
    Kpe = 5;
    E0m = 0.6;
%     Fmax = 624.3;
    % Lm = linspace(0,2*Lceopt, 41);
    % Vm = linspace(-4*0.055/0.1,4*0.055/0.1,41);
    width = 0.888;
    c = -1/width^2;
    gamma = 0.45;
    tau = 0.1;
    f_assymp = 1.6;
    Kce1 = 0.25;
    Kce2 = (f_assymp-1)/2*Kce1/(1+Kce1);

    
   
    for i=1:length(Lm)
        lm_n(i) = Lm(i) / Lceopt;
        vm_n(i) = Vm(i)/(Lceopt/tau);
        f_iso(i) = c * lm_n(i)^2 - 2 * c * lm_n(i) + c +1;
    end

    switch type
        case 'Soest'
            for i=1:length(Lm)
                % Force-Lenght relationship
                FL_Ak(i) = (a(i)*Br*(f_iso(i)+Ar)-a(i)*Ar*(Br))/(Br);
                
                % Force-Velocityb2=-f_iso*f_assymp;    b2V=-f_assymp;
                b2V = -f_assymp;
                b1V = (Br*(1 + b2V)^2)/(a(i) + Ar);
                b3V = b1V/(1 + b2V);
                
                if (vm_n(i) <=0 )
                    FV_Ak(i)=(a(i)*Br*(1+Ar)-a(i)*Ar*(Br-vm_n(i)))/(Br-vm_n(i));
                else
                    FV_Ak(i)=(b1V*a(i)-b2V*a(i)*(b3V-vm_n(i)))/(b3V-vm_n(i));
                end
                
                f (i) = Fmax * FV_Ak(i) * a * FL_Ak(i);
            end
            
        case 'Silva'
            for i=1:length(Lm)
                % Force-Length
                FL_sk (i) = exp(-((9/4*(lm_n(i)-19/20))^4+1/4*(-9/4*(lm_n(i)-19/20))^2));
                
                % Force-Velocity
                if (vm_n(i)<-1)
                    FV_sk(i) = 0;
                elseif (vm_n(i) >= -1 && vm_n(i)<=0.2)
                    FV_sk(i) = 1 -1 / atan(5) * atan(-5 * vm_n(i));
                else
                    FV_sk(i) = pi/(4*atan(5))+ 1;
                end
                
                f(i) = Fmax * a(i) * FL_sk (i) * FV_sk(i);
                
            end
            
        case 'Thelen'
            for i=1:length(Lm)
                % Force-Length
                
                FL_ou(i) = exp(-(lm_n(i)-1)^2/gamma);
                
                % Force-Velocity
                
                if (vm_n(i) <= -1)
                    FV_ou(i) = 0;
                elseif (vm_n(i) >= -1 && vm_n(i) <= 0)
                    FV_ou(i) = (1 + vm_n(i))/(1 - vm_n(i)/Kce1);
                else
                    FV_ou(i) = (1 + vm_n(i) * f_assymp/Kce2)/(1 + vm_n(i)/Kce2);
                end
                f(i) = Fmax * a(i) * FL_ou (i) * FV_ou(i);
            end
            
    end





