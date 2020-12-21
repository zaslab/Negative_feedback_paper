function [t , Ca, L, I, R_a] = Simulation(varargin)
    %% parameters    
    
    %General consts   
        %light mode for debug
        ez                                  =   true;
        
        %1 = only our model. 2 = including previous H.H model for AWA
        model_type                          =   1;
        
        %ending time of simulation in mili seconds
        t_end                               =   100000;
        
        %print time tracking
        time_track                          =   true;
    
        %decay time of calcium concentration in mili seconds 
        t_c                                 =   4e3;
        
        %initial value of inhibition
        I_initial_val                       =   0;
 
        %steady state calcium concentration in mollar
        Ca_0                                =   100e-9;
        
    %Receptor activation function.
        K2                                  =   10;   
        K1                                  =   25;
        L_0                                 =   1; % In uM
        
    %Self activation function
        R_t                                 =   0.95;
        S_min                               =   1e-6;
        S_max                               =   1;
        K4                                  =   1e-07;  

    %InputLigandProfile
    %gradType: 1 = linear
    %           2 = step
    %           3 = sigmoid
    %           4 = stairs
    %           5 = multi_steps
    %           6 = exponential
        gradType                            =   2;    
        baseline                            =   0.1150; %in uM   
        ts                                  =   10000; 
        tf                                  =   100000; 
        lin_slope                           =   1150/1000000;
        amplitude                           =   1150;
        step_duration                       =   30000;
        break_duration                      =   60000;
    %Modification function. 
        K5                                  =   5;
        K6                                  =   2e-6;
        t_I                                 =   3e5; %in ms 
        
    % only relevant if intergrating with previous model for AWA H.H
        c_current                           =   35;
    %% parse
    p = inputParser;
        
    addOptional(p,'ez',ez);
    addOptional(p,'model_type',model_type);
    addOptional(p,'t_end',t_end);
    addOptional(p,'time_track',time_track);
    addOptional(p,'calcium_decay_time_const',t_c);
    addOptional(p,'I_initial_val',I_initial_val);
    addOptional(p,'Ca_0',Ca_0);

    addOptional(p,'C_I',K2);
    addOptional(p,'C_L',K1);
    addOptional(p,'L_0',L_0);
    
    addOptional(p,'critic_receptor_activation',R_t);
    addOptional(p,'S_min',S_min);
    addOptional(p,'S_max',S_max);
    addOptional(p,'S2Ca_coeff',K4);
   
    addOptional(p,'gradType',gradType);
    addOptional(p,'baseline',baseline);
    addOptional(p,'ts',ts);
    addOptional(p,'tf',tf);
    addOptional(p,'lin_slope',lin_slope);
    addOptional(p,'amplitude',amplitude);
    addOptional(p,'step_duration',step_duration);
    addOptional(p,'break_duration',break_duration);
    
    addOptional(p,'K5',K5);
    addOptional(p,'K6',K6);
    addOptional(p,'inhibition_decay_time_const',t_I);

    addOptional(p,'c_current',c_current);
    
    parse(p,varargin{:});
    ourparams       = zeros(1,100);
    
    ez              = p.Results.ez;
    model_type      = p.Results.model_type;
    t_end           = p.Results.t_end;
    time_track      = p.Results.time_track;
    ourparams(3)    = p.Results.calcium_decay_time_const; 
    ourparams(4)    = p.Results.I_initial_val;   
    ourparams(5)    = p.Results.Ca_0;
   
    ourparams(21)   = p.Results.C_I;
    ourparams(22)   = p.Results.C_L;
    ourparams(25)   = p.Results.L_0;
    
    ourparams(31)   = p.Results.gradType;
    ourparams(32)   = p.Results.baseline;
    ourparams(33)   = p.Results.ts;
    ourparams(34)   = p.Results.tf;  
    ourparams(35)   = p.Results.lin_slope; 
    ourparams(36)   = p.Results.amplitude; 
    ourparams(37)   = p.Results.step_duration; 
    ourparams(38)   = p.Results.break_duration; 
    
    ourparams(41)   = p.Results.K5; 
    ourparams(42)   = p.Results.K6;
    ourparams(43)   = p.Results.inhibition_decay_time_const; 

    ourparams(51)   = p.Results.critic_receptor_activation;
    ourparams(52)   = p.Results.S_min;
    ourparams(53)   = p.Results.S_max;
    ourparams(54)   = p.Results.S2Ca_coeff;
    
    ourparams(60)   = p.Results.c_current;
    %% general parameters of the simulation
    %frequency of data saving
    savingResolution            =   50;
    if ez, savingResolution     =   100; end

    %minimum allowed dt in mili seconds
    min_dt                      =   0.0001;
    if ez, min_dt               =   0.001; end

    %maximum allowed dt in mili seconds
    max_dt                      =   0.5;

    %upper bound for the vectors size
    vec_size                    =   floor(t_end/min_dt/savingResolution);

    %initialize the times vector
    t                           =   nan(vec_size,1);

    %initialize the self-amplification motif vector
    S                           =   nan(vec_size,1);
    
    %initialize Ligand gradient
    L                           =   nan(vec_size,1);

    %initialize receptor activation 
    R_a                         =   nan(vec_size,1);

    %initialize Calcium concentration 
    Ca                          =   nan(vec_size,1);

    %initialize inhibition level 
    I                           =   nan(vec_size,1);
    
    %initializing parameters which required for intergrated model
    if (model_type == 2)
        %number of physical parametes in model
        gating_param_num        =   7;   
        %initialize matrix which holds the gating parameters along time
        gating_params           =   nan(vec_size,gating_param_num);
        %initialize the Voltage vector
        V                       =   nan(vec_size,1);
        %initialize the input current vector
        E_C                       =   nan(vec_size,1);
    end
    
    %% setting starting parameters
    %initializing time
    t(1)                        =   0;

    %initializing odorant concentration
    L(1)                        =   InputLigandProfile(t(1), p.Results.gradType, p.Results.baseline, p.Results.ts,...
                                    p.Results.tf, p.Results.lin_slope, p.Results.amplitude, p.Results.step_duration);

    %initializing Calcium concentration             
    Ca(1)                       =   ourparams(5);

    %initializing inhibition level
    I(1)                        =    ourparams(4);
    
    %initializing receptor activation
    R_a(1)                      =   Receptor_activation(L(1) ,I(1),ourparams(21),ourparams(22), ourparams(25));

    %initializing self-amplification value
    S(1)                        =   ourparams(52);
    
     %starting parameters which required for intergrated model
    if (model_type == 2)
        %initial Votage in mV
        V(1)                    =   -75.03;      

        %initial Calicum 1 gating parameter (C_1)
        gating_params(1,1)      =   1.13e-5;                   

        %initial Calicum 1 gating parameter (C_2)
        gating_params(1,2)      =   1.13e-5;                      

        %initial pottassium SHK1 gating parameter (C_2)
        gating_params(1,3)      =   2.04e-7;                      

        %initial pottassium Kr gating parameter (Kr)
        gating_params(1,4)      =   0;                            

        %initial High V 1 pottasoim channels gating parameter
        gating_params(1,5)      =   3.48e-5;                      

        %initial pottassium Kr gating parameter (Kra)
        gating_params(1,6)      =   0;                            

        %initial High V 2 pottasoim channels gating parameter
        gating_params(1,7)      =   1.24e-9;   
        
        %initial current
        E_C(1)                    =    0;
    end
    %% running the numeric simulation
    curr_t                      = t(1);
    curr_Ca                     = Ca(1);
    curr_L                      = L(1);
    curr_I                      = I(1);
    curr_R_a                    = R_a(1);
    curr_S                      = S(1);
    if (model_type == 2)
        curr_V = V(1);
        curr_E_C = E_C(1);
        curr_gating_params = gating_params(1,:);
    end
    
    
    %counter
    i = 1;
    while (curr_t < t_end)
        
        if (model_type == 1) 
            [d_Ca_dt, d_I_dt, d_S_dt] = fullModel(curr_Ca,curr_R_a,curr_S,curr_I,ourparams);
        
            %get the needed resolution of dt for this step
            min_needed_dt       =   1/max(abs([ d_I_dt,d_S_dt]));
            dt                  =   max(min(min_needed_dt,max_dt),min_dt);
        elseif (model_type == 2)
            [d_gating_params_dt, d_V_dt, d_Ca_dt, d_I_dt, d_S_dt] = fullIntegratedModel(curr_gating_params,...
                        curr_V, curr_E_C, curr_Ca, curr_R_a, curr_S, curr_I, ourparams);
            min_needed_dt       =   1/max(abs([d_gating_params_dt, d_V_dt, d_I_dt, d_S_dt]));
            dt                  =   max(min(min_needed_dt,max_dt),min_dt);
        end
            
        %get the changes in params for this step
        d_Ca                        = dt * d_Ca_dt;
        d_I                         = dt * d_I_dt;
        d_S                         = dt * d_S_dt;
        
        if (model_type == 2)
            d_gating_params         = dt * d_gating_params_dt;
            d_V                     = dt * d_V_dt;
        end
        
        %update curr data for this step
        curr_t                      =   curr_t + dt;
        curr_Ca                     =   curr_Ca + d_Ca;
        curr_L                      =   InputLigandProfile(curr_t,ourparams(31),ourparams(32),ourparams(33),...
                                        ourparams(34),ourparams(35),ourparams(36),ourparams(37),ourparams(38));
        curr_I                      =   curr_I + d_I; 
        curr_R_a                    =   Receptor_activation(curr_L,curr_I,ourparams(21),ourparams(22), ourparams(25));
        curr_S                      =   curr_S + d_S;
        curr_S                      =   min(max(curr_S,ourparams(52)),ourparams(53)); 
        if (model_type == 2)
            curr_V                  =   curr_V + d_V;
            curr_gating_params      =   curr_gating_params + d_gating_params;
            curr_E_C                =   ourparams(60) * curr_S;
        end
        i                           =   i + 1;
       
        
        
        %saving the data for this step
        if (mod(i,savingResolution) == 0)
            ii                      =   i/savingResolution;
            t(ii+1)                 =   curr_t;
            Ca(ii+1)                =   curr_Ca;
            L(ii+1)                 =   curr_L;
            I(ii+1)                 =   curr_I;
            R_a(ii+1)               =   curr_R_a;
            S(ii+1)                 =   curr_S; 
            if (model_type == 2)
                V(ii+1)             =   curr_V;
                %gating_params(ii+1,:) =   curr_gating_params;
                E_C(ii+1)           =   curr_E_C;
            end
        end

        %follow the process
        if (mod(i,10^5) == 0 && time_track)
            disp(num2str(curr_t)); 
        end
    end

    %erasing unused space in all vectors and matrices
    t                               =   t(~isnan(t)); 
    Ca                              =   Ca(~isnan(Ca));
    L                               =   L(~isnan(L));
    I                               =   I(~isnan(I));
    R_a                             =   R_a(~isnan(R_a));
    S                               =   S(~isnan(S));
    if (model_type == 2)
        V                           =   V(~isnan(V));
        %gating_params              =   gating_params(1:length(t),:);
        E_C                         =   E_C(~isnan(E_C));
    end
    
    %% drawing the data from simulation
    subplot(12,1,1:2);
    plot(t/1000,L*10^-3,'LineWidth', 1.5, 'Color', 'r'); 
    ylabel({'Ligand' , '[mM]'});
    ylim([-0.2*max(L*10^-3),1.2*max(L*10^-3)]);
    xlim([min(t/1000) max(t/1000)]);
    set(gca,'FontSize',16);
    set(gca,'XTickLabel',[]);
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'of';

    subplot(12,1,3:7); 
    plot(t/1000,Ca*10^6,'LineWidth', 1.5, 'Color', 'b'); 
    ylabel({'Calcium conc.','[\muM]'});
    ylim([-0.1*max(Ca*10^6),1.1*max(Ca*10^6)]);
    xlim([min(t/1000) max(t/1000)]);
    set(gca,'FontSize',16);
    set(gca,'XTickLabel',[]);
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'of';

    subplot(12,1,8:12); 
    yyaxis right; plot(t/1000,I,'LineWidth', 1.5); 
    ylabel('Inhibition');
    ylim([-0.03*max(I),1.03*max(I)]);
    yyaxis left; plot(t/1000,R_a,'LineWidth', 1.5);
    ylabel('Receptor Activation');
    ylim([-0.05,1.05]);
    xlabel('Time [Sec]');
    xlim([min(t/1000) max(t/1000)]);
    set(gca,'FontSize',16);
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'of';
    %% main model function
    function [d_Ca_dt, d_I_dt, d_S_dt] = fullModel(curr_Ca, curr_R_a, curr_S,curr_I, ourparams)
    
        % tracking the calcium concentration derivative in M/ms. 
        d_Ca_dt     = ourparams(54)*curr_S - (1/ourparams(3))*(curr_Ca - ourparams(5));

        %change in inhibition level based on calcium level in cell    
        d_I_dt      = Inhibition_rate_equation(curr_Ca, ourparams(41), ourparams(42) ,ourparams(43), ourparams(5),curr_I,curr_R_a); 

        %change in apmlifier
        d_S_dt      = curr_S * (curr_R_a - ourparams(51));
    end
    
    %% function which relevant only for integrated model
    function [d_gating_params_dt, d_V_dt, d_Ca_dt, d_I_dt, d_S_dt] = fullIntegratedModel(curr_gating_params,curr_V, curr_E_C, curr_Ca, curr_R_a, curr_S, curr_I, ourparams)
        %parameters for the main model function

        %avogadro number
        mol         =   6.022e23;

        %electron charge for computing the charge of Ca ion
        e_charge    =   1.602e-19;

        %AWA cell volume in Liters
        AWA_vol     =   1.5e-13;

        %capacitance of the cell in nano farad
        cap         =   1.5e-3;

        %contuctence of Calcium channels (egl-19) in nS 
        gCa         =   0.1;

        %calicum reversal potential in mV
        vCa         =   120.0;

        %calcium inactivation factor
        f           =   0.4;

        %conductence of potasium, SHK-1 channels in nS.
        gK1         =   1.5;

        %conductence of inward rectifyng current-dependent K channel in nS  
        gKir        =   0.8;

        %inward rectifying potassium channel saturated currnt. in pA
        I_0         =   5.0;

        %inward rectifying potassium channel smoothing factor
        beta        =   5.0;

        %conductence of slow inactivating potasium channels Kr in nS
        gKr         =   0.39;

        %conductence of first high-V potasium channel (slo-1 or slo-2) in nS
        gK2         =   1.0;

        %conductence of second high-V potasium channel (slo-1 or slo-2) in nS
        gK3         =   0.1;

        %conductence of slow activating potasium channels Kr_a in nS
        gKa         =   0.91;

        %potasium reversal potential in mV
        vK          =   -84.0;

        %conductence for leak current  in nS
        gL          =   0.25;

        %leak current reversal potential in mV
        vL          =   -65.0;  

        % parameters for all gating dynamics

        %SHK1 potassium channels decay time in mili seconds ("switching time")
        Tau_SHK1            =   30.0;

        %SHK1 potassium channels V0 ("activation voltage") in milivolt
        V0_SHK1             =   2;

        %SHK1 potassium channels delta V ("activation slope") in milivolt
        deltaV_SHK1         =   10;

        %Kr slowly inactivating potassium channels decay time in mili seconds ("switching time")
        Tau_Kr              =   1800.0;

        %Kr slowly inactivating potassium channels V0 ("activation voltage") in milivolt
        V0_Kr1              =   -42.0;

        %Kr slowly inactivating potassium channels delta V ("activation slope") in milivolt
        deltaV_Kr1          =   5.0;

        %Kr slowly inactivating potassium channels V0 ("inactivation voltage") in milivolt
        V0_Kr2              =   -42.0;

        %Kr slowly inactivating potassium channels delta V ("inactivation slope") in milivolt
        deltaV_Kr2          =   5.0;

        %First high-V potassium channels decay time in mili seconds ("switching time")
        Tau_K_highV1        =   1000.0;

        %First high-V potassium channels V0 ("activation voltage") in milivolt
        V0_K_highV1         =   13.0;

        %First high-V potassium channels delta V ("activation slope") in milivolt
        deltaV_K_highV1     =   20.0;

        %Second high-V potassium channels decay time in mili seconds ("switching time")
        Tau_K_highV2        =   1000.0;

        %Second high-V potassium channels V0 ("activation voltage") in milivolt
        V0_K_highV2         =   -25.0;

        %Second high-V potassium channels delta V ("activation slope") in milivolt
        deltaV_K_highV2     =   5.0;
        % model itself

        %calcium current in pA
        I_Ca = gCa*(curr_gating_params(1)+f*curr_gating_params(2))*(curr_V-vCa);  

        %leakage current  in pA
        I_leak = gL*(curr_V - vL);

        %inward rectifying potassium current in pA
        I_Kir = - gKir*(-I_0 + I_0*log(1 + exp(-1*(curr_V - vK - I_0)/beta)));

        %3 kind of pottasium channels with simple gating dynamics in pA
        I_K_123 = (curr_V - vK)*(gK1*curr_gating_params(3) + gK2*curr_gating_params(5) + gK3*curr_gating_params(7));

        % current of Kr in pA
        I_Kr = gKr*(1-curr_gating_params(4))*m_inf(curr_V,V0_Kr1,deltaV_Kr1)*(curr_V - vK);

        % current of Kra in pA
        I_Ka = gKa*curr_gating_params(6)*(curr_V - vK);

        %get the new rates (d_params_dt) for the next step
        d_gating_params_dt = nan(1,7);
        d_V_dt =  (curr_E_C - (I_Ca + I_leak + I_Kir + I_K_123 + I_Kr + I_Ka))/cap;
        [d_gating_params_dt(1), d_gating_params_dt(2)] = d_C12_dt(curr_V,curr_gating_params(1),curr_gating_params(2));
        d_gating_params_dt(3) = simple_gating_dynamics(curr_V,curr_gating_params(3),Tau_SHK1,V0_SHK1,deltaV_SHK1);
        d_gating_params_dt(4) = simple_gating_dynamics(curr_V,curr_gating_params(4),Tau_Kr,V0_Kr2,deltaV_Kr2);
        d_gating_params_dt(5) = simple_gating_dynamics(curr_V,curr_gating_params(5),Tau_K_highV1,V0_K_highV1,deltaV_K_highV1);
        d_gating_params_dt(6) = dm_Kra_dt(curr_V,curr_gating_params(6));
        d_gating_params_dt(7) = simple_gating_dynamics(curr_V,curr_gating_params(7),Tau_K_highV2,V0_K_highV2,deltaV_K_highV2);

        % tracking the calcium concentration derivative in M/ms. we convert from C/ms to M/ms
        d_Ca_dt = -(I_Ca*1e-15/(2*e_charge*mol*AWA_vol) + (1/ourparams(3))*(curr_Ca - ourparams(5)));

        %change in modification level of receptros based on calcium level in cell    
        d_I_dt = Inhibition_rate_equation(curr_Ca, ourparams(41), ourparams(42), ourparams(43),ourparams(5), curr_I, curr_R_a); 

        d_S_dt = curr_S * (curr_R_a - ourparams(51));
    end

    function [ mInf ] = m_inf(V,V0,deltaV )                                               
        %this is function mInf which is used for K channels and also Ca channels (for Ca there is product of 2 mInf)
        mInf =  0.5 * (1 + tanh((V-V0)/deltaV));
    end

    function [dm_dt] = simple_gating_dynamics(V,curr_m,Tau,V0,deltaV)
        %the change in gating parametes for the case of simple dynamicss
        mInf = m_inf(V,V0,deltaV);
        dm_dt = (mInf-curr_m)/Tau;
    end

    function [d_Ka_dt] = dm_Kra_dt (V,curr_Ka)
        %slow activation time in mili seconds ("Tau_Ka_1")
        Tau_Ka_1 = 27000;

        %activation voltage in mili volt ("V_Ka_1")
        V_Ka_1 = -52;

        %activation slope in mili volt ("deltaV_Ka_1")
        deltaV_Ka_1 = 20;

        %fast activation time in mili seconds ("Tau_Ka_2")
        Tau_Ka_2 = 3000;

        %activation voltage in mili volt ("V_Ka_2")
        V_Ka_2 = -32;

        %activation slope in mili volt ("deltaV_Ka_2")
        deltaV_Ka_2 = 4;

        %this is how it supposed to be according to the paper
            %Tau_Ka = Tau_Ka_1 + (Tau_Ka_2 - Tau_Ka_1)*m_inf(V,V_Ka_2,deltaV_Ka_2);
            %d_Ka_dt = simple_gating_dynamics(V,curr_Ka,Tau_Ka,V_Ka_1,deltaV_Ka_1);  

        %this is how it appears in the old simulation code
            Tau_Ka = Tau_Ka_1 + (Tau_Ka_2 - Tau_Ka_1)*0.5*(1+tanh(V-V_Ka_1)/deltaV_Ka_1);    
            d_Ka_dt = simple_gating_dynamics(V,curr_Ka,Tau_Ka,V_Ka_2,2);  
    end

    function [ d_C1_dt , d_C2_dt ] = d_C12_dt(V,curr_C1, curr_C2)
        %calcium channel steady state normalization
        C0 = 0.55;

        %Calcium channel inactivation time voltage scale ("V_C3")
        V_C3 = 23;                                                                       

        %Calcium channel inactivation time voltage slope ("deltaV_C3")
        deltaV_C3 = 24; 

        %Calcium channel activation voltage ("V_C1") 
        V_C1 = -21.6;    

        %Calcium channel activation slope ("deltaV_C1") 
        deltaV_C1 = 9.17;                                                                     

        %Calcium channel inactivation voltage ("V_C2") 
        V_C2 = 16.2;    

        %Calcium channel inactivation slope ("deltaV_C2") 
        deltaV_C2 = -16.1;  

        %calcium channel activation time scale ("Tau_C1")
        Tau_C1 = 1.0; 

        %calcium channel activation time scale ("Tau_C2")
        Tau_C2 = 160.0; 

        %this is how it appears in the paper
        %     Tau_C2_of_V = Tau_C2/cosh((V-V_C3)/deltaV_C3);
        %     d_C2_dt = (curr_C1 - curr_C2) / Tau_C2_of_V;
        % 
        %     Ca_inf = m_inf(V,V_C1,deltaV_C1);
        %     Ci_inf = (2/C0)*m_inf(V,V_C2,deltaV_C2);
        %     d_C1_dt = (Ca_inf*Ci_inf*(1-curr_C2) - curr_C1)/Tau_C1 - (curr_C1 - curr_C2)/Tau_C2_of_V;
        %this is how it appears in the old simulation code
            T2 = 80*(1/cosh((V-23)/(2*24)));                                                    
            mx = 3.612/4;
            m1 = m_inf(V,-21.6,9.17);
            m2 = m_inf(V,16.2,-16.1);
            d_C1_dt = (m1*m2/mx - curr_C1)/1 - m1*m2*curr_C2/(mx*1) - curr_C1/(2*T2) + curr_C2/(2*T2);
            d_C2_dt = (curr_C1-curr_C2)/(2*T2);   
    end
    %% other functions
    function  [L_curr] = InputLigandProfile(t_curr,gradType, baseline, ts, tf, lin_slope,amplitude,step_duration, break_duration)
        %this is the input Odorant concentration as function of time, in micro molar

        if ((t_curr < tf) && (t_curr > ts))
            switch gradType 
                case 1 %linear
                    L_curr      =   (t_curr - ts)*lin_slope + baseline;
                case 2 %step
                    L_curr      =   amplitude + baseline;
                case 3 %tanh. sigmoid wil start at 1% at ts and each 99.9 precent at tf
                    t_center    =   ts + (tf - ts)/2;
                    L_curr      =   baseline + amplitude*(1 + tanh(-3.8*(t_curr - t_center)/(ts - t_center)));
                case 4 %growing steps
                    growRatio   =   7;
                    %ligand_curr = baseline + amplitude*(growRatio^floor((t_curr - ts)/step_duration));
                    L_curr      =   baseline*(growRatio^(1+floor((t_curr - ts)/step_duration)));
                case 5 %steps - notice that the breaks duration is the same as steps duration
                    if     mod(floor((t_curr - ts)/step_duration),2) == 0, L_curr = baseline + amplitude;
                    elseif mod(floor((t_curr - ts)/step_duration),2) == 1, L_curr = baseline;
                    end
                case 6 %exponential
                    K           =   log(amplitude/0.1)/(tf-ts);
                    L_curr      =   baseline + 0.1*exp(K*(t_curr-ts));
                case 7 %2-steps
                    if ((t_curr-ts) < step_duration), L_curr = baseline + amplitude;
                    elseif ((t_curr-ts - step_duration) < break_duration), L_curr = baseline;
                    else, L_curr = baseline + amplitude;
                    end
            end
        else
            L_curr = baseline;
        end  
    end

    function [R_a] = Receptor_activation(curr_L, curr_I,C_I,C_L,L_0) 
        %the activation of the ceceptor based on the ligand conc. and the inhibition level 
        R_a = C_L*log10(curr_L/L_0) - C_I*curr_I;
        R_a = 1/(1+exp(-R_a));
    end

    function [d_I_dt] = Inhibition_rate_equation(curr_Ca, calcium_dependent_inhibition_coeff,...
                                                calcium_independent_inhibition_coeff,inhibition_decay_time_const , Ca_0,curr_I,curr_R_a)
        d_I_dt = calcium_dependent_inhibition_coeff * (curr_Ca - Ca_0)*(curr_R_a) + ...
                 calcium_independent_inhibition_coeff*(curr_R_a) - (1/inhibition_decay_time_const)*curr_I*(1-(curr_R_a));
    end
end