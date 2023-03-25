%% Function Name: set_parameters
%
% Description: Set parameters for solve_iaf.
% This file is to be modified manually.
%
% Inputs: None
%
% Returns:
%     P: struct containing all parameters solve_iaf needs
%
% $Revision: R2022b$
%---------------------------------------------------------

function P = set_parameters()

%% Network settings

    P.N = 20;                         % number of neurons
    P.w0 = 1;
    P.W = P.w0*ones(P.N,P.N) / P.N;     % synaptic weights
%     P.W = P.w0*hilb(P.N) / P.N;     % synaptic weights

%% ODE
    tau = 1;                        %(s) Time constant
    P.tau = tau;
    u_rest = 0;                     %(V) Resting potential
    delta_T = 0.5;                  %( ) Sharpness parameter
    theta_rh = 5;                   %(V) Rheobase threshold
    R = 1;                          %(Ω) Resistance


    %(A) External stimulus
    I = @(t) 4.75;

    tPulse = 1;                     %(s) time of pulse
    dt = 1e-2;                      %(s) time interval of pulse input
    dV = 5.5;                         %(V) voltage increase due to pulse
%     I = @(t) (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt));


%     I = @(t) 15*exp(-t);

    P.I = I;

    % Note: for certain odes elementwise operations may be needed
    P.ode = @(t,u) (-(u-u_rest) + delta_T * exp((u-theta_rh) / delta_T)...
                  + R*I(t))/tau;

    % Nondimensionalized parameters
    % lam = c1/c0; ome1 = w1/c0; ome2 = w2/c0;
    % rho0 = r0/(c0*duration); rho1 = r1/(c0*duration);
    % gam0 = delta0/(c0*duration); gam1 = delta1/(c0*duration); mu = M/c0;
    
    % Nondimensionalized initial conditions
    % tau_span = [tstart*c0 tend*c0];
    % s0 = S0/N; s1 = S1/N; s2 = S2/N; i1 = I1/N; i2 = I2/N;
    % y0 = [s0; s1; s2; i1; i2];


%% Reset conditions
% Conditions that determine what happens when a spike occurs
    P.theta_reset = 15;    %(V) Reset threshold
    P.u_r = -10;           %(V) Potential after reset

    % Event function for ode solver to stop at spikes
    function [value,isterminal,direction] = ResetEvent(~,u)
    
        value = min(P.theta_reset - u);
        
        % Stop if the reset potential is reached
        isterminal = true;
        
        % Should only be approaching theta_reset from below
        direction = -1;
    end


%% Exit conditions
% Run until either the maximum number of spikes or the maximum time 
    P.maxSpikes = 100*P.N;    %(#) Number of spikes before exit
    P.maxTime = 40;        %(s) Maximum simulation time


%% Initial conditions

    P.tStart = 0;                   %(s) Starting time

    %(V) Starting potentials
    P.r = 5;                        %(V) (Expected) radius of the initial potentials
    P.u0 = linspace(P.r,-P.r,P.N);
%     P.u0 = 2*P.r*(rand(P.N,1)-0.5);
    %     P.u0 = P.r*randn(P.N,1);


%% ode23s options struct
    P.Options = odeset('RelTol', 1e-4,...
                       'AbsTol', 1e-8,...
                       'Events', @ResetEvent,...
                       'Vectorized', true,...
                       'InitialStep', 1e-3,...
                       'MaxStep', 1e-2);
    
    % Initial step is needed in case there is a pulse (using heaviside)
    % near the start of the time interval

end

