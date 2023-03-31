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

    P.N = 100;                         % number of neurons
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
%     I = @(t) exp(t/10);

    tPulse = 1;                     %(s) time of pulse
    dt = 1e-2;                      %(s) time interval of pulse input
    dV = 5.5;                         %(V) voltage increase due to pulse
%     I = @(t) (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt));


%     I = @(t) 15*exp(-t);

    P.I = I;

    % Note: for certain odes elementwise operations may be needed
    % Exponential
    P.ode = @(t,u) (-(u-u_rest) + delta_T * exp((u-theta_rh) / delta_T)...
                  + R*I(t))/tau;
    % Leaky
%     P.ode = @(t,u) (-(u-u_rest) + R*I(t))/tau;



    % Symbolic ode (temporary manual implementation)
%     syms odeSym tSym uSym;
%     odeSym(tSym,uSym) = (-(uSym-0) + 1 * exp((uSym-5) / 1)...
%                   + 1*4.75)/1;
%     jacobian(odeSym);

    % Precompute Jacobian, assuming constant I
%     P.J = @(t,u) [0, -1 + exp((u - theta_rh) / delta_T)] / tau;


%% Reset conditions
% Conditions that determine what happens when a spike occurs
    P.V_F = 15;             %(V) Firing threshold
    P.V_R = -10;           %(V) Potential after reset

    % Event function for ode solver to stop at spikes
    function [value,isterminal,direction] = ResetEvent(~,u)
    
        value = min(P.V_F - u);
        
        % Stop if the reset potential is reached
        isterminal = true;
        
        % Should only be approaching V_F from below
        direction = -1;
    end


%% Exit conditions
% Run until either the maximum number of spikes or the maximum time 
    P.maxSpikes = 100*P.N;    %(#) Number of spikes before exit
    P.maxTime = 20;        %(s) Maximum simulation time


%% Initial conditions

    P.tStart = 0;                   %(s) Starting time

    %(V) Starting potentials
    P.r = 5;                        %(V) (Expected) radius of the initial potentials
    P.u0 = linspace(P.r,-P.r,P.N);
%     P.u0 = 2*P.r*(rand(P.N,1)-0.5);   % ~Unif[-r,r]
%     P.u0 = P.r*randn(P.N,1);          % ~N(0,r^2)


%% ode23s options struct

    P.Options = odeset('RelTol', 1e-4,...
                       'AbsTol', 1e-6,...
                       'Events', @ResetEvent,...
                       'Vectorized', true,...
                       'InitialStep', 1e-3,...
                       'MaxStep', 1e-2,...
                       'JPattern', speye(P.N));
%                        'Jacobian', P.J);
    
    % Initial step is needed in case there is a pulse (using heaviside)
    % near the start of the time interval

end

