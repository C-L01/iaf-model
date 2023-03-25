%% Function Name: solve_iaf
%
% Description: Solve a supplied Integrate-and-Fire model for N neurons
% Uses ode23s, since spikes likely lead to stiffness.
%
% Inputs:
%     P: parameter struct, containing all required parameters
% Returns:
%     sol: ode solution struct
%
% $Revision: R2022b$
%---------------------------------------------------------

function [sol] = solve_iaf(P)
    
    % Solve until the first spike (or maxTime if there are none)
    sol = ode23s(P.ode, [P.tStart P.maxTime], P.u0, P.Options);

    % Get the values at time of ode termination (could be maxTime)
    tEnd = sol.x(end);
    uEnd = sol.y(:,end);

    if P.maxTime <= tEnd
        disp("Termination condition: Max Time")
        return
    end

    for spikecnt=1:(P.maxSpikes-1)

%         if rem(spikecnt,10) == 0
%             fprintf("Spike %i\n", spikecnt)
%         end

        %% Perform post-spike potential updates

        % Find the neuron(s) that fired
        firedNeurons = (uEnd >= P.theta_reset);  % logical array
        
        % Add presynaptic spikes to potentials, based on weights
        u0 = uEnd + P.W * firedNeurons;
        
        % Reset fired neurons, explicitly done after the spike arrivals
        % to avoid immediate refiring
        u0(firedNeurons) = P.u_r;
        

        %% Restart the ode solver with updated initial condition
        sol = odextend(sol,[],P.maxTime,u0);

        % Update values at time of ode termination (could be maxTime)
        tEnd = sol.x(end);
        uEnd = sol.y(:,end);

        % Check termination condition
        if P.maxTime <= tEnd
            disp("Termination condition: Max Time")
            return
        end
    end
    
    % Loop exited without reaching maxTime
    disp("Termination condition: Max Spikes")
end

