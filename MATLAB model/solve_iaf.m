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

        % uRestart is the new starting vector
        uRestart = uEnd;

        % Find the neuron(s) that fired
        firedNeurons = (uRestart >= P.V_F);  % logical array
        unfiredNeurons = ones(P.N,1);

        while any(firedNeurons)
            % Update the list of unfired neurons
            unfiredNeurons = unfiredNeurons & ~firedNeurons;

            % Increment potential of unfired neurons based on weights
            uRestart(unfiredNeurons) = uRestart(unfiredNeurons)...
                                + P.W(unfiredNeurons,:) * firedNeurons;
            
            % Reset fired neurons
            uRestart(firedNeurons) = P.V_R;

            % The arrival of spikes might have caused new neurons to fire
            firedNeurons = (uRestart >= P.V_F);  % logical array
        end
        

        %% Restart the ode solver with updated initial condition
        sol = odextend(sol,[],P.maxTime,uRestart);

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

