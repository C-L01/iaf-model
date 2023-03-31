%% Script Name: run_iaf
%
% Description: Use solve_iaf with parameters from set_parameters to run
% a model of N interacting neurons, and plot potentials,
% population activity and potential density. Optionally save these
% visualizations.
%
% $Revision: R2022b$
%---------------------------------------------------------

clear; close all;

%% Run model

% Set seed (needed when i.c.'s are random)
rng(0);

% Get parameters
P = set_parameters();

% Solve ode
sol = solve_iaf(P);

% Whether to compute density evolution (somewhat expensive)
density = false;

% Whether to save the produced figures
save = false;


%% Extract and prepare the data for plotting

% Get time and potential values
% dt = 1e-8;
% time = linspace(P.tStart,P.maxTime,(P.maxTime - P.tStart)/dt);
% uSol = deval(sol,time);

% Remove duplicates to avoid double counting spikes
% The kept indices are the last (second) ones, making everything Càdlàg
[time,keptIndices,~] = unique(sol.x,'last');
uSol = sol.y(:,keptIndices);

% Get spike times
spikeTimes = sol.xe;                % could be multiple spikes at one time
numSpikeTimes = length(spikeTimes);
% Because of Càdlàgness, we check for equality to the reset value
% to detect all spikes at a spike time.
% Note that, theoretically, this could detect spikes that are not spikes
% if a neuron happens to be exactly at V_R at a spike time.
% TODO: check that at t_{i-1} the potential is close to V_F
firingNeurons = (uSol(:,ismember(time,spikeTimes)) == P.V_R);

% Approximate population activity A(t)
spikeCnt = arrayfun(@(t) ismember(t,spikeTimes)...
            * sum(sum(firingNeurons(:,t==spikeTimes))),time);
timeWindow = P.tau*min(5,(P.maxTime - P.tStart)/10);   % average activity time window
A = movmean(spikeCnt,timeWindow,'SamplePoints',time) / (P.N*timeWindow);

if density
    % Approximate density evolution p(t,u)
    h = figure;
    axis([P.V_R P.V_F 0 1]);    % fix axis limits
    xlabel("$u$ (V)",'Interpreter','latex','FontSize',16)
    ylabel("Density",'FontSize',16)
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    h.Visible = 'off';
    
    timeBinSize = 1;                 % the time length to group together
    timeBins = P.tStart:timeBinSize:P.maxTime;
    timeBinInc = find(diff(discretize(time,timeBins))); % increment indices
    timeBinInc = [1 timeBinInc];                        % start at time index 1
    
    loops = numel(timeBinInc)-1;
    M(loops) = struct('cdata',[],'colormap',[]);
    for j=1:loops
    %     fprintf("Loop %i\n",j)
        histogram(uSol(:,timeBinInc(j):timeBinInc(j+1)),'Normalization','pdf')
        drawnow
        M(j) = getframe(h);
    end
end

% Estimate variance evolution
uMean = sum(uSol,1) / P.N;
uVar = sum((uMean - uSol).^2,1) / P.N;  % or P.N - 1 to estimate actual var


if P.N <= 5
    % Print firing times
    disp("Firing times:")
    fmt = ['\t\t\t\t\t' repmat('%5.2f ',1,numSpikeTimes-1) '%5.2f\n\n'];
    fprintf(fmt,spikeTimes)
    
    disp("Firing neurons:")
    fmt = ['\t\t\t\t\t' repmat('%5i ',1,numSpikeTimes-1) '%5i\n'];
    fprintf(fmt,firingNeurons')
end


%% Plot results

% Separate potential line plots over time
if P.N <= 20
    f1 = figure;
    plot(time,uSol,'LineWidth',2.0)
    ylim([P.V_R P.V_F])
    xlabel("$t$ (s)",'Interpreter','latex','FontSize',16)
    ylabel("Potential (V)",'FontSize',16)
    title(sprintf("Potential over time for %i coupled neurons",P.N),'FontSize',20)
    legend(sprintfc("$i=%i$",1:P.N),'Interpreter','latex')
    legend('boxoff')
end


% Activity
f2 = figure;
plot(time,A,'LineWidth',2.0)
xlabel("$t$ (s)",'Interpreter','latex','FontSize',16)
ylabel("$A(t)$ (\#/s)",'Interpreter','latex','FontSize',16)
title(sprintf("Population activity for %i coupled neurons",P.N),...
       'FontSize',20)


% Variance
f3 = figure;
plot(time,uVar,'LineWidth',2.0)
% ylim([0 1])
xlabel("$t$ (s)",'Interpreter','latex','FontSize',16)
ylabel("$V(t)$ (V$^2$)",'Interpreter','latex','FontSize',16)  % 2 V's is confusing
title(sprintf("Variance evolution for %i coupled neurons",P.N),...
       'FontSize',20)



% Play potential density movie
if density
    figure;
    movie(gcf,M,1,3);
end

%% Save figures

if save

    % Set background color
    % f = gcf;
    % ax = gca;
    % set(f,'color',[241 239 239]/255);
    % set(ax,'color',[241 239 239]/255);
    
    paramInfo = sprintf('N=%i_w0=%.1f_I0=%.2f_r=%i',P.N,P.w0,P.I(0),P.r);
    
    % Potential figure
    if exist('f1','var')
        exportgraphics(f1,strcat('u_',paramInfo,'.pdf'),'ContentType', ...
                        'vector','BackgroundColor','none')
    end
    
    % Activity figure
    exportgraphics(f2,strcat('A_',paramInfo,'.pdf'),'ContentType', ...
                    'vector','BackgroundColor','none')
    
    % Variance figure
    exportgraphics(f3,strcat('V_',paramInfo,'.pdf'),'ContentType', ...
                    'vector','BackgroundColor','none')
    
    % Potential density evolution movie
    v = VideoWriter(strcat('p(u)_',paramInfo,'.avi'),'Uncompressed AVI');
    v.FrameRate = 3;
    open(v)
    writeVideo(v,M)
    close(v)
    disp("Finished saving files")
end
