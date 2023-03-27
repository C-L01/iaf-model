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
spikeTimes = sol.xe;
numSpikes = length(spikeTimes);
firingNeurons = sol.ye >= P.V_F;

% Approximate population activity A(t)
spikeOccurs = arrayfun(@(t) ismember(t,spikeTimes),time);
timeWindow = 1e-1*P.tau*(P.maxTime - P.tStart);   % average activity time window
A = movmean(spikeOccurs,timeWindow,'SamplePoints',time) / (P.N*timeWindow);

density = false;

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

if P.N <= 5
    % Print firing times
    disp("Firing times:")
    fmt = ['\t\t\t\t\t' repmat('%5.2f ',1,numSpikes-1) '%5.2f\n\n'];
    fprintf(fmt,spikeTimes)
    
    disp("Firing neurons:")
    fmt = ['\t\t\t\t\t' repmat('%5i ',1,numSpikes-1) '%5i\n'];
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


% Spikes
f2 = figure;
plot(time,A,'LineWidth',2.0)
% ylim([0 P.V_F])
xlabel("$t$ (s)",'Interpreter','latex','FontSize',16)
ylabel("$A(t)$ (\#/s)",'Interpreter','latex','FontSize',16)
title(sprintf("Population activity for %i coupled neurons",P.N),...
       'FontSize',20)

% Play potential density movie
if density
    figure;
    movie(gcf,M,1,3);
end

%% Save figures

save = false;

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
    
    % Potential density evolution movie
    v = VideoWriter(strcat('p(u)_',paramInfo,'.avi'),'Uncompressed AVI');
    v.FrameRate = 3;
    open(v)
    writeVideo(v,M)
    close(v)
    disp("Finished saving files")
end
