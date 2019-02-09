%
% This script processes the results of the multiple swingby interplanetary
% transfer optimization run by the multiObjectiveEarthMarsTransferExample.cpp
% Tudat/Pagmo2 example.
%

set(0, 'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultTextInterpreter','latex');

clc
clear all
close all

saveFolder = '../SimulationOutput/';

% Define plot settings
numberOfOptimizers = 1;
numberOfGenerations = 100;
plotInterval = 98;

% Create data storage containers
fitness = cell(numberOfOptimizers,numberOfGenerations);
population = cell(numberOfOptimizers,numberOfGenerations);

% Iterate over requested generations, and plot results
counter = 1;
i = 1;



for j = 1:plotInterval:numberOfGenerations
    
    % Load data for current generation
    fitness{i,counter} =  load(strcat(saveFolder,'fitness_mo_mga_EVEMJ_',num2str(j),'.dat'));
    population{i,counter} =  load(strcat(saveFolder,'population_mo_mga_EVEMJ_',num2str(j),'.dat'));
    
    % Plot data, colored by departure date
    figure(1)
    subplot(4,4,counter)
    scatter(fitness{i,counter}(:,1)/1000,fitness{i,counter}(:,2), 10,'*')
    if(rem(counter,4)==1)
        ylabel('Travel time [days]')
    end
    if(counter > 12 )
        xlabel('\Delta V [km/s]')
    end
    title(strcat('Iteration ',{' '},num2str(j)));
    grid on
    
    % Plot data, colored by leg duration (per leg)
    for k=1:size(population{i,counter},2)
        figure(k+1)
        subplot(4,4,counter)
        if( k == 1 )
            scatter(fitness{i,counter}(:,1),fitness{i,counter}(:,2)/1000, 10,(population{i,counter}(:,k)-2451545)/365,'*')
        else
            scatter(fitness{i,counter}(:,1),fitness{i,counter}(:,2)/1000, 10,population{i,counter}(:,k),'*')
        end
        if(rem(counter,4)==1)
            ylabel('Travel time [days]')
        end
        if(counter > 12 )
            xlabel('Delta V [km/s]')
        end
        title(strcat('Gen.=',{' '},num2str(j)));
        colorbar
        grid on
    end
    [Minimum,Index]=min(fitness{i,counter}(:,1))
    bestPopulation = population{i,counter}(Index,:)
    counter = counter + 1;
    grid on
    
    
end

% Resize figure windows
for k=1:size(population{1,1},2)+1
    set(figure(k), 'Units', 'normalized', 'Position', [0,0,0.75 0.75]);
    set(figure(k),'PaperUnits','centimeters','PaperPosition',[0 0 45 30]);
    set(figure(k),'PaperPositionMode','auto');
end

% Add figure titles
for k=1:size(population{i,counter},2)
    figure(k+1)
    if( k== 1)
        suptitle('Color scale: departure date [years since J2000]')
    elseif( k== 2)
        suptitle('Color scale: leg 1 (E-V) travel time [days]')
    elseif( k== 3)
        suptitle('Color scale: leg 2 (V-E) travel time [days]')
    elseif( k== 4)
        suptitle('Color scale: leg 3 (E-E) travel time [days]')
    elseif( k== 5)
        suptitle('Color scale: leg 4 (E-J) travel time [days]')
    end
end
%%
pause(0.1)
figure(1)
saveas(gcf,strcat('swingbyOptimizationEVEEJ_paretoFront'),'png');

for k=1:size(population{1,1},2)
    pause(0.1)
    figure(k+1)
    saveas(gcf,strcat('swingbyOptimizationEVEEJ_parameter',num2str(k)),'png');
end
