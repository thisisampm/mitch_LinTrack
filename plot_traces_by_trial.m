load('C:\Users\mitch\OneDrive - University of Toronto\PhD\Experiments\LinearTrack\Processed\152-2\AMPM_Analyzed\2021-04-16.mat')

close all;
% Generate a plot of the behavioural trials.
%%
nTrials = numel(unique(behavior_mtx(~isnan(behavior_mtx(:,4)),4)));

fig1 = figure;
plot(behavior_mtx(:,1)./60, behavior_mtx(:,2))
ylim([-0.05 1.05])

hold on
for i = 1:nTrials
    ind = find(behavior_mtx(:,4) == i);
    plot(behavior_mtx(ind,1)./60 ,behavior_mtx(ind,2),'LineWidth',2);
end    
xlabel('Time (min)')
ylabel('Relative Track Position')
hold off
%%
figure;
fig2 = subplot(2,1,1);
plot(fig2, behavior_mtx(:,1)./60,behavior_mtx(:,2))
hold(fig2,'on')
ylim([-0.05 1.05])
% Generate a plot of Ca activity per trial.
fig3 = subplot(2,1,2);
hold(fig3,'on')
col = colororder;
time = cell(nTrials,1);
Ca = cell(nTrials,1);

for i = 1:10
        ind = find(behavior_mtx(:,4) == i);
        col = [rand,rand,rand];
        plot(fig2,behavior_mtx(ind,1)./60,behavior_mtx(ind,2),'LineWidth',2,'Color',col);

        time{i} = (behavior_mtx(ind,1) - behavior_mtx(ind(1),1))/max(behavior_mtx(ind,1) - behavior_mtx(ind(1),1));
        %time = discretize(time,45);
        Ca{i} = (traces(24,ind,2) - min(traces(24,ind,2) ))./max(traces(24,ind,2));
        plot(fig3,time{i},Ca{i},'Color',col)
        pause(1/2)
end
%%

j = 1:2:21;
for i = 1:11
    ind = find(behavior_mtx(:,4) == j(i));
    col = [rand,rand,rand];
    time{i} = (behavior_mtx(ind,1) - behavior_mtx(ind(1),1))/max(behavior_mtx(ind,1) - behavior_mtx(ind(1),1));
    %time = discretize(time,45);
    Ca{i} = (traces(24,ind,2) - min(traces(24,ind,2) ))./max(traces(24,ind,2));
    subplot(8,2,i)
    plot(time{i},Ca{i},'Color',col)
    title(['Trial ', num2str(i)])
    set(gca,'XTickLabel',[])
end
subplot(8,2,[11 12 13 14])
plot(time{8},Ca{8},'Color','k')
title('Average Over Trials')
xlabel('Relative Position')
ylabel('Normalized dF/F')

%%
time_bin = discretize(time{8},45);
Ca11 = Ca{8};
Ca_bin = nan(size(unique(time_bin)));

for i = 1:numel(unique(time_bin))
    
    Ca_bin(i) = mean(Ca11(time_bin == i));
end

subplot(8,2,[15 16])
imagesc(Ca_bin')
set(gca,'YTickLabels',[])
axis equal
axis tight
xlabel('Spatial Bin')




    
    
    
    
    
    