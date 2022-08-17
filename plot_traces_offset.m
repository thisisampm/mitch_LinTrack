function [ax, lines] = plot_traces_offset(traces_mtx)
% Take traces and plot them with an offset 

C = traces_mtx; 

n_neurons = size(C,1); 

offset = linspace(0,n_neurons/2,n_neurons); % Create a vector to add to the traces causing them to be offset by a 3rd when plotting
C_offset = C + offset';

%axis(hAx_trace,[0 size(C_offset,1) 0 max(max(C_offset))],'padded');
figure;
ax = axes;
lines = plot(ax, C_offset');

