%% Create colorbars for figures.

% First - open the figure to create an isolated colorbar for

% Get the range from figure
lims = caxis;

% Make a new figure and recreate desired colorbar
figure;
colorbar;
caxis(gca, lims);

% Make the figure transparent
set(gca, 'color', 'none');
set(gca, 'FontSize', 16);

% Save out colorbar
print('be-cf-colorbar','-dpdf')

close all