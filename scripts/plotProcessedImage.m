%% function to display image processing results

function [] = plotProcessedImage(prImagePath, Masks, AdjacencyMatrix, CellLocations, MaskBoundaries, LocalMetrics)

% plot graph overlaid on image
f = figure('Visible', 'Off');
imshow(Masks, 'Border', 'Tight'); hold on;

if ~isempty(AdjacencyMatrix)
    gplot(AdjacencyMatrix, CellLocations, 'y');
    
end

for k = 1:length(MaskBoundaries)
    boundary_i = MaskBoundaries{k};
    plot(boundary_i(:, 2), boundary_i(:, 1), 'r');
    text(CellLocations(k, 1), CellLocations(k, 2), num2str(k), 'FontSize', 12, 'Color',  'r');
end

set(findall(gcf,'type','line'), 'LineWidth', 2);

saveas(f, prImagePath);
close(f);

% plot local metrics as heatmaps overlaid on image
localMetrics2Plot = {'degree', 'avgeNeighborDegree', 'clusteringCoefficient', 'closenessCentrality', 'betweennessCentrality'};
nLocalMetrics2Plot = numel(localMetrics2Plot);
for i = 1:nLocalMetrics2Plot
    currentLocalMetric = LocalMetrics.(localMetrics2Plot{i});
    f = plotMetricOnImage(AdjacencyMatrix, CellLocations, currentLocalMetric, Masks);
    saveas(f, strrep(prImagePath, '.tif', strcat('_', localMetrics2Plot{i}, '.tif')));
    close(f);
end