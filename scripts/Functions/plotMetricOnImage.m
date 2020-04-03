function f = plotMetricOnImage(AdjacencyMatrix, CellLocations, metric, Masks)

CC = bwconncomp(Masks);
nNodes = size(AdjacencyMatrix, 1);
colors = parula(numel(unique(metric)));
cmap = zeros(nNodes, 3);
L = labelmatrix(CC);
UniqueMetricValues = unique(metric);

for i = 1:numel(UniqueMetricValues)
    index = ismember(metric, UniqueMetricValues(i));
    cmap(index, :) = repmat(colors(i, :), size(cmap(index, :), 1), 1);
end

metricHeatmap = label2rgb(L, cmap);

f = figure('Visible', 'Off'); set(gcf, 'color', 'w');
imshow(metricHeatmap);
hold on;
gplot(AdjacencyMatrix, CellLocations, 'r');
c = colorbar;
c.Ticks = [0 1];
c.TickLabels = [min(metric), max(metric)];
c.FontSize = 18;

end