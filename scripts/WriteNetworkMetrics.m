
function [] = WriteNetworkMetrics(OutputPath, AdjacencyMatrices, GlobalMetrics, GlobalMetricNames, LocalMetrics, LocalMetricNames)

nGlobalMetrics = numel(GlobalMetricNames);
nLocalMetrics = numel(LocalMetricNames);
nImages = size(GlobalMetrics, 1);

%% Print Global Metrics

% opening file for writing global metrics
WritePath = strcat(OutputPath, filesep, 'GlobalMetrics.csv');

fid = fopen(WritePath, 'w');

% printing metric names
fprintf(fid, ',');

for i = 1:nGlobalMetrics
    fprintf(fid, '%s,', GlobalMetricNames{i});
end

fprintf(fid, '\n');
formatSpec = '%f,';

% printing image names and metrics
for i = 1:nImages
    currentMaskName = GlobalMetrics{i, 1};
    GlobalMetricsCurrent = GlobalMetrics{i, 2};
    imageName = sprintf('%s', currentMaskName);
    fprintf(fid, '%s,', imageName);
    fprintf(fid, formatSpec, GlobalMetricsCurrent(:));
    fprintf(fid, '\n');
end

fclose(fid);

%% Print Local Metrics

for k = 1:nImages
    
    currentMaskName = LocalMetrics{k, 1};
    
    % opening file for writing local metrics
    WritePath = strcat(OutputPath, filesep, 'LocalMetrics_', currentMaskName, '.csv');
    writetable(LocalMetrics{k, 2}, WritePath);
end

%% Print Adjcacency Matrix

for k = 1:nImages
    
    currentMaskName = AdjacencyMatrices{k, 1};
    
    % opening file for writing local metrics
    WritePath = strcat(OutputPath, filesep, 'AdjacencyMatrix_', currentMaskName, '.csv');
    csvwrite(WritePath, AdjacencyMatrices{k, 2});
end