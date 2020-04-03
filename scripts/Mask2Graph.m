
function [errorReport, imageProcessingParams] = Mask2Graph(MaskInputPath, OutputPath, thresholdType, threshold, maskLayerNumber, errorReport)

s = createImageProcessingParams(0, 0, 0, 0, 0, 0, 0, 0);
imageProcessingParams = repmat(s, 0);

flag = false;

%% starting parallel pool
myCluster = parcluster('local');
myCluster.JobStorageLocation = pwd;
poolObj = parpool(myCluster);

%% Loading parameters and initializing variables

% reading all files in folder specified by InputPath
MaskReadPath = strcat(MaskInputPath, filesep, '*');
d = dir(MaskReadPath);
fnme = {d.name};

% accepting only files and not folders
maskIdx = cellfun(@isfolder, fnme);
fnmeMasks = fnme(~maskIdx);

nMaskFiles = numel(fnmeMasks);

GlobalMetrics = {}; % global metric cell array
LocalMetrics = {}; % local metric cell array
AdjacencyMatrices = {}; % cell array to hold adjacency matrices

%% processing all images in folder
ctr = 1;
for i = 1:nMaskFiles
    currentMaskFile = strcat(MaskInputPath, filesep, fnmeMasks{i});
    currentMaskBaseName = createBaseName(fnmeMasks{i});
    
    % finding layer number and base file name for already segmented images
    idxLayer = strfind(currentMaskBaseName, 'layer');
    idxMASKS = strfind(currentMaskBaseName, 'MASKS');
    
    s = '';
    if isempty(idxLayer)
        layerNumberParam = '';
    else
        layerNumberParam = currentMaskBaseName(idxLayer+5);
        s = strcat('-layer', layerNumberParam);
        
        if ~isempty(idxMASKS)
            s = strcat(s, '-MASKS');
        end
    end

    originalImageName = strrep(currentMaskBaseName, s, '');
    
    % making sure currentMaskFile is a readable image file
    try
        info = imfinfo(currentMaskFile);
    catch
        errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' is not a readable image file'], 1);
        continue;
    end
    
    imageDimensions = [info.Width info.Height];
    nLayers = numel(info); % number of layers in image

    % note that maskLayerNumber will be an empty vector if segmentation has
    % been requested
    if isempty(maskLayerNumber)
        maskLayerNumber = 1;
    end
    
    for k = 1:numel(maskLayerNumber) % looping over number of layers in image to be processed
        if maskLayerNumber == 1
            CurrentMask = imread(currentMaskFile);
            layerSuffix = '';
        else
            if nLayers<maskLayerNumber(k)
                errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' does not have ', num2str(maskLayerNumber(k)), ' layers'], 1);
                continue;
            end
            
            CurrentMask = imread(currentMaskFile, maskLayerNumber(k));
            layerSuffix = strcat('layer', num2str(maskLayerNumber(k)));
            layerNumberParam = num2str(maskLayerNumber(k));
        end
        
        flag = true; % this means we have usable input
            
        % if image is RGB, converting to grayscale
        if size(CurrentMask, 3) > 1
            CurrentMask = rgb2gray(CurrentMask(:, :, 1:3));
            errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' is an RGB image: converting to grayscale'], 1);
        end
        
        % making sure image is logical
        if ~islogical(CurrentMask)
            CurrentMask = CurrentMask > 0;
            errorReport(end+1) = makeErrorStruct(['file ', fnmeMasks{i}, ' is not a binary image: converting to binary'], 1);
        end
        
        fprintf('\t creating graph representation for image %s\n', fnmeMasks{i});

        currentMaskName = strcat(currentMaskBaseName, layerSuffix);
        GlobalMetrics{ctr, 1} = currentMaskName;
        LocalMetrics{ctr, 1} = currentMaskName;
        AdjacencyMatrices{ctr, 1} = currentMaskName;
        
        % calling CalculateNetworkMetrics to compute graph-based metrics
        [GlobalMetrics{ctr, 2}, LocalMetrics{ctr, 2}, AdjacencyMatrices{ctr, 2}, MaskBoundaries, ...
            CellLocations, GlobalMetricNames, LocalMetricNames, processParams]...
            = CalculateNetworkMetrics(CurrentMask, threshold, thresholdType);
        
        if processParams.processCompleted
            % call function to plot processed image figure
            prImagePath = strcat(OutputPath, filesep, currentMaskName, '-PROCESSED.tif');
            plotProcessedImage(prImagePath, CurrentMask, AdjacencyMatrices{ctr, 2}, CellLocations, MaskBoundaries, LocalMetrics{ctr, 2});
        else
            errorReport(end+1) = makeErrorStruct(['processing for file ', fnmeMasks{i}, ' exceeds computing capacity; graph not created'], 1);
        end
        
        ctr = ctr + 1;
        imageProcessingParams(end+1) = createImageProcessingParams(originalImageName, ...
            layerNumberParam, imageDimensions, processParams.nNodes, processParams.nEdges, ...
            processParams.totalWorkspaceSize, processParams.processingTime, processParams.processCompleted);
    end
end

%% writing metrics to file

fprintf('\t\t\twriting metrics to file\n');

if flag
    WriteNetworkMetrics(OutputPath, AdjacencyMatrices, GlobalMetrics, GlobalMetricNames, LocalMetrics, LocalMetricNames);
end

%% cleaning up
delete(poolObj);

end