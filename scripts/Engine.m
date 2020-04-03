function [errorReport, imageProcessingParams] = Engine(InputPath, OutputPath, segmentationStruct, Mask2GraphStruct)

s = makeErrorStruct('', '');
errorReport = repmat(s, 0);

segmentationRequested = segmentationStruct.segmentationRequested;
segmentationType = segmentationStruct.segmentationType;
segmentationLayerNumber = segmentationStruct.segmentationLayerNumber;

Mask2GraphRequested = Mask2GraphStruct.Mask2GraphRequested;
thresholdType = Mask2GraphStruct.thresholdType;
threshold = Mask2GraphStruct.threshold;
maskLayerNumber = Mask2GraphStruct.maskLayerNumber;

if segmentationRequested
    
    % make new folder to store masks
    MaskInputPath = 'MASKS';
    mkdir(MaskInputPath);
    
    % reading all files in folder specified by InputPath
    ImageReadPath = strcat(InputPath, filesep, '*');
    d = dir(ImageReadPath);
    fnme = {d.name};
    
    % accepting only files and not folders
    imageIdx = cellfun(@isdir, fnme);
    fnmeImages = fnme(~imageIdx);
    
    nImages = numel(fnmeImages);
    SingleCellMetrics = {}; % cell array containing cell shape and intensity metrics
    SingleCellMetricNames = {}; % cell array containing single-cell metric names
    ctr = 1;
    for i = 1:nImages
        currentImageFile = strcat(InputPath, filesep, fnmeImages{i});
        
        % making sure currentImageFile is a readable image file
        try
            info = imfinfo(currentImageFile);
        catch
            errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' is not a readable image file; skipping this file'], 1);
            continue;
        end
        
        % using segmentationLayerNumber to extract relevant images for
        % segmentation
        nLayers = numel(info); % total number of layers in image
        for k = 1:numel(segmentationLayerNumber)
            if nLayers == 1 || all(segmentationLayerNumber == 1) % nLayers will be 1 if image is a non-layered image
                CurrentImage = imread(currentImageFile);
                layerSuffix = '';
                
                if nLayers<segmentationLayerNumber(k)
                    errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' does not have ', num2str(segmentationLayerNumber(k)), ' layers; skipping this file'], 1);
                    continue;
                end
            else
                if nLayers<segmentationLayerNumber(k)
                    errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' does not have ', num2str(segmentationLayerNumber(k)), ' layers; skipping this file'], 1);
                    continue;
                end
                
                CurrentImage = imread(currentImageFile, segmentationLayerNumber(k));
                layerSuffix = strcat('-layer', num2str(k));
            end

            % if image is RGB, converting to grayscale
            if size(CurrentImage, 3) >= 3
                CurrentImage = rgb2gray(CurrentImage(:, :, 1:3));
                errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' is an RGB image: converting to grayscale'], 1);
            else 
            % Check for unusual situation of a two channel image
                if size(CurrentImage, 3) == 2
                    CurrentImage = CurrentImage(:, :, 1);
                    errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' has only two channels; using only the first channel'], 1);
                end
            end
            
            currentImageBaseName = createBaseName(fnmeImages{i});
            currentImageName = strcat(currentImageBaseName, layerSuffix);
            
            % performing segmentation if image is non-logical
            if ~islogical(CurrentImage)
                
                fprintf('\tperforming segmentation for image %s\n', fnmeImages{i});
                
                switch segmentationType
                    case 't'
                        Masks = segmentThreshold(CurrentImage);
                    case 'tw'
                        Masks = segmentWatershed(CurrentImage);
                    otherwise
                        error(['[Engine] Unexpected segmentation type: ', segmentationType]);
                end
                
                SingleCellMetrics{ctr, 1} = currentImageName;
                
                % calculating single-cell shape and intensity metrics
                [SingleCellMetricNames, SingleCellMetrics{ctr, 2}] = calculateSingleCellMetrics(CurrentImage, Masks);
                
            else
                errorReport(end+1) = makeErrorStruct(['file ', fnmeImages{i}, ' is a binary image: did not perform segmentation'], 1);
                Masks = CurrentImage;
            end
            
            % writing mask files to mask folder
            MaskFileName = strcat(MaskInputPath, filesep, currentImageBaseName, layerSuffix, '-MASKS.tif');
            imwrite(Masks, MaskFileName);
            
            ctr = ctr + 1;
        end
    end
    
    % writing single-cell metrics to .csv files
    WriteSingleCellMetrics(OutputPath, SingleCellMetricNames, SingleCellMetrics);
    
else
    MaskInputPath = InputPath;
end

if Mask2GraphRequested
    
    % call Mask2Graph to convert mask info to graph output
    [errorReport, imageProcessingParams] = Mask2Graph(MaskInputPath, OutputPath, thresholdType, threshold, maskLayerNumber, errorReport);
    
else
    
    imageProcessingParams = struct; % empty struct indicates no graph analysis was done
    % transfer masks to output directory
    maskFiles = strcat(MaskInputPath, filesep, '*-MASKS.tif');
    d = dir(maskFiles);
    
    if numel(d) > 0
        movefile(maskFiles, OutputPath);
    end
end
