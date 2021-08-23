function [PRvalues, GTHRvalues, RMSE] = GREEN_VERKRUYSSE(VideoFile, FS, TxtFile)
% GREEN_VERKRUYSSE The Green-Channel Method from: Verkruysse, W., Svaasand, L. O., & Nelson, J. S. (2008). Remote plethysmographic imaging using ambient light. Optics express, 16(26), 21434-21445. DOI: 10.1364/OE.16.021434
% 
%   Inputs:
%       VideoFile               = Video file path.
%       FS                      = Video framerate (fps).
%       StartTime               = Timepoint at which to start process (default = 0 seconds).
%       Duration                = Duration of the time window to process (default = 60 seconds).
%       ECGFile                 = File path to corresponding ECG data file (.mat) containing: 1) The waveform - ECGData.data, 2) The ECG sampling rate - ECGData.fs, 3) The ECG peak locations (in samples) - ECGData.peaks.
%       PPGFile                 = File path to corresponding PPG data file (.mat) containing: 1) The waveform - PPGData.data, 2) The PPG sampling rate - PPGData.fs, 3) The PPG peak locations (in samples) - PPGData.peaks.
%       PlotTF                  = Boolean to turn plotting results on or off.
%
%   Outputs:
%       BVP                     = Processed Blood Volume Pulse (BVP).
%       PR                      = Estimated Pulse Rate (PR) from processed BVP timeseries using peak in periodogram.
%       HR_ECG                  = Gold standard Heart Rate (HR) measured from the ECG timeseries R-waves for the window.
%       PR_PPG                  = Pulse Rate (PR) measured from the PPG timeseries systolic onsets for the window.
%       SNR                     = Blood Volume Pulse Signal-to-Noise Ratio (SNR) calculated based on the ECG HR frequency using a method adapted from the method by G. de Haan, TBME, 2013.
%
%   Requires - Signal Processing Toolbox
%
% Daniel McDuff, Ethan Blackford, January 2019
% Copyright (c)
% Licensed under the MIT License and the RAIL AI License.

addpath(genpath('tools'))

%% Parameters
LPF = 0.7; %low cutoff frequency (Hz) - 0.8 Hz in reference
HPF = 2.5; %high cutoff frequency (Hz) - both 6.0 Hz and 2.0 Hz used in reference

% %% Plot Control
% if(PlotTF)
%     PlotPRPSD = true;
%     PlotSNR = true;
% else
%     PlotPRPSD = false;
%     PlotSNR = false;
% end
% % 
%% Load Video:
VidObj = VideoReader(VideoFile);
Duration = floor(VidObj.Duration);
FramesToRead=ceil(30*VidObj.FrameRate); %video may be encoded at slightly different frame rate
gtdata=dlmread(TxtFile);

PRvalues = [];
GTHRvalues = [];

% Create a cascade detector object.
faceDetector = vision.CascadeObjectDetector();

for i = 1:(Duration-30)
    VidObj.CurrentTime = i;
    FN = 0;
    %% Read Video and Spatially Average:
    T = zeros(FramesToRead,1);%initialize time vector
    RGB = zeros(FramesToRead,3);%initialize color signal
    
    while hasFrame(VidObj) && (VidObj.CurrentTime <= i+30)
        FN = FN+1;
        T(FN) = VidObj.CurrentTime;
        VidFrame = readFrame(VidObj);
    
        %% face detection

        bbox         = step(faceDetector, VidFrame);
       
        if ~isempty(bbox)
             % Draw the returned bounding box around the detected face.
            VidFrame = insertShape(VidFrame, 'Rectangle', bbox);

            VidFrame = imcrop(VidFrame, bbox(1, :));
        end
        
        %fprintf("i = %d and FN = %d\n", i ,FN);
        
        
        %position for optional face detection/tracking - originally specified in reference as a manual segmentation.
        VidROI = VidFrame;
    
        % skin segmentation 
        YCBCR = rgb2ycbcr(VidROI);
        Yth = YCBCR(:,:,1)>80;
        CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
        CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
        ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
        %RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
        RGB(FN,:) = squeeze(sum(sum(ROISkin,1), 2));%if different size regions are used for different frames, the signals should be normalized by the region size, but not necessary for whole frame processing or constant region size
    end
    
    %% Select BVP Source:
    % Green channel
    BVP = RGB(:,2);

    %% Filter, Normalize
    NyquistF = 1/2*FS;
    [B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter - originally specified in reference with a 4th order butterworth using filtfilt function
    
    display(double(BVP)-mean(BVP));
 
    BVP_F = filtfilt(B,A,(double(BVP)-mean(BVP)));
   

    BVP = BVP_F;

    % Estimate Pulse Rate from periodogram
    PR = prpsd(BVP,FS,40,240);

    PRvalues = [PRvalues PR];
    
    txtend = i+30+1;
    gtHR = gtdata(2,i+1:txtend); 
    GTHRvalues = [GTHRvalues mean(gtHR)];
    
end

RMSE = sqrt(mean((PRvalues-GTHRvalues).^2));






