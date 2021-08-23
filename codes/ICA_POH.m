function [PRvalues, GTHRvalues, RMSE] = ICA_POH(VideoFile, FS, TxtFile)
% ICA_POH The Independent Component Analysis (ICA) Method from: Poh, M. Z., McDuff, D. J., & Picard, R. W. (2010). Non-contact, automated cardiac pulse measurements using video imaging and blind source separation. Optics express, 18(10), 10762-10774. DOI: 10.1364/OE.18.010762
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
LPF = 0.7; %low cutoff frequency (Hz) - 0.7 Hz in reference
HPF = 2.5; %high cutoff frequency (Hz) - 4.0 Hz in reference

% %% Plot Control
% if(PlotTF)
%     PlotPRPSD = true;
%     PlotSNR = true;
% else
%     PlotPRPSD = false;
%     PlotSNR = false;
% end

%% Load Video:
VidObj = VideoReader(VideoFile);
Duration = floor(VidObj.Duration);
FramesToRead=ceil(30*VidObj.FrameRate); %video may be encoded at slightly different frame rate
gtdata=dlmread(TxtFile);

PRvalues = [];
GTHRvalues = [];

% Create a cascade detector object.
faceDetector = vision.CascadeObjectDetector();
BVP_F = [0 0];

for i = 0:(Duration-30)
    VidObj.CurrentTime = i;
    FN = 0;
    
    %% Read Video and Spatially Average:
    T = zeros(FramesToRead,1);%initialize time vector
    RGB=zeros(FramesToRead,3);%initialize color signal
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
        
        VidROI = VidFrame;

    
        % skin segmentation 
        YCBCR = rgb2ycbcr(VidROI);
        Yth = YCBCR(:,:,1)>80;
        CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
        CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
        ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
        RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
    
        %RGB(FN,:) = sum(sum(VidROI));%if different size regions are used for different frames, the signals should be normalized by the region size, but not necessary for whole frame processing or constant region size
    end

    %% Detrend & ICA:
    NyquistF = 1/2*FS;
    RGBNorm=zeros(size(RGB));
    Lambda=100;
    for c=1:3
        RGBDetrend= spdetrend(RGB(:,c),Lambda); %M. P. Tarvainen, TBME, 2002
        RGBNorm(:,c) = (RGBDetrend - mean(RGBDetrend))/std(RGBDetrend); %normalize to zero mean and unit variance
    end
    [W,S] = ica(RGBNorm',3); %JADE ICA - J. F. Cardoso 1997, G. D. Clifford, MIT, 2004

    %% Select BVP Source:
    % Component with maximum normalized (by total power) power
    MaxPx=zeros(1,3);
    for c=1:3
        FF = fft(S(c,:));
        F=(1:length(FF))/length(FF)*FS*60;
        FF(1)=[];
        N=length(FF);
        Px = abs(FF(1:floor(N/2))).^2;
        Fx = (1:N/2)/(N/2)*NyquistF;
        Px=Px/sum(Px);
        MaxPx(c)=max(Px);
    end
    [M,MaxComp]=max(MaxPx(:));
    BVP_I = S(MaxComp,:);

    %% Filter, Normalize
    %originally specified in reference with 5-point moving average, bandpass
    %filter, and cubic-spine interpolation to 256Hz
    [B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter
    
    if(isfinite(double(BVP_I)))
        BVP_F = filtfilt(B,A,double(BVP_I));
    else
        imshow(VidROI);
    end

    BVP=BVP_F;

    % Estimate Pulse Rate from periodogram
    PR = prpsd(BVP,FS,40,240);
    PRvalues = [PRvalues PR];

    txtend = i+30+1;
    gtHR = gtdata(2,i+1:txtend); 
    GTHRvalues = [GTHRvalues mean(gtHR)];
end

RMSE = sqrt(mean((PRvalues-GTHRvalues).^2));


