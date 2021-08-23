function [PRvalues, GTHRvalues, RMSE] = CHROM_DEHAAN(VideoFile, FS, TxtFile)
% CHROM_DEHAAN The Chrominance Method from: De Haan, G., & Jeanne, V. (2013). Robust pulse rate from chrominance-based rPPG. IEEE Transactions on Biomedical Engineering, 60(10), 2878-2886. DOI: 10.1109/TBME.2013.2266196
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

%% Parameters
SkinSegmentTF=false;

LPF = 0.7; %low cutoff frequency (Hz) - specified as 40 bpm (~0.667 Hz) in reference
HPF = 2.5; %high cutoff frequency (Hz) - specified as 240 bpm (~4.0 Hz) in reference

WinSec=1.6; %(was a 32 frame window with 20 fps camera)

%% Add Backup Functions
if(~license('test', 'image_toolbox')&&SkinSegmentTF)
    addpath([cd '\optional\rgb2ycbcr.m']);%GNU GPL rgb2ycbcr.m function
end

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
FramesToRead=ceil(30*VidObj.FrameRate);
gtdata=dlmread(TxtFile);

PRvalues = [];
GTHRvalues = [];

% Create a cascade detector object.
faceDetector = vision.CascadeObjectDetector();

for i = 0:(Duration-30)
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
        
        VidROI = VidFrame;

        
        YCBCR = rgb2ycbcr(VidROI);
        Yth = YCBCR(:,:,1)>80;
        CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
        CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
        ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
        RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
        
%         RGB(FN,:) = sum(sum(VidROI,2)) ./ (size(VidROI,1)*size(VidROI,2));

    end%endwhile video

    if(~license('test', 'image_toolbox')&&SkinSegmentTF)%remove path if added
        rmpath([cd '\optional\']);
    end

    %% CHROM:
    NyquistF = 1/2*FS;
    [B,A] = butter(3,[LPF/NyquistF HPF/NyquistF]);%Butterworth 3rd order filter - originally specified as an a FIR band-pass filter with cutoff frequencies 40-240 BPM

    %Window parameters - overlap, add with 50% overlap
    WinL = ceil(WinSec*FS);
    if(mod(WinL,2))%force even window size for overlap, add of hanning windowed signals
        WinL=WinL+1;
    end
    
    NWin = floor((FN-WinL/2)/(WinL/2));
    S = zeros(NWin,1);
    WinS = 1;%Window Start Index
    WinM = WinS+WinL/2;%Window Middle Index
    WinE = WinS+WinL-1;%Window End Index

    for j = 1:NWin
        TWin = T(WinS:WinE,:);

        RGBBase = mean(RGB(WinS:WinE,:));
        RGBNorm = bsxfun(@times,RGB(WinS:WinE,:),1./RGBBase)-1;

        % CHROM
        Xs = squeeze(3*RGBNorm(:,1)-2*RGBNorm(:,2));%3Rn-2Gn
        Ys = squeeze(1.5*RGBNorm(:,1)+RGBNorm(:,2)-1.5*RGBNorm(:,3));%1.5Rn+Gn-1.5Bn

        if(isfinite(double(Xs)))
            Xf = filtfilt(B,A,double(Xs));
        else
            imshow(VidROI);
        end
        
        if(isfinite(double(Ys)))
            Yf = filtfilt(B,A,double(Ys));
        else
            imshow(VidROI);
        end
        
        Alpha = std(Xf)./std(Yf);

        SWin = Xf - Alpha.*Yf;

        SWin = hann(WinL).*SWin;
        %overlap, add Hanning windowed signals
        if(j==1)
            S = SWin;
            TX = TWin;
        else
            S(WinS:WinM-1) = S(WinS:WinM-1)+SWin(1:WinL/2);%1st half overlap
            S(WinM:WinE) = SWin(WinL/2+1:end);%2nd half
            TX(WinM:WinE) = TWin(WinL/2+1:end);
        end

        WinS = WinM;
        WinM = WinS+WinL/2;
        WinE = WinS+WinL-1;
    end

    BVP=S;
    T=T(1:length(BVP));

    % Estimate Pulse Rate from periodogram
    PR = prpsd(BVP,FS,40,240);
    PRvalues = [PRvalues PR];
    
    txtend = i+30+1;
    gtHR = gtdata(2,i+1:txtend); 
    GTHRvalues = [GTHRvalues mean(gtHR)];

end

RMSE = sqrt(mean((PRvalues-GTHRvalues).^2));


