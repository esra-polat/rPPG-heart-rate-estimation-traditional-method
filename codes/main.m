function [RMSE_green, RMSE_pos, RMSE_ica, RMSE_chrom] = main(VideoFile, FS, TxtFile)
    addpath(genpath('tools'))
    
    %% Add Backup Functions
    if(~license('test', 'image_toolbox'))
        addpath([cd '\optional\rgb2ycbcr.m']);
    end
    %% Parameters
    LPF = 0.7; 
    HPF = 2.5; 
    WinSec=1.6; 
    %% Load Video:
    VidObj = VideoReader(VideoFile);
    Duration = floor(VidObj.Duration);
    FramesToRead=ceil(30*VidObj.FrameRate); 
    gtdata=dlmread(TxtFile);
    
    PRvalues_green = [];    
    PRvalues_ica = [];
    PRvalues_chrom = [];
    PRvalues_pos = [];
    GTHRvalues = [];
    
    faceDetector = vision.CascadeObjectDetector();
    for i = 0:(Duration-30)
        VidObj.CurrentTime = i;
        FN = 0;
        %% Read Video and Spatially Average:
        T = zeros(FramesToRead,1);
        RGB = zeros(FramesToRead,3);
        fprintf("i = %d  ", i);
        while hasFrame(VidObj) && (VidObj.CurrentTime <= i+30)
            FN = FN+1;
            T(FN) = VidObj.CurrentTime;
            VidFrame = readFrame(VidObj);
            %% face detection
            bbox = step(faceDetector, VidFrame);
            if ~isempty(bbox)
                VidFrame = insertShape(VidFrame, 'Rectangle', bbox);
                VidFrame = imcrop(VidFrame, bbox(1, :));
            end
            %fprintf("i = %d and FN = %d\n", i ,FN);
            VidROI = VidFrame;
            %% skin segmentation 
            YCBCR = rgb2ycbcr(VidROI);
            Yth = YCBCR(:,:,1)>80;
            CBth = (YCBCR(:,:,2)>77).*(YCBCR(:,:,2)<127);
            CRth = (YCBCR(:,:,3)>133).*(YCBCR(:,:,3)<173);
            ROISkin = VidROI.*repmat(uint8(Yth.*CBth.*CRth),[1,1,3]);
            RGB(FN,:) = squeeze(sum(sum(ROISkin,1),2)./sum(sum(logical(ROISkin),1),2));
            % RGB(FN,:) = sum(sum(VidROI));
            RGB_green_ica(FN,:) = squeeze(sum(sum(ROISkin,1), 2));
        end
        
        %% GREEN
        
        BVP_green = RGB_green_ica(:,2);
        NyquistF_green = 1/2*FS;
        [B_green,A_green] = butter(3,[LPF/NyquistF_green HPF/NyquistF_green]);
        if(isfinite(double(BVP_green)-mean(BVP_green)))
            BVP_F_green = filtfilt(B_green,A_green,(double(BVP_green)-mean(BVP_green)));
        end
        BVP_green = BVP_F_green;
        PR_green = prpsd(BVP_green,FS,40,240);
        PRvalues_green = [PRvalues_green PR_green];
        
        %% CHROM
 
        if(~license('test', 'image_toolbox'))
            rmpath([cd '\optional\']);
        end
        
        NyquistF_chrom = 1/2*FS;
        [B_chrom,A_chrom] = butter(3,[LPF/NyquistF_chrom HPF/NyquistF_chrom]);
        WinL_chrom = ceil(WinSec*FS);
        if(mod(WinL_chrom,2))
            WinL_chrom=WinL_chrom+1;
        end
        NWin_chrom = floor((FN-WinL_chrom/2)/(WinL_chrom/2));
        S_chrom = zeros(NWin_chrom,1);
        WinS_chrom = 1;
        WinM_chrom = WinS_chrom + WinL_chrom / 2;
        WinE_chrom = WinS_chrom + WinL_chrom - 1;
        T_chrom = T;
        for j = 1:NWin_chrom
            TWin_chrom = T_chrom(WinS_chrom:WinE_chrom,:);
            RGBBase_chrom = mean(RGB(WinS_chrom:WinE_chrom,:));
            RGBNorm_chrom = bsxfun(@times,RGB(WinS_chrom:WinE_chrom,:),1./RGBBase_chrom)-1;
            Xs_chrom = squeeze(3*RGBNorm_chrom(:,1)-2*RGBNorm_chrom(:,2));
            Ys_chrom = squeeze(1.5*RGBNorm_chrom(:,1)+RGBNorm_chrom(:,2)-1.5*RGBNorm_chrom(:,3));
            if(isfinite(double(Xs_chrom)))
                Xf_chrom = filtfilt(B_chrom,A_chrom,double(Xs_chrom));
            end
            if(isfinite(double(Ys_chrom)))
                Yf_chrom = filtfilt(B_chrom,A_chrom,double(Ys_chrom));
            end
            Alpha_chrom = std(Xf_chrom)./std(Yf_chrom);
            SWin_chrom = Xf_chrom - Alpha_chrom.*Yf_chrom;
            SWin_chrom = hann(WinL_chrom).*SWin_chrom;
            if(j==1)
                S_chrom = SWin_chrom;
                TX_chrom = TWin_chrom;
            else
                S_chrom(WinS_chrom:WinM_chrom-1) = S_chrom(WinS_chrom:WinM_chrom-1)+SWin_chrom(1:WinL_chrom/2);
                S_chrom(WinM_chrom:WinE_chrom) = SWin_chrom(WinL_chrom/2+1:end);
                TX_chrom(WinM_chrom:WinE_chrom) = TWin_chrom(WinL_chrom/2+1:end);
            end
            WinS_chrom = WinM_chrom;
            WinM_chrom = WinS_chrom+WinL_chrom/2;
            WinE_chrom = WinS_chrom+WinL_chrom-1;
        end
        BVP_chrom = S_chrom;
        T_chrom=T_chrom(1:length(BVP_chrom));
        PR_chrom = prpsd(BVP_chrom,FS,40,240);
        PRvalues_chrom = [PRvalues_chrom PR_chrom];
        
        %% POS
        
        useFGTransform_pos=true;
        if useFGTransform_pos
            RGBBase_pos = mean(RGB);
            RGBNorm_pos = bsxfun(@times,RGB,1./RGBBase_pos)-1;
            FF_pos = fft(RGBNorm_pos);
            F_pos = (0:size(RGBNorm_pos,1)-1)*FS/size(RGBNorm_pos,1);
            H_pos = FF_pos*[-1/sqrt(6);2/sqrt(6);-1/sqrt(6)];
            W_pos = (H_pos.*conj(H_pos))./sum(FF_pos.*conj(FF_pos),2);
            FMask_pos = (F_pos >= LPF)&(F_pos <= HPF);
            % FMask_pos(length(FMask_pos)/2+1:end)=FMask_pos(length(FMask_pos)/2:-1:1);
            FMask_pos = FMask_pos + fliplr(FMask_pos);
            W_pos=W_pos.*FMask_pos';
            FF_pos = FF_pos.*repmat(W_pos,[1,3]);
            RGBNorm_pos=real(ifft(FF_pos));
            RGBNorm_pos = bsxfun(@times,RGBNorm_pos+1,RGBBase_pos);
            RGB_pos=RGBNorm_pos;
        end
        N_pos = size(RGB,1);
        H_pos = zeros(1,N_pos);
        l_pos = ceil(WinSec*FS);
        C_pos = zeros(length(l_pos),3);
        for n = 1:N_pos-1
            m = n - l_pos + 1;
            if(m > 0)
                Cn_pos = ( RGB(m:n,:) ./ mean(RGB(m:n,:)) )';
                S_pos = [0, 1, -1; -2, 1, 1] * Cn_pos;
                h_pos = S_pos(1,:) + ((std(S_pos(1,:)) / std(S_pos(2,:))) * S_pos(2,:));
                H_pos(m:n) = H_pos(m:n) + (h_pos - mean(h_pos));
            end
        end
        BVP_pos=H_pos;
        T_pos = T;
        T_pos=T_pos(1:length(BVP_pos));
        
        PR_pos = prpsd(BVP_pos,FS,40,240);
       
        PRvalues_pos = [PRvalues_pos PR_pos];
        
        %% ICA
        NyquistF_ica = 1/2*FS;
        RGBNorm_ica=zeros(size(RGB_green_ica));
        Lambda_ica=100;
        for c=1:3
            RGBDetrend_ica= spdetrend(RGB_green_ica(:,c),Lambda_ica); 
            RGBNorm_ica(:,c) = (RGBDetrend_ica - mean(RGBDetrend_ica))/std(RGBDetrend_ica); 
        end
        if(isfinite(double(RGBNorm_ica)))
            [W_ica,S_ica] = ica(RGBNorm_ica',3); 
        end
        MaxPx_ica=zeros(1,3);
        for c=1:3
            FF_ica = fft(S_ica(c,:));
            F_ica=(1:length(FF_ica))/length(FF_ica)*FS*60;
            FF_ica(1)=[];
            N_ica=length(FF_ica);
            Px_ica = abs(FF_ica(1:floor(N_ica/2))).^2;
            Fx_ica = (1:N_ica/2)/(N_ica/2)*NyquistF_ica;
            Px_ica=Px_ica/sum(Px_ica);
            MaxPx_ica(c)=max(Px_ica);
        end
        [M_ica,MaxComp_ica]=max(MaxPx_ica(:));
        BVP_I_ica = S_ica(MaxComp_ica,:);
        [B_ica,A_ica] = butter(3,[LPF/NyquistF_ica HPF/NyquistF_ica]);
        if(isfinite(double(BVP_I_ica)))
            BVP_F_ica = filtfilt(B_ica,A_ica,double(BVP_I_ica));
        end
        BVP_ica=BVP_F_ica;
        PR_ica = prpsd(BVP_ica,FS,40,240);
        PRvalues_ica = [PRvalues_ica PR_ica];
        
        
        %% ------------------------------------
        txtend = i+30+1;
        gtHR = gtdata(2,i+1:txtend); 
        GTHRvalues = [GTHRvalues mean(gtHR)];
        
    end
RMSE_green = sqrt(mean((PRvalues_green-GTHRvalues).^2));
RMSE_chrom = sqrt(mean((PRvalues_chrom-GTHRvalues).^2));
RMSE_pos = sqrt(mean((PRvalues_pos-GTHRvalues).^2));
RMSE_ica = sqrt(mean((PRvalues_ica-GTHRvalues).^2));
