clc, clear;
%% Set variables:
DataDirectory	= [cd '/test_data/'];
VideoFile       = [DataDirectory '43.avi'];
FS              = 30;
TxtFile         = [DataDirectory '43.txt'];
%% Check Files
if(~exist(VideoFile,'file'))
    error('VideoFile: ''%s'' does not exist.',VideoFile1)
end
[RMSE_green, RMSE_pos, RMSE_ica, RMSE_chrom] = main(VideoFile, FS, TxtFile);
display(RMSE_green);
display(RMSE_pos);
display(RMSE_ica);
display(RMSE_chrom);
%%===========================================================================================
% %% Green - Verkruysse, Svaasand, Nelson (2008)
% [PRvalues, GTHRvalues, RMSE] = GREEN_VERKRUYSSE(VideoFile, FS, TxtFile);
% fprintf('GREEN_VERKRUYSSE:\n')
% display(PRvalues);
% display(GTHRvalues);
% display(RMSE)
% 
% %% ICA - Poh, McDuff, Picard (2010)
% [PRvalues, GTHRvalues, RMSE] = ICA_POH(VideoFile, FS, TxtFile);
% fprintf('ICA_POH:\n')
% display(PRvalues);
% display(GTHRvalues);
% display(RMSE)
% 
% %% CHROM - De Haan & Jeanne (2013)
% [PRvalues, GTHRvalues, RMSE] = CHROM_DEHAAN(VideoFile, FS, TxtFile);
% fprintf('CHROM_DEHAAN:\n')
% display(PRvalues);
% display(GTHRvalues);
% display(RMSE)
% 
% %% POS - Wang, den Brinker, Stuijk & de Haan (2017)
% [PRvalues, GTHRvalues, RMSE] = POS_WANG(VideoFile, FS, TxtFile);
% fprintf('POS_WANG:\n')
% display(PRvalues);
% display(GTHRvalues);
% display(RMSE)
