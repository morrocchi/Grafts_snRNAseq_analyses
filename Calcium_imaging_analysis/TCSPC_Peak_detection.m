%%% Script for the analysis of Calcium Signal obtained from Fiber
%%% Photometry technique.

%%We load the LookupTable

% Here Add the path where to find all the recordings of the dataset and the
% relative xlsx file with the related MetaData to be used for the analysis.
%The .xlsx file should be structured such as the following variables has
%been defined (also see the example in the shared folder):
% - Rat ID
% - File Name (of calcium recordings)
% - Calcium Frame Rate
% - Camera Frame Rate
% - Camera starting frame
% - Onset frame for first shock start
% - Offset frame for first shock end
% - Onset frame for second shock start
% - Offset fram for second shock end

Foldername=['...'];
filename=['...'];

DataTable=readtable(fullfile(Foldername,filename));

%Allow to decide which type of baseline subtraction to be used
BaselineSubType=2; %Put 0 if you want to use Low/High pass filters; put 1 if you want to use the smoothing and subtraction mode, otherwise put 2 to use analysis from Bruno et al., 2021
STDThr=2.5; %Define the number of standard deviation you want to consider for peak detection

%Add option to plot raw data or the filtered one with the peaks
%detected.

PlotRaw=0; %Put 1 if you want to plot raw data, otherwise set 0.

if BaselineSubType==0
    %Set the filtering parameters
    LowPassFilter=5; %Set the value of the low pass filter
    HighPassFilter=0.02; %Set the value of the high pass filter
elseif BaselineSubType==1

    MovingWindow_Correction= 10000; %set the width of the moving window used to adjust for photobleaching
    MovingWindow_Raw=20; %set the width of the moving window used to filter raw data
else
    MovingWindowFil=20; % width of the moving window used for low-pass filter data (at 20 Hz)
    CutoffFiltr=5; % low-pass cutoff frequency (Hz, if CutoffFiltr=0 Filter is disabled)
end


%We add the option to decide if we want to evaluate the peak, the sum or
%the mean when counting photons from the photometry system.
CountingMode=0; % - 0=mean ; 1=sum; 2=max -
BaselineTimeWindow=400; %Time in seconds to evaluate the baseline frequency
StimulusTimeWindow=2; %Time in seconds to evaluate the post-stimulus frequency
PeakWindow=2.5; %Define the time window for the transient shape evaluation in seconds
PeakDistance=0.7; %Define the minimum peak distance in seconds
Prominence=0.5; %Define peak prominence

%Initialize Variables
PeaksFrequency_baseline=zeros(size(DataTable,1),1);
PeaksFrequency_stimulus=zeros(size(DataTable,1),2);
PeakNumber_baseline=zeros(size(DataTable,1),1);
PeaksNumber_stimulus=zeros(size(DataTable,1),2);

%%Now we cycle for each recording
for thisrecording=1:size(DataTable,1)

    %%We load the calcium trace

    CaFilename=[DataTable.FileName{thisrecording} '.asc'];
    FileType='fc';
    AcquisitionFreq=DataTable.CaFR(thisrecording);
    Data2=ascimport(Foldername,CaFilename,FileType,AcquisitionFreq,CountingMode);
    
    if BaselineSubType==0
        %Filter the raw data
        Data_filtered=highpass(Data2(:,2),HighPassFilter,AcquisitionFreq);
        Data_filtered=lowpass(Data_filtered,LowPassFilter,AcquisitionFreq);
    elseif BaselineSubType==1
        %Smooth and subtract the raw trace
        Data2(:,2)=smooth(Data2(:,2),MovingWindow_Raw,'moving');
        Correction=smooth(Data2(:,2),MovingWindow_Correction,'moving');
        Data_Corrected=Data2;
        Data_Corrected(:,2)=Data2(:,2)-Correction; %Subtract the drifting curve from the data
        Data_filtered=Data_Corrected(:,2);

    else
        %We apply the analysis with a low pass filter only and an overall smoothing
        Data_filtered=Filtr(Data2(:,2),MovingWindowFil,CutoffFiltr,AcquisitionFreq);
        Correction=smooth(Data_filtered,length(Data_filtered),'moving'); % correction for photobleaching (span=entire trace)
        %Store the Data_filtered variable before baseline subtraction
        Data_filtered_unsub=Data_filtered;

        Data_filtered=Data_filtered-Correction+nanmean(Correction);

    end

    DataAVG=mean(Data_filtered);
    DataSTD=std(Data_filtered);
    Thr=DataAVG+STDThr*DataSTD;

    [pks,PeaksFrames] = findpeaks(Data_filtered,'MinPeakHeight',Thr,'MinPeakProminence',Prominence,'MinPeakDistance',PeakDistance*DataTable.CaFR(thisrecording));
    

    %We want to evaluate the Spikes Frequencies at baseline in the first
    %400s segment of recordings
    
    %%% TO FIND THE STIMULUS ONSET WE HAVE TO CONSIDER THAT THE FREQUENCY FOR THE
    %%% VIDEO IS 40 WHILE THE ONE FOR THE Ca2+ IS 20.
    %Set stimuli onsets
    
    Shock1_frame=(DataTable.Shock1_st(thisrecording)-DataTable.FrameStart(thisrecording))/DataTable.VideoFR(thisrecording); %We get the time of the first shock / or evaluated segment for spontaneuos activity
    Shock2_frame=(DataTable.Shock2_st(thisrecording)-DataTable.FrameStart(thisrecording))/DataTable.VideoFR(thisrecording); %We get the time of the second shock
    
    Shock1_frame=Shock1_frame*DataTable.CaFR(thisrecording); %We transform time in frame based on Calcium Acquisition Frequency
    Shock2_frame=Shock2_frame*DataTable.CaFR(thisrecording); %We transform time in frame based on Calcium Acquisition Frequency

    %Evaluate the baseline Frequency of peaks
    tt=find(PeaksFrames>=(Shock1_frame-BaselineTimeWindow*DataTable.CaFR(thisrecording))&PeaksFrames<=Shock1_frame);
    PeaksFrequency_baseline(thisrecording)=length(tt)/(BaselineTimeWindow);
    PeakNumber_baseline(thisrecording)=length(tt);
    
    %Now find the frequency of peaks after the first stimulus / applied  %just to FSK analysis
    tt2=find(PeaksFrames>Shock1_frame&PeaksFrames<=(Shock1_frame+StimulusTimeWindow*DataTable.CaFR(thisrecording)));
    PeaksFrequency_stimulus(thisrecording,1)=length(tt2)/(StimulusTimeWindow);
    PeaksNumber_stimulus(thisrecording,1)=length(tt2);

    %Now find the frequency of peaks after the second stimulus  %just to FSK analysis
    tt3=find(PeaksFrames>Shock2_frame&PeaksFrames<=(Shock2_frame+StimulusTimeWindow*DataTable.CaFR(thisrecording)));
    PeaksFrequency_stimulus(thisrecording,2)=length(tt3)/(StimulusTimeWindow);
    PeaksNumber_stimulus(thisrecording,2)=length(tt3);
    
    %Create a time variable
    TimeVariable=1/DataTable.CaFR(thisrecording):1/DataTable.CaFR(thisrecording):length(Data_filtered)/DataTable.CaFR(thisrecording);

    figure
    if PlotRaw==1
        plot(TimeVariable,Data2(:,2),'b','LineWidth',1);
    else
        plot(TimeVariable,Data_filtered,'b','LineWidth',1);
        hold on
        plot(PeaksFrames/DataTable.CaFR(thisrecording),pks,'k*');% black= no classified peaks
        plot(PeaksFrames(tt)/DataTable.CaFR(thisrecording),pks(tt),'m*');% in blue= baseline peaks
        plot(PeaksFrames(tt2)/DataTable.CaFR(thisrecording),pks(tt2),'r*'); %red= peaks after first stimulus
        plot(PeaksFrames(tt3)/DataTable.CaFR(thisrecording),pks(tt3),'g*'); %green= peaks after second stimulus
        yLim=ylim;
        yLim=[200,-10];
        plot([Shock1_frame/DataTable.CaFR(thisrecording) Shock1_frame/DataTable.CaFR(thisrecording)],[yLim(1) yLim(2)],'r--','LineWidth',1)
        plot([Shock2_frame/DataTable.CaFR(thisrecording) Shock2_frame/DataTable.CaFR(thisrecording)],[yLim(1) yLim(2)],'r--','LineWidth',1)
        title(CaFilename)
    end

    xlabel('Time (s)')
    ylabel('dF/F)')

    %Collect average Peak Shape
    %We have to consider sometimes we might have peaks at the real
    %beginning or at the real end, so the time window considered will
    %exceed the limits of the recordings. We decide to remove those peaks
    %from the analysis
    PeakShape=NaN(length(PeaksFrames),2*PeakWindow*DataTable.CaFR(thisrecording)+1);PeakClassifier=NaN(1,length(PeaksFrames));
    for ii=1:length(PeaksFrames)
        if PeaksFrames(ii)+PeakWindow*DataTable.CaFR(thisrecording)>=length(Data_filtered)
        elseif PeaksFrames(ii)-PeakWindow*DataTable.CaFR(thisrecording)<=0
        else
            PeakShape(ii,:)=Data_filtered((PeaksFrames(ii)-PeakWindow*DataTable.CaFR(thisrecording)):(PeaksFrames(ii)+PeakWindow*DataTable.CaFR(thisrecording)));

            if ismember(ii,tt2)||ismember(ii,tt3)
                PeakClassifier(ii)=1;
            else
                PeakClassifier(ii)=0;
            end
        end
    end
    
    %We now have NaNs where the peak was close to the beginning or the end
    %of the recording so we can remove it.

    PeakShape(isnan(PeakClassifier),:)=[];
    PeakClassifier(isnan(PeakClassifier))=[];

    All_PeakShape{thisrecording}=PeakShape;
    All_PeakClassifier{thisrecording}=PeakClassifier;
    All_PeaksFrames{thisrecording}=PeaksFrames;

end

%% Analyse peak shape

%%% Cycle across recordings and substract the overall baseline - calculate the average
%%% across all peaks for each recording and subtract the first value.

%Get the number of recordings
RecordingNumber=size(All_PeaksFrames,2);

%Cycle across recordings and normalize each peak for the baseline
%calculated across all peaks for baseline and stimulus condition
ALlPeaks_Stim=[]; AllPeaks_Baseline=[];AllPeaks=[];
for thisrec=1:RecordingNumber
    peaks=All_PeakShape{thisrec};
    classifier=All_PeakClassifier{thisrec};

    stim_peaks=peaks(classifier==1,:);
    basal_peaks=peaks(classifier==0,:);

    %Normalize Peaks for their AVG baseline 
    %removing the first value of the AVG shape
    bsl=nanmean(stim_peaks,1);
    stim_peaks=stim_peaks-bsl(1);

    bsl=nanmean(basal_peaks,1);
    basal_peaks=basal_peaks-bsl(1);
    
    %Now we subistitute in the variable containing all peaks the traces
    %normalized for their basal
    peaks(classifier==1,:)=stim_peaks;
    peaks(classifier==0,:)=basal_peaks;

    if thisrec==1
        ALlPeaks_Stim=stim_peaks;
        AllPeaks_Baseline=basal_peaks;
        AllPeaks=peaks;
    else
        ALlPeaks_Stim=[ALlPeaks_Stim;stim_peaks];
        AllPeaks_Baseline=[AllPeaks_Baseline;basal_peaks];
        AllPeaks=[AllPeaks;peaks];
    end
end

% %First concatenate all detected peaks
AllPeakClassifier=cat(2,All_PeakClassifier{:})';
% 
ALlPeaks_Stim=AllPeaks(AllPeakClassifier==1,:);
AllPeaks_Baseline=AllPeaks(AllPeakClassifier==0,:);
% 
% %Subtract baseline according to Claudio Script "Ca_Ana_Peak.m"-Based on Bruno et al,. 2021  - Calculate
% %the mean across all peaks and subtract the value of the first time point
bsl=nanmean(AllPeaks_Baseline);
AllPeaks_Baseline=AllPeaks_Baseline-bsl(1);
% 
bsl=nanmean(ALlPeaks_Stim);
ALlPeaks_Stim=ALlPeaks_Stim-bsl(1);

%Generate the time axis for plotting the clusters
timeShape=-PeakWindow:1/DataTable.CaFR:PeakWindow;

%Calculate the 95% CI for baseline peaks- Spontaneous activity
AVG=mean(AllPeaks_Baseline,1);
STD=std(AllPeaks_Baseline,1);

N=size(AllPeaks_Baseline,1);
ySEM = STD/sqrt(N);  
ts = tinv([0.025  0.975],N-1);
Baseline_CI95 = bsxfun(@times, ySEM, ts(:));

%Calculate the 95% CI for stimuli peaks - exclusive for FSK
AVG=mean(ALlPeaks_Stim,1);
STD=std(ALlPeaks_Stim,1);
N=size(ALlPeaks_Stim,1);
ySEM = STD/sqrt(N);  
ts = tinv([0.025  0.975],N-1);
STIM_CI95 = bsxfun(@times, ySEM, ts(:));

%Now plot the shape for each cluster for the peak detected in spontaneous activty 
figure
subplot(2,2,1) % spontaneous activty
plot(timeShape,AllPeaks_Baseline,'k');
xlim([-PeakWindow PeakWindow])
title('All Spontaneous Peaks')
subplot(2,2,2)
boundedline(timeShape,mean(AllPeaks_Baseline,1),Baseline_CI95(1,:)','k')
xlim([-PeakWindow PeakWindow])
ylim([-0.5 2])
title('AVG Spontaneous Shape- 95CI')

%% Function for raw data post-processing based on Bruno et al. 2021 "pMAT: An Open-Source Software Suite for the Analysis of Fiber Photometry Data"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xLP = Filtr(x,dimFilt,cutoff,sampling_freq)

numberData=length(x);

Wn=cutoff/(sampling_freq/2);

data=x';
delayFilt=dimFilt/2;

b=fir1(dimFilt,Wn);

a=zeros(1,length(b));

a(1)=1;
meanInitialX=mean(data(1:delayFilt));
meanFinalX=mean(data(numberData-delayFilt:numberData));

data=[data, [meanFinalX*ones(1,delayFilt)]]';
dataFiltr=filter(b,a,data);
dataFiltr=[dataFiltr((delayFilt+1:length(data)),:)];
xLP=dataFiltr;
xLP(1:delayFilt)=meanInitialX;
xLP(length(xLP)-delayFilt:length(xLP))=meanFinalX;

end