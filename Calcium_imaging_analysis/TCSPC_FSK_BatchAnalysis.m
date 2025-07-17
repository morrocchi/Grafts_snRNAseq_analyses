%%%This script is to analyze data from fiber photometry loading files from
%%%an excel file containing all the information regarding the video and the
%%%calcium trace.

%%% 01/03/2023 - UPDATED SZ - Added the possibility to select the mode for
%%% photocounting + added the option to select the type of analysis for
%%% response characterization.

%%% 21/07/2023 - UPDATED SZ - Removed the zScore calculation and the dFF as
%%% it introduced artificat due to negative values - Added the area
%%% calculation on the signal corrected for photobleaching only

close all
clear all

%%%First we define the path
dffmode=0; %Set 0 if you want to calculate f0 as average of full trace, otherwise 1 for manual detection
MovingWindow_Correction= 10000; %set the width of the moving window used to adjust for photobleaching
MovingWindow_Raw=20; %set the width of the moving window used to filter raw data

%Define the timewindow for the statistical comparison Pre vs Post-stim *based on
%frequency of the calcium recording!!
StatWindow=40; %Define Window for statistical comparison before and after stimulus
BlankWindow=10; %Option to add a blank window to get further apart from the stimulus onset (both ways)
FrameWindow=60; %Number of frames before and after the event 


%We add the option to decide if we want to evaluate the peak, the sum or
%the mean when counting photons from the photometry system.
CountingMode=0; % - 0=mean ; 1=sum; 2=max 

%We add the option to select the type of analysis to be performed GG= WE
%USED ALL OPTION.
Analysis_Type='ALL';  %Write 'ALL' if you want to perform the analysis for the zscore signal calculated over the whole trace, otherwise write 'STIM'

%For the STIM analysis - the dFF trace is first cut and then the zscore is
%calculated over the selected recording period. The average zscore value in
%the pre stimulus window is then calculated and subtracted to the zscored
%trace
PreStimWindow=3; %Define the pre-stimulus window in SECONDS


%%Then we load the LookupTable
Foldername='...'; %Add the folder name where all calcium recordings are
filename='....xlsx';%Add the folder name where the .xlsx file with the metadata is


DataTable=readtable(fullfile(Foldername,filename));


%%Now we cycle for each recording
for thisrecording=1:size(DataTable,1)

    %%We load the calcium trace

    CaFilename=[DataTable.FileName{thisrecording} '.asc'];
    FileType='fc';
    AcquisitionFreq=DataTable.CaFR(thisrecording);
    Data2=ascimport(Foldername,CaFilename,FileType,AcquisitionFreq,CountingMode);
    Data2(:,2)=smooth(Data2(:,2),MovingWindow_Raw,'moving');

    %%Calculate the smoothed signal to be subtracted to correct for
    %%photobleaching
    Correction=smooth(Data2(:,2),MovingWindow_Correction,'moving');
    Data_Corrected=Data2;
    Data_Corrected(:,2)=Data2(:,2)-Correction; %Subtract the drifting curve from the data


    figure
    subplot(2,1,1)
    plot(Data2(:,1),Data2(:,2),'k');
    hold on
    plot(Data2(:,1),Correction,'r','LineWidth',1);
    subplot(2,1,2)
    plot(Data_Corrected(:,1),Data_Corrected(:,2),'b');


    StimTimestamps=[];

    Events=[];
    b=DataTable.StartMov{thisrecording}; % Data with FSK aligned events
    bb=strsplit(b,',');
    for ii=1:length(bb);tt=bb{ii};tt2=str2double(tt);Events(ii)=tt2;end
    Events=Events';
    Events(isnan(Events(:,1)))=[];
    Events(:,2)=0;
    Events(:,3)=thisrecording;
    Events(:,4)=DataTable.Rat(thisrecording);
    Events(:,5)=DataTable.Anesthesia(thisrecording); %Update-1/2/23 - Add the info for anesthesia for each event

    Events2=[];
    b=DataTable.Stationary{thisrecording}; % Data with RANDOM aligned events
    bb=strsplit(b,',');
    for ii=1:length(bb);tt=bb{ii};tt2=str2double(tt);Events2(ii)=tt2;end
    Events2=Events2';
    Events2(isnan(Events2(:,1)))=[];
    Events2(:,2)=1;
    Events2(:,3)=thisrecording;
    Events2(:,4)=DataTable.Rat(thisrecording);
    Events2(:,5)=DataTable.Anesthesia(thisrecording); %Update-1/2/23 - Add the info for anesthesia for each event

    %%We look at all timestamps events
    StimTimestamps=[Events;Events2];
    DataTable.Events{thisrecording}=StimTimestamps;


    %We calculate the dFF
    DataTable.Fluo{thisrecording}=Data_Corrected(:,2);
    DataTable.TimeFluo{thisrecording}=Data_Corrected(:,1);
    %Calculate the dff
    if dffmode==1
        figure
        plot(Data_Corrected(:,2));
        [x,~]=ginput(2);
        f0=mean(Data_Corrected(x(1):x(2),2));
        close all
    else
        f0=mean(Data_Corrected(:,2));
    end


    dff_zscore=Data_Corrected(:,2); %%% === TO BE ADJUSTED AND TO DECIDE WHICH IS THE FINAL WAY TO MEASURE THE DFF.  
    DataTable.dFFzScore{thisrecording}=dff_zscore;

    %Now we cycle across all events;
    StimResponse=zeros(size(StimTimestamps,1),2*FrameWindow+1);
    for thisevent=1:size(StimTimestamps,1)

        b=(StimTimestamps(thisevent,1)-DataTable.FrameStart(thisrecording))/DataTable.VideoFR(thisrecording);
        time=max(find(Data_Corrected(:,1)<=b));
        
        if time-FrameWindow>0
            %%Analysis mode zscore ALL - It first calculates the zscore
            %%over the full dFF trace and then it cuts for the time window
            %%selected
            if strcmp(Analysis_Type,'ALL')

                StimResponse(thisevent,:)=dff_zscore(time-FrameWindow:time+FrameWindow);% to choose what to analyse dff_zscore 

            else
                %%Analysis mode zscore STIM - It first cuts the responses over
                %%the selected time window from the dFF. It then calculates the
                %%zscore over the cut trace and it finally subtract the average
                %%zscore in the pre-time window
                dFF_window=dff(time-FrameWindow:time+FrameWindow); %first we cut the dFF trace for the selected time window
                %We then calculate the zScore
                zdFF_window=(dFF_window-mean(dFF_window))/std(dFF_window);
                %We calculate the mean in the pre-stim window
                %Transform the pre time window from seconds to frame
                PreStimWindow_FR=PreStimWindow*AcquisitionFreq;
                pre_zdff=mean(zdFF_window(FrameWindow-PreStimWindow_FR:FrameWindow));

                %We calculate the response by subtracting the mean and we
                %store it
                StimResponse(thisevent,:)=zdFF_window-pre_zdff;
            end

        end

    end

    AllStimResponse{thisrecording}=StimResponse;
    AllStimType{thisrecording}=StimTimestamps;

end


%% Reorganize data - Generate a matrix with all recordings and a matrix with all the related information

AllRecordings=cat(1,AllStimResponse{:});
AllInfo=cat(1,AllStimType{:});
timeline=-FrameWindow/DataTable.CaFR:1/DataTable.CaFR:FrameWindow/DataTable.CaFR;

%% Remove AVG Baseline calculated before stimulus (2.5s before Stim)
PreWindow=50; %Define here the number of frames before stimulus
AVGBaseline= mean(AllRecordings(:,FrameWindow-PreWindow:FrameWindow),2);
STDBaseline=std(AllRecordings(:,FrameWindow-PreWindow:FrameWindow),0,2);
AllRecordings_Norm=(AllRecordings-AVGBaseline);


%% Plot and Analyse data - Group data based on animals

% - 1 - Generate a plot where you compare dF/F across all recordings for
% FSK or RANDOM alignment for each animal.

%For each animal we need to find all the events and split them based on
%FSK or RANDOM alignment

%First find the total number of animals - AnimalID
AnimalNumber=unique(AllInfo(:,4));

%We cycle across animals
for thisanimal=1:length(AnimalNumber)
    %Find all FSK Events
    ttx=find(AllInfo(:,2)==0&AllInfo(:,5)==0&AllInfo(:,4)==AnimalNumber(thisanimal)); %change value colum 5, 0=awake, 1=anteshesia
    %Find all RANDOM Events
    ttx2=find(AllInfo(:,2)==1&AllInfo(:,5)==0&AllInfo(:,4)==AnimalNumber(thisanimal));%change value colum 5, 0=awake, 1=anteshesia

    %Generate plot for each animal where you see all the EVOKED events across all sessions
    figure
    subplot(2,2,1) %PLOT FSK EVENTS
    imagesc(timeline,1:1:length(ttx),AllRecordings(ttx,:));
    hold on
    plot([0 0],[1 length(ttx)],'r--','LineWidth',2) 
    xlabel('time to stimulus (s)')
    ylabel('DF/F');
    title('Rat FSK-All events')

    subplot(2,2,2) %PLOTS RANDOM ALIGNED EVENTS 
    imagesc(timeline,1:1:length(ttx2),AllRecordings(ttx2,:));
    hold on
    plot([0 0],[1 length(ttx2)],'r--','LineWidth',2) 
    xlabel('time to stimulus (s)')
    ylabel('DF/F');
    title('Rat  Random-All events')

    subplot(2,2,3)%PLOT FSK AVERAGE PEAKE-SHAPE 
    boundedline(timeline,nanmean(AllRecordings(ttx,:),1),nansem(AllRecordings(ttx,:),1),'k')
    xlim([-FrameWindow/DataTable.CaFR(1) FrameWindow/DataTable.CaFR(1)])
    ylim=([-5 10]);
    tt=ylim;
    hold on
    plot([0 0],[tt(1) tt(2)],'r--','LineWidth',2)
    xlabel('time to stimulus (s)')
    ylabel('DF/F')
 
    subplot(2,2,4)%PLOT RANDOM AVERAGE PEAKE-SHAPE 
    boundedline(timeline,nanmean(AllRecordings(ttx2,:),1),nansem(AllRecordings(ttx2,:),1),'k')
    xlim([-FrameWindow/DataTable.CaFR(1) FrameWindow/DataTable.CaFR(1)])
    ylim=([-5 10]);
    tt=ylim;
    hold on
    plot([0 0],[tt(1) tt(2)],'r--','LineWidth',2)
    xlabel('time to stimulus (s)')
    ylabel('DF/F')
    
end

%% Plot and Analyse data - Group data based on !!trials!!

% - 1 - Generate a HEATMAP plot where you compare dF/F across all recordings for
% stationary and movement for all trials across all animals - AWAKE

%Find all start of FSK events
ttx=find(AllInfo(:,2)==0&AllInfo(:,5)==0);
%Find all Start Random Events
ttx2=find(AllInfo(:,2)==1&AllInfo(:,5)==0);

%Generate HEATMAP plot for all events where you see all the FSK and
%stationary events across all sessions
figure
subplot(2,2,1)
imagesc(timeline,1:1:length(ttx),AllRecordings(ttx,:));
hold on
plot([0 0],[1 length(ttx)],'r--','LineWidth',2) %PLOTS FSK-ALIGNED EVENTS
title('Graft-FSK-All Events' )
xlabel('time to stimulus (s)')
ylabel('dF/F');

subplot(2,2,2) 
imagesc(timeline,1:1:length(ttx2),AllRecordings(ttx2,:));
hold on
plot([0 0],[1 length(ttx2)],'r--','LineWidth',2) %PLOTS RANDOM-ALIGNED EVENTS
title('Graft-Random-All Events' )
xlabel('time (s)')
ylabel('dF/F');

subplot(2,2,3) %PLOTS FSK-ALIGNED EVENTS
boundedline(timeline,nanmean(AllRecordings(ttx,:),1),nansem(AllRecordings(ttx,:),1),'k')
xlim([-FrameWindow/DataTable.CaFR(1) FrameWindow/DataTable.CaFR(1)])
ylim=([-5 10]);
tt=ylim;
hold on
plot([0 0],[tt(1) tt(2)],'r--','LineWidth',2)
xlabel('time to stimulus (s)')
ylabel('dF/F')
title('Average trace')

subplot(2,2,4) %PLOTS RANDOM-ALIGNED EVENTS
boundedline(timeline,nanmean(AllRecordings(ttx2,:),1),nansem(AllRecordings(ttx2,:),1),'k')
xlim([-FrameWindow/DataTable.CaFR(1) FrameWindow/DataTable.CaFR(1)])
ylim=([-5 10]);
tt=ylim;
hold on
plot([0 0],[tt(1) tt(2)],'r--','LineWidth',2)
xlabel('time to stimulus (s)')
ylabel('dF/F')
title('Average trace')

%%Run statistical analysis comparing Pre window vs Post-stim Window 
%FSK AND RANDOM - All variables
PreVar_1=mean(AllRecordings(ttx,(FrameWindow+1-StatWindow-BlankWindow):(FrameWindow-BlankWindow)),2); % Pre-Stim FSK
PreVar_2=mean(AllRecordings(ttx2,(FrameWindow+1-StatWindow-BlankWindow):(FrameWindow-BlankWindow)),2); % Pre-stim RANDOM

PostVar_1=mean(AllRecordings(ttx,(FrameWindow+1+BlankWindow):(FrameWindow+StatWindow+BlankWindow)),2); %Post-Stim FSK
PostVar_2=mean(AllRecordings(ttx2,(FrameWindow+1+BlankWindow):(FrameWindow+StatWindow+BlankWindow)),2); %Post-stim RANDOM

pval(1,:)=signrank(PreVar_1,PostVar_1); %Pre vs Post, condition 1 FSK
pval(2,:)=signrank(PreVar_2,PostVar_2); %Pre vs Post, condition 2 RANDOM


pval2(1,:)=signrank(PreVar_1,PreVar_2); %Pre condition 1 vs 2
pval2(2,:)=signrank(PostVar_1,PostVar_2); %Post condition 1 vs 2

FluoCondition1=AllRecordings(ttx,:);
FluoCondition2=AllRecordings(ttx2,:);

pval_time=NaN(3,size(FluoCondition2,2));
for ii=1:size(FluoCondition1,2)
    pval_time(1,ii)=signrank(FluoCondition1(:,ii),FluoCondition2(:,ii)); %Condition1 vs Condition 2

end