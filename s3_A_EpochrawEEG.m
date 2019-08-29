%% File one: making the raw EEG file
%based on IG's Raw_EEG but with corrections for catch timing.

cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment')
basefol=pwd;
%%
for ippant = 1:29
    %%
    cd(basefol)
    % find ppant folder
    fileis = dir([pwd filesep num2str(ippant) '_*' ]);
    %change into that directory
    cd(fileis.name);
    
% allcatchEEG = [];
    allEEG = dir([pwd filesep '*.vhdr']);
for ieeg = 1:length(allEEG)
    
    clear trialtrack catchinfo trialstart adjustedcatch EpochedEEG EpochEEG
    loadfile = allEEG(ieeg).name;
    
    EEG= pop_loadbv(pwd, num2str(loadfile));
    
    eegchans = EEG.data(1:64,:);
    PhotoDiode = EEG.data(65,:);
    AudioChan=EEG.data(66,:);
    
    Timing = (1:length(PhotoDiode))/1000;
    %%
    trialtrack = [];
    catchinfo = [];
    for ievent = 1:length(EEG.event)
        eventype = EEG.event(ievent).type;
        if strcmp(eventype, 'S 66'); %catchmarker
            catchinfo =[catchinfo  EEG.event(ievent).latency];
        end
        if strcmp(eventype, 'S 88') %trial end
            trialtrack =  [trialtrack EEG.event(ievent).latency];
        end
    end
    trialstart = trialtrack - (60*EEG.srate); %60 second trials
    adjustedcatch = catchinfo-trialstart;
    
   
    
    
    %% Now epoch
    EpochedEEG=zeros(length(trialstart), 67, 60*EEG.srate); %60 seconds trials.
    
    for itrial = 1:length(trialstart)
        trial= EEG.data(:, trialstart(itrial):trialtrack(itrial));
        EpochEEG(itrial,:,:) = trial;
    end
    switch ieeg
        case 1
            data1=EpochEEG;
            allcatchEEG = [allcatchEEG, adjustedcatch];
        case 2
            data2=EpochEEG;
            allcatchEEG = [allcatchEEG, adjustedcatch];
            
    end
end
    
    %EpochedEEGdata = cat(1, data1, data2);
    EpochedEEGdata = cat(1, data1, data2);
   savename =  ['P' num2str(ippant) 'RawEEG'];
    save(savename, 'EpochedEEGdata')
    
    save(['P' num2str(ippant) 'Trigger_Catch_Timing'], 'allcatchEEG')
end