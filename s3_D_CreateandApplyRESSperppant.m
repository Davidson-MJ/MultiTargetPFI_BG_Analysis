% Take 3 . Motivation is to create 16 spatial filters per participant
% 4 x locations , 4 x Hz.

%we need all the trials in which the targets were present,




clear all
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/MatthewDavidson/Desktop/FromIrenes Mac/Data Experiment copy')
    addpath('/Users/matthewdavidson/Documents/MATLAB/SSVEP-PFI_BackgroundTag')
end
%%
basefol=pwd;
clearvars -except basefol allppants
dbstop if error


job.concatLiketrialsperppant=1; %also performs RESS, we are left with 16 types of spatial filter, per ppant. 

% job.applyRESSpertrialtype=1; %after constructing RESS filters in the previous step, here we apply to all relevant trials, 
% %also reduces the size of the 


getelocs;

cd(basefol)
cd('newplots-MD')
%% load data to determine physical catch timing.

%remaining participants after behavioral data analysis/exclusion
allppants=1:29;
badppants = [8, 15,28, 4 , 7,5 6];
allppants(badppants)=[];

% window=[-2 4];
window=[-3 3];

srate=250;

epochdur = sum(abs(window))*srate;

timeid = [0:1/srate:epochdur];
timeid= timeid-3;

onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing

%%
if job.concatLiketrialsperppant==1
    
    
    
    % SET UP RESS params.
%     peakfreqsare=[8,13,15,18,20,40];
    peakfreqsare=[16,26,30,36];
    
    neighbour=[];
    for ihz = 1:length(peakfreqsare)
    
        hzAT = peakfreqsare(ihz);
        
        %needs to be narrow to avoid other stim freqs.  
    neighbour(ihz).centrefreqs=[hzAT-1 , hzAT+1];   % a low and high centre freq for covR
    neighbour(ihz).sd=[.5, .5]; %FWHM of those...
    neighbour(ihz).snrlow = [hzAT-1.5 hzAT-.5];    %not so important, were being used in SNR calcs. now in a laterscript
    neighbour(ihz).snrhigh = [hzAT+.5 hzAT+1.5];
    
    
    end
%     
%     %13 new
%     neighbour(2).centrefreqs=[11 15];   % a low and high centre freq for covR %
%     neighbour(2).sd=[1, 1];
%     neighbour(2).snrlow = [10.5 11.5];
%     neighbour(2).snrhigh = [14.5 15.5];
%     
%     %15 Hz
%     neighbour(3).centrefreqs=[14 16];   % a low and high centre freq for covR %
%     neighbour(3).sd=[1, 1];
%     neighbour(3).snrlow = [13.5 14.5];
%     neighbour(3).snrhigh = [15.5 16.5];
%     
%     %18 Hz
%     neighbour(4).centrefreqs=[17 22];   % a low and high centre freq for covR %
%     neighbour(4).sd=[1, 1];
%     neighbour(4).snrlow = [16.5 17.5];
%     neighbour(4).snrhigh = [18.5 22.5];
%     
%     
%     %  %% 20
%     neighbour(5).centrefreqs=[18 , 22]; %+1  % a low and high centre freq for covR
%     neighbour(5).sd=[1, 1];
%     neighbour(5).snrlow = [17.5 18.5]; %+1
%     neighbour(5).snrhigh = [21.5 22.5]; %+1
%     
% 
%     
%     %
%     %         %% 40
%     neighbour(6).centrefreqs=[38 , 42];   % a low and high centre freq for covR
%     neighbour(6).sd=[1, 1];
%     neighbour(6).snrlow = [37.5 38.5];
%     neighbour(6).snrhigh = [41.5 42.5];
%     %
%     
    
    %%
    windowsmall=[];
    windowsmall(1,:) = [-3 -1] ;%  window targets present
    windowsmall(2,:) = [1 3] ;%  window targets present after BP
    %%
%     usechans=[51:64, 17:32]; %occipital
    usechans=[1:64]; %occipital
    
    
    for ifol =10%allppants
        
        
        %%
        cd(basefol)
        cd(num2str(ifol))
        
        
        %% load the relevant PFI data.
        load('ppant_PFI_Epoched');
        load('ppant_Catch_Epoched')
        load('TrialIndicesbyLocationandHz.mat')
        
        %%
        
        ressEVEC_byHzxLoc = zeros(length(peakfreqsare),4,length(usechans));
        ressEVEC_byHzxLoc_MAPS=zeros(length(peakfreqsare),4,length(usechans));
        
        for ifreq=1:length(peakfreqsare)
            
            for iloc=1:4
                
                
                
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
              
                %%
                for id=1:10% Use all epochs so as not to bias condition comparisons.
                    
                    switch id
                        case 1
                            dataIN=ppant_SNREEG_PFI_0_1;
                            searchtrials = [Freqwas.dir0_1.Trialind];
%                             durscheck=durs0_1;
                            usewindow=1:2; 
                        case 2
                            dataIN=ppant_SNREEG_PFI_1_0;
                            searchtrials = [Freqwas.dir1_0.Trialind];
%                             durscheck=durs1_0;
                            usewindow=1:2; 
                        case 3
                            dataIN=ppant_SNREEG_PFI_1_2;
                            searchtrials = [Freqwas.dir1_2.Trialind];
%                             durscheck=durs1_2;
                            usewindow=1:2; 
                        case 4
                            dataIN=ppant_SNREEG_PFI_2_1;
                            searchtrials = [Freqwas.dir2_1.Trialind];
%                             durscheck=durs2_1;
                            usewindow=1:2; 
                        case 5
                            dataIN=ppant_SNREEG_PFI_3_2;
                            searchtrials = [Freqwas.dir3_2.Trialind];
%                             durscheck=durs3_2;
                            usewindow=1:2;
                        case 6
                            dataIN=ppant_SNREEG_PFI_2_3;                            
                            searchtrials = [Freqwas.dir2_3.Trialind];
%                             durscheck=durs2_3;
                            usewindow=1:2;
                        case 7
                            dataIN=ppant_SNREEG_disapBPwithincatch;
                            searchtrials = [1:24]; 
                            usewindow=1; 
                        case 8
                            dataIN=ppant_SNREEG_reapBPaftercatch;
                            searchtrials = [1:24]; %which catch was disappearing or not?
                            usewindow=2; 
                        case 9
                            dataIN=ppant_SNREEG_catchramponset;
                            searchtrials = [1:24]; %which catch was disappearing or not?
                            usewindow=1;
                        case 10
                            dataIN=ppant_SNREEG_catchrampoffset;
                            searchtrials = [1:24]; %which catch was disappearing or not?
                            usewindow=2;
                        
                        
                            
                    end
                    
                    % which trials can we look into?
                    
                    if ifreq>4 %|| id>6 % we can use all trials if we had 20, 40hz of interest (or catch)
                    
                        relevanttrials=1:24;
                    else
                        switch iloc
                            case 1
                                relevanttrials= TopLeftTrialindexbyHz(ifreq,:);
                            case 2
                                relevanttrials= TopRightTrialindexbyHz(ifreq,:);
                            case 3
                                relevanttrials= BottomLeftTrialindexbyHz(ifreq,:);
                            case 4
                                relevanttrials= BottomRightTrialindexbyHz(ifreq,:);
                        end
                    end
                
                    
                    
                    
                    % extractthe equivalent types:
                    trialtypesperTargPresent = find(ismember(searchtrials, relevanttrials));
                    % only check durs for PFI.
%                     if id<6
%                     durscheck = durscheck(trialtypesperTargPresent);
%                     else
                      durscheck=[];
%                     end
                        
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    for windch=1:length(usewindow)
                        wind=usewindow(windch);
                    windcheck= windowsmall(wind,:);
                    
                    tidx=dsearchn(timeid', [windcheck]');
                    
                    
                    %reduce size.
                    datast= squeeze(dataIN(trialtypesperTargPresent,:,tidx(1):tidx(2)));
                    
                    
                    %% check for bad trials (noisy)
                    %std per trial(average over electrodes)
                    tmp=[];
                    if ndims(datast)<3
                        tmp(1,:,:)= datast;
                        datast=tmp;
                    end
                    datastSD = nanstd(squeeze(nanmean(datast,2)),0,2);
                    
                    %remove those with 2.5*std from mean.
                    trialSD=nanstd(datastSD);
                    mSD=nanmean(datastSD);
                    keeptrials=1:size(datast,1);
                    
                    %remove the trials with excessive noise.
                    badtrials=find(datastSD>mSD+2.5*trialSD)';
                    
                    % also skip the trials which were only transient button
                    % presses. (less than one second).                    
                    shorttrials = find(durscheck<60);
                    badtrials = [badtrials, shorttrials];
                    
                    % remove these from consideration.
                    keeptrials(badtrials)=[];
                    datast=datast(keeptrials,:,:);
                    
                    
                    % concatenate all these types to increase quality of cov matrix.
                    if id==1 && windch==1%start here
                        if ndims(datast)<3
                        dataOUT(1,:,:)=datast;
                        else
                            dataOUT=datast;
                        end
                        
                    else
                        % append
                         if ndims(datast)<3
                        tmp(1,:,:)=datast;
                        dataOUT=cat(1,dataOUT,tmp);
                        
                        else
                            dataOUT=cat(1,dataOUT,datast);
                         end                        
                        
                    end
                    
                    end
                    
                end
                
                % We want a single component for this type of config, per ppant.
                % trials as last dimension
                %%
                ressd=permute(dataOUT,[2 3 1]);
                peakfreq1= peakfreqsare(ifreq);
                
                
                
                [ressEVEC, ~, ~, ressMAPS, hz]= constrRESSacrosstrialsGROUPEDCOMBO(ressd, peakfreq1, usechans, neighbour, window,srate);
                
                %
                
                % now we can multiple each type of data by the appropriate
                %                 % ressEVEC
                %%                
%                 figure(2)
%                                 hold on
%                                 subplot(211)
%                                 topoplot(ressMAPS, elocs(usechans));
% %                                 subplot(212)
%                                 plot(hz, ressSNR);
%                                 xlim([0 45])
%                                 hold on
%                                 plot([peakfreq1 peakfreq1], ylim, ['r:'])
%                                 shg
% %                                 
                %%
                ressEVEC_byHzxLoc(ifreq,iloc,:) = ressEVEC;
                ressEVEC_byHzxLoc_MAPS(ifreq,iloc,:)= ressMAPS;
            end
            
            
        end
        if length(peakfreqsare)~=6
            savename ='RESSfilterdata_2fTG';
        else
            savename='RESSfilterdata';
        end
        save(savename, 'ressEVEC_byHzxLoc', 'ressEVEC_byHzxLoc_MAPS', 'usechans', 'neighbour', 'windowsmall')
        disp(['fin for ppant' num2str(ifol)])
    end
end



