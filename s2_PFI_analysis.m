% Calculate new participant PFI data, outside of catch events.
clear all; %close all; clc;
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy/newplots-MD')
catch 
    cd('/Users/mdavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy/newplots-MD')
%     addpath('/Users/matthewdavidson/Documents/GitHub/SSVEP-PFI_BackgroundTag-master');
end
basedir=pwd;



dbstop if error
%%
%load data structure


% goodppants = [1:3, 7:9, 11:14, 16:27,29];
allppants=1:29;
% badppants = [8, 15,28, 4 , 7,5 6,10];
badppants = [8, 15,28, 4 , 7,5 6]; % note that 8,15,28 failed catch time-courseBP, 4:7 had no PFI.
allppants(badppants)=[];

% save('MD_AllBP_perppant.mat', 'goodppants', '-append');
load('MD_AllBP_perppant.mat');
%load Catch Performance structure
load('Catchperformance.mat')
%work with retained participants.


%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% 
excludeTransientlength=0; %frames! 60 = 1 second.
%%%% %%%% %%%% %%%% %%%% %%%% %%%% %%%% 

fontsize=25;

job.removeCatchperiodsfromBPtrace=0;

%Three types of analysis (IVs( RMfactors)). Number absent, location absent, and flicker
%speed.
job.calcPFIdataperNum=1; %resaves into PFI_only with new fields in structure..
job.calcPFIdataperFreqandLoc=0; %resaves into PFI_only with new fields in structure..


%PFI analysis (nPFI)
job.concatPFIacrossPpants_num=1;
job.concatPFIacrossPpants_freqandLoc=0;

% shuffling analysis. ( n=1000 shuffles).
job.createShufflePFIdata_pernum=1;
job.calcPFIdataperNum_shuffled=1;
job.concatPFIacrossPpants_num_shuffled=1;

job.calcPFIdataperNum=1; %resaves into PFI_only with new fields in structure..
job.calcPFIdataperFreqandLoc=0; %resaves into PFI_only with new fields in structure..


%archived:
job.concatPFIacrossPpants_slope_nulldistribution=0; %not in use.
% export ppant mean to excel file for analysis in JASP.
job.exportDatatoExcelfolonDesktop=0;


%get ready for LME.
job.exportDatatoExcelfolonDesktop_jamovi =0; %can save as matlab or on desktop.

%perform LME in matlab
job.LMEonHzxLOC=0;
job.LMEnPFIwShuff=0;



%%%%%% plotting
%plot the results for PFI, in separate bar graphs.
job.plotBehaviouraldata=0;

%
job.plotBehaviouraldata_num_with_shuffled=1;

job.compareSlopesofShuffledvsObserveddata=0;



% supplementary analyses.
job.plotPFIdata_asfunctionofTrialNumber=0; %to see if increase in PFI with time on task, by type of PFI.
job.concatPpantDurationsbyPFI_num=1; %for plotting histogram of disap durations, by number gone in PFI

job.plothistofDisapdurationbynumbergone=1; %see above.
job.plothistofDisapdurationbynumbergone_perppants=1; %cut off duration?


%not using anymore:
%redo the stats.
job.RMANOVA_PFI_byNum=0;
job.RMANOVA_PFI_byNum_vs_shuffled=0;


job.plotlikelihoodofPFIovertime=0;

%%
if job.removeCatchperiodsfromBPtrace==1
    %%
    PFI_only=[];
    for ippant=1:29 %do all ppants, just analyze / plot retained, in case we change criteria at some point.
        
        BPtmp= ppantTrialDatawithDetails(ippant).AllBPData;
        catchremovedallBPtmp=BPtmp;
        for itrial=1:24
            %some trials had catch start / end ,yet no actual target removal.
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TotalCatchTargetsRemoved > 0
                %new analysis, replacing the catch target with NaNs, to remove from averaging.
                
                
                catchstart= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchstart_frames;
                catchend=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).Catchend_frames;
                
                
                
                for iloc=1:4
                    switch iloc
                        case 1 %check TL
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Catch;
                        case 2
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Catch;
                        case 3
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Catch;
                        case 4
                            checkloc=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Catch;
                    end
                    
                    if checkloc==1 %if this BP location contained catch event
                       %convert section of that BP trace to nan
                       nantrain = nan(1,(catchend-catchstart)+1);
                       tmptrace=squeeze(BPtmp(itrial,iloc,:))';
                       tmptrace(1,catchstart:catchend)=nantrain;
                       catchremovedallBPtmp(itrial,iloc, :)=tmptrace;
                    end
                end
                
                catchdur= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).totalCatchdur;	
                
                PFI_only(ippant).Trial(itrial).allBPs= squeeze(catchremovedallBPtmp(itrial,:,:));
                PFI_only(ippant).Trial(itrial).resultTriallength= length(BPtmp) - catchdur;
                
            else
                
                PFI_only(ippant).Trial(itrial).allBPs= squeeze(BPtmp(itrial,:,:));
                PFI_only(ippant).Trial(itrial).resultTriallength= length(BPtmp);
            end
        end
    end
    save('PFI_data', 'PFI_only')
end
%%

%%
if job.calcPFIdataperNum==1
    %%
    load('PFI_data.mat')
    %%
    
%     excludeTransientlength=30; %frames minimum for real PFI.
    for ippant = 1:29
        
        
        for itrial=1:24
            allBPtmp= PFI_only(ippant).Trial(itrial).allBPs;
            accumBP = squeeze(nansum(allBPtmp,1)); %combine each locations BP.
            

            
            
            
            
            if excludeTransientlength>0;%; %ms
                % then adjust for short transient button presses:
                
                figure(1); clf
                plot(accumBP, 'k', 'linew', 2); %plot original.
                hold on
                % find all diffs
                sws = find(diff(accumBP));
                
                %how long between BPs?
                diff_sws = diff(sws);
                
                rmv = find(diff_sws<excludeTransientlength);
                % now work through and replace.
                for irm= rmv
                    
                    tremovesw =sws(irm);
                    %replace values for this period:
                    vals = accumBP(1,tremovesw);
                    %replace to = sws
                    trep = sws(irm+1);
                    %replace now:
                    accumBP(1,tremovesw:trep)=vals;
                    
                    
                end
                
%                 hold on;
                plot(accumBP, 'r','linew', 3); %check method. 
                legend('original', 'new'); ylabel('accumBP');
                ylim([-.1 4])
  PFI_only(ippant).Trial(itrial).newaccumBPs= accumBP;
            end
            
            
            
                        
            
            ZeroTargetDisap=0;
            OneTargetDisap=0;
            TwoTargetDisap=0;
            ThreeormoreTargetDisap=0;
            
            framestamp0_begin=[];
            framestamp0_end=[];
            framestamp1_begin=[];
            framestamp1_end=[];
            framestamp2_begin=[];
            framestamp2_end=[];
            framestamp3_begin=[];
            framestamp3_end=[];
            
            
            
            rowis = (itrial+1) + 24*(ippant-1);
            
            %6th column is logical index for trial reject
            if strcmp(catchStruct{rowis,6},'0') && sum(accumBP)>0
                
                for timeDur=2:(length(accumBP)-1)
                    
                    if accumBP(1,timeDur)==1 %if single button is pressed.
                        
                        
                        if accumBP(1,timeDur-1)~=1 %single button onset
                            if framestamp1_begin>0
                                framestamp1_begin= [framestamp1_begin, timeDur];
                            else
                                framestamp1_begin= timeDur;
                            end
                            OneTargetDisap=OneTargetDisap+1;
                            
                        end
                        
                        
                        if accumBP(1,timeDur+1)~=1 || timeDur==length(accumBP)-1
                            %singlebutton release
                            
                            if framestamp1_end>0
                                framestamp1_end= [framestamp1_end, timeDur+1];
                            else
                                framestamp1_end= timeDur+1;
                            end
                        end
                        
                        
                    elseif accumBP(1,timeDur)==2
                        
                        
                        
                        if accumBP(1,timeDur-1)~=2 %double button onset
                            if framestamp2_begin>0
                                framestamp2_begin= [framestamp2_begin, timeDur];
                            else
                                framestamp2_begin= timeDur;
                            end
                            TwoTargetDisap=TwoTargetDisap+1;
                            
                        end
                        
                        if accumBP(1,timeDur+1)~=2 || timeDur==length(accumBP)-1
                            %double button offset
                            
                            if framestamp2_end>0
                                framestamp2_end= [framestamp2_end, timeDur+1];
                            else
                                framestamp2_end= timeDur+1;
                            end
                        end
                        
                        
                    elseif  accumBP(1,timeDur)==3
                        
                        
                        
                        if accumBP(1,timeDur-1)~=3 %double button onset
                            if framestamp3_begin>0
                                framestamp3_begin= [framestamp3_begin, timeDur];
                            else
                                framestamp3_begin= timeDur;
                            end
                            
                            ThreeormoreTargetDisap=ThreeormoreTargetDisap+1;
                            
                        end
                        
                        if accumBP(1,timeDur+1)~=3 || timeDur==length(accumBP)-1
                            %triple button offset
                            
                            if framestamp3_end>0
                                framestamp3_end= [framestamp3_end, timeDur+1];
                            else
                                framestamp3_end= timeDur+1;
                            end
                        end
                        
                        
                    
                 elseif accumBP(1,timeDur)==0 % now included for new figs.
                    if accumBP(1,timeDur-1) ~=0 %start of a disappearance.
                        if framestamp0_begin>0
                            framestamp0_begin = [framestamp0_begin, timeDur];
                        else
                            framestamp0_begin = timeDur;
                        end
                        
                        ZeroTargetDisap=ZeroTargetDisap+1;
                        
                        
                    end
                    
                    if accumBP(1,timeDur+1)~=0 || timeDur==length(accumBP)-1 %accounting for end of trial.
                        %singlebutton offset
                        
                        if framestamp0_end>0
                            framestamp0_end= [framestamp0_end, timeDur+1];
                        else
                            framestamp0_end= timeDur+1;
                        end
                        
                    end
                       
                end
                end
                %in case a button was held until the end of a trial, we
                %need to remove these last disapperances
                
                
                %adjust for trial beginning.
                if accumBP(1,1)==0
                    framestamp0_begin = [1 framestamp0_begin];                    
                    
%                     ZeroTargetDisap = length(framestamp0_begin);
                end
                
%                 if length(framestamp0_begin) > length(framestamp0_end)
%                     %adjust length
%                     framestamp0_begin= framestamp0_begin(1:length(framestamp0_end));
%                     
%                 end
                
                
                
                if length(framestamp1_begin) > length(framestamp1_end)
                    %adjust length
                    framestamp1_begin= framestamp1_begin(1:length(framestamp1_end));
                    OneTargetDisap = length(framestamp1_begin);
                end
                %same for 2 and three target instances.
                if  length(framestamp2_begin) > length(framestamp2_end)
                    
                    framestamp2_begin= framestamp2_begin(1:length(framestamp2_end));
                    TwoTargetDisap = length(framestamp2_begin);
                    
                end
                
                if  length(framestamp3_begin) > length(framestamp3_end)
                    
                    framestamp3_begin= framestamp3_begin(1:length(framestamp3_end));
                    ThreeormoreTargetDisap = length(framestamp3_begin);
                end
                
                alldurs0 = framestamp0_end - framestamp0_begin;
                alldurs1 = framestamp1_end - framestamp1_begin;
                alldurs2 = framestamp2_end - framestamp2_begin;
                alldurs3= framestamp3_end - framestamp3_begin;
                
                dur0Disap = sum(alldurs0);
                dur1Disap = sum(alldurs1);
                dur2Disap = sum(alldurs2);
                dur3Disap = sum(alldurs3);
                
                totalrec= sum([dur0Disap, dur1Disap, dur2Disap, dur3Disap]);
                
                
                
                   %for each type, remove the 'disappearances' which were only
                %transients.
                for in=1:4
                    switch in
                        case 1
                            ddurs=alldurs0;
                            nPFI=ZeroTargetDisap;
                            frame_end=framestamp0_end;
                            frame_begin=framestamp0_begin;
                        case 2
                            ddurs=alldurs1;
                            nPFI=OneTargetDisap;
                            frame_end=framestamp1_end;
                            frame_begin=framestamp1_begin;
                        case 3
                            ddurs=alldurs2;
                            nPFI=TwoTargetDisap;
                            frame_end=framestamp2_end;
                            frame_begin=framestamp2_begin;
                        case 4
                            ddurs=alldurs3;
                            nPFI=ThreeormoreTargetDisap;
                            frame_end=framestamp3_end;
                            frame_begin=framestamp3_begin;
                        
                    end
                    
                    keep= find(ddurs>=excludeTransientlength);
%                   
                    
                    %now adjust outgoings to only keep the PFI that were
                    %not transients.
                    durDisap = sum(ddurs(keep));
                    
                    nPFI=length(keep);
                    alldurs = ddurs(keep);
                    
                    frame_begs = frame_begin(keep);
                    frame_ends = frame_end(keep);
                    
                    
                    if nPFI==0 && in==1
                        shg
                    end
                
                    
                    %store
                    switch in
                        case 1
                        
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_0target = durDisap/nPFI/60;
%                       
                            %  total duration per trial (in seconds)
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_0target = durDisap/60;
                            %frequency of PFI per trial (how many times PFI occurs)
                            PFI_only(ippant).Trial(itrial).PFI_num0target_unadjusted = nPFI;
                             %record time stamps (in frames) of these events.
                             PFI_only(ippant).Trial(itrial).PFI_disap_0target_framestart= frame_begs;
                             PFI_only(ippant).Trial(itrial).PFI_disap_0target_frameend= frame_ends;
                             PFI_only(ippant).Trial(itrial).PFI_disap_0target_durs= alldurs;
                             
                             % sanitycheck
%                              figure(1);clf; plot(accumBP); ylim([-1 3]);
%                              total0=durDisap/60;
%                              m0=durDisap/nPFI/60;
%                              title({['total 0 is ' num2str(total0) 'sec'];['crossprod =' num2str(m0*nPFI) '(' num2str(nPFI) '*' num2str(m0) ')']})
%                             shg
%                             shg
                        case 2
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_1target = durDisap/nPFI/60;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_1target = durDisap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num1target_unadjusted = nPFI;
                                PFI_only(ippant).Trial(itrial).PFI_disap_1target_framestart= frame_begs;
                             PFI_only(ippant).Trial(itrial).PFI_disap_1target_frameend= frame_ends;
                             PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs= alldurs;
                        case 3
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_2target = durDisap/nPFI/60;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_2target = durDisap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num2target_unadjusted = nPFI;
                                PFI_only(ippant).Trial(itrial).PFI_disap_2target_framestart= frame_begs;
                             PFI_only(ippant).Trial(itrial).PFI_disap_2target_frameend= frame_ends;
                             PFI_only(ippant).Trial(itrial).PFI_disap_2target_durs= alldurs;
                        case 4
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_3ormoretarget = durDisap/nPFI/60;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_3ormoretarget = durDisap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num3ormoretarget_unadjusted = nPFI;
                                PFI_only(ippant).Trial(itrial).PFI_disap_3ormoretarget_framestart= frame_begs;
                             PFI_only(ippant).Trial(itrial).PFI_disap_3ormoretarget_frameend= frame_ends;
                             PFI_only(ippant).Trial(itrial).PFI_disap_3ormoretarget_durs= alldurs;                        
                        
                end
                end
           
                
                
                %collect
                PFI_only(ippant).Trial(itrial).Goodtrial=1;
            else
                PFI_only(ippant).Trial(itrial).Goodtrial=0;
            end
            
            
        end
        
    end
    
    save('PFI_data', 'PFI_only')
end
%%

if job.calcPFIdataperFreqandLoc==1
    %%
    load('PFI_data.mat')
    %%
%         excludeTransientlength=30; 
    for ippant = 1:29
        
        
        for itrial=1:24
            allBPtmp= PFI_only(ippant).Trial(itrial).allBPs;
            
            
            rowis = (itrial+1) + 24*(ippant-1);
            
            %6th column is logical index for trial reject
            if strcmp(catchStruct{rowis,6},'0') && PFI_only(ippant).Trial(itrial).Goodtrial==1;
                
                for iloc=1:4
                    switch iloc
                        case 1
                            freq=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Freq;
                        case 2
                            freq=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Freq;
                        case 3
                            freq=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Freq;
                        case 4
                            freq=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Freq;
                    end
                    
                    %channel for certain location.
                    BPuse = allBPtmp(iloc,:);
                    
                    %collect stats.
                    
                    %only dealing with single button disaps.
                    dur1Disap=0; %total time PFI
                    OneTargetDisap=0;
                    
                    framestamp1_begin=[];
                    framestamp1_end=[];
                    
                    
                    for timeDur=2:(length(BPuse)-1)
                        
                        if BPuse(1,timeDur)==1
                            
                            
                            if BPuse(1,timeDur-1)~=1 %single button onset
                                if framestamp1_begin>0
                                    framestamp1_begin= [framestamp1_begin, timeDur];
                                else
                                    framestamp1_begin= timeDur;
                                end
                                OneTargetDisap=OneTargetDisap+1;
                                
                            end
                            
                            
                            if BPuse(1,timeDur+1)~=1 || timeDur==length(BPuse)-1
                                %singlebutton offset
                                
                                if framestamp1_end>0
                                    framestamp1_end= [framestamp1_end, timeDur+1];
                                else
                                    framestamp1_end= timeDur+1;
                                end
                            end
                            
                        end
                        
                    end
                    
%                     %in case a button was held until the end of a trial, we
%                     %need to remove these last disapperances
                    
                    if length(framestamp1_begin) > length(framestamp1_end)
                        %adjust length
                        framestamp1_begin= framestamp1_begin(1:length(framestamp1_end));
                        OneTargetDisap = length(framestamp1_begin);
                    end
                    
                    
                    alldurs1 = framestamp1_end - framestamp1_begin;
                    
                    
                 
                    %adjust for very short length windows.
                    ddurs=alldurs1;
                    
                    keep= find(ddurs>=excludeTransientlength);
%                     keep= find(ddurs);
                    
                    %now adjust outgoings to only keep the PFI that were
                    %not transients.
                    durDisap = sum(ddurs(keep));
                    nPFI=length(keep);
                    alldurs = ddurs(keep);
                    
                    framestamp1_begin = framestamp1_begin(keep);
                    framestamp1_end = framestamp1_end(keep);
                    
                    %rename for easy storage:
                    dur1Disap = sum(alldurs);
                    OneTargetDisap=length(framestamp1_begin);
                    alldurs1=alldurs;
                    
                    switch iloc
                        case 1
                            
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_1target_TL = dur1Disap/60/OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_1target_TL = dur1Disap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num1target_unadjusted_TL = OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs_TL= alldurs1;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_TL_framestart= framestamp1_begin;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_TL_frameend= framestamp1_end;
                            
                        case 2
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_1target_TR = dur1Disap/60/OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_1target_TR = dur1Disap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num1target_unadjusted_TR = OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs_TR= alldurs1;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_TR_framestart= framestamp1_begin;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_TR_frameend= framestamp1_end;
                        case 3
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_1target_BL = dur1Disap/60/OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_1target_BL = dur1Disap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num1target_unadjusted_BL = OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs_BL= alldurs1;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_BL_framestart= framestamp1_begin;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_BL_frameend= framestamp1_end;
                        case 4
                            PFI_only(ippant).Trial(itrial).meanPFIduration_perdisap_1target_BR = dur1Disap/60/OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFItotalduration_disap_1target_BR = dur1Disap/60;
                            PFI_only(ippant).Trial(itrial).PFI_num1target_unadjusted_BR = OneTargetDisap;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs_BR= alldurs1;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_BR_framestart= framestamp1_begin;
                            PFI_only(ippant).Trial(itrial).PFI_disap_1target_BR_frameend= framestamp1_end;
                    end
                    %
                    %collect
                    PFI_only(ippant).Trial(itrial).Goodtrial=1;
                    
                end
                PFI_only(ippant).Trial(itrial).TL_Freq= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TL_Freq;
                PFI_only(ippant).Trial(itrial).TR_Freq= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TR_Freq;
                PFI_only(ippant).Trial(itrial).BL_Freq= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BL_Freq;
                PFI_only(ippant).Trial(itrial).BR_Freq= ppantTrialDatawithDetails(ippant).TrialDetails(itrial).BR_Freq;
                
                
               
            else
                PFI_only(ippant).Trial(itrial).Goodtrial=0;
            end
            
        end
    end
    save('PFI_data', 'PFI_only')
end

%%

%%
if job.concatPFIacrossPpants_num==1
    load('PFI_data')
    Freq_NumPFI_acrossTrials = nan(29,24,4);
    mDurperNumPFI_acrossTrials=nan(29,24,4);
    totalDurperNumPFI_acrossTrials=nan(29,24,4);
    %now including all durs also.
    alldurspernPFI_acrossTrials=[];
    
    for ippant = 1:29
        
        for itrial=1:24
            
            if PFI_only(ippant).Trial(itrial).Goodtrial==1 %not rejected
                usedata=PFI_only(ippant).Trial(itrial);
                %zero targets first
                Freq_NumPFI_acrossTrials(ippant,itrial,1) = usedata.PFI_num0target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,1)=usedata.meanPFIduration_perdisap_0target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,1)=usedata.PFItotalduration_disap_0target;              
                %include all durations for histograms
                try alldurspernPFI_acrossTrials(ippant).zeroTarg= [alldurspernPFI_acrossTrials(ippant).zeroTarg, usedata.PFI_disap_0target_durs];
                catch
                    alldurspernPFI_acrossTrials(ippant).zeroTarg= [usedata.PFI_disap_0target_durs];
                end
                
                
                %single
                Freq_NumPFI_acrossTrials(ippant,itrial,2) = usedata.PFI_num1target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,2)=usedata.meanPFIduration_perdisap_1target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,2)=usedata.PFItotalduration_disap_1target;
                %include all durations for histograms
                try alldurspernPFI_acrossTrials(ippant).oneTarg= [alldurspernPFI_acrossTrials(ippant).oneTarg, usedata.PFI_disap_1target_durs];
                catch
                    alldurspernPFI_acrossTrials(ippant).oneTarg= [usedata.PFI_disap_1target_durs];
                end
                
                %double
                
                Freq_NumPFI_acrossTrials(ippant,itrial,3) = usedata.PFI_num2target_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,3)=usedata.meanPFIduration_perdisap_2target;
                totalDurperNumPFI_acrossTrials(ippant,itrial,3)=usedata.PFItotalduration_disap_2target;
                try alldurspernPFI_acrossTrials(ippant).twoTarg= [alldurspernPFI_acrossTrials(ippant).twoTarg, usedata.PFI_disap_2target_durs];
                catch
                    alldurspernPFI_acrossTrials(ippant).twoTarg= [usedata.PFI_disap_2target_durs];
                end
                
                %'2 or more'               
                Freq_NumPFI_acrossTrials(ippant,itrial,4) = usedata.PFI_num3ormoretarget_unadjusted;
                mDurperNumPFI_acrossTrials(ippant,itrial,4)=usedata.meanPFIduration_perdisap_3ormoretarget;
                totalDurperNumPFI_acrossTrials(ippant,itrial,4)=usedata.PFItotalduration_disap_3ormoretarget;
                try alldurspernPFI_acrossTrials(ippant).threeTarg= [alldurspernPFI_acrossTrials(ippant).threeTarg, usedata.PFI_disap_3ormoretarget_durs];
                catch
                    alldurspernPFI_acrossTrials(ippant).threeTarg= [usedata.PFI_disap_3ormoretarget_durs];
                end
                
            end
            
        end
        
        %% sanity check for debugging
%         disp('ppant trials')
%         [Freq_NumPFI_acrossTrials(ippant,:,1)', mDurperNumPFI_acrossTrials(ippant,:,1)', totalDurperNumPFI_acrossTrials(ippant,:,1)']
%         disp('means are ')
%         [nanmean(Freq_NumPFI_acrossTrials(ippant,:,1),2)', nanmean(mDurperNumPFI_acrossTrials(ippant,:,1),2)', nanmean(totalDurperNumPFI_acrossTrials(ippant,:,1),2)']
%         disp('crossprod')
%         
%         [nanmean(Freq_NumPFI_acrossTrials(ippant,:,1),2)'.* nanmean(mDurperNumPFI_acrossTrials(ippant,:,1),2)']
%         
    end
    
 
    %%
    
    
    save('PFI_data_concat', 'Freq_NumPFI_acrossTrials', 'mDurperNumPFI_acrossTrials', 'totalDurperNumPFI_acrossTrials', 'alldurspernPFI_acrossTrials')
    
    
end
%%

if job.concatPFIacrossPpants_freqandLoc==1
    load('PFI_data', 'PFI_only')
    load('PFI_data_concat')
    
    %sort by flicker frequency of target
    Freq_FreqPFI_acrossTrials = nan(29,24,4);
    mDurperFreqPFI_acrossTrials=nan(29,24,4);
    totalDurperFreqPFI_acrossTrials=nan(29,24,4);
    
    %and location
    Freq_LocPFI_acrossTrials = nan(29,24,4);
    mDurperLocPFI_acrossTrials=nan(29,24,4);
    totalDurperLocPFI_acrossTrials=nan(29,24,4);
    
    
    
    
    %and COMBINED!
    Freq_HzxLocPFI_acrossTrials = nan(29,24,4,4);
    mDurper_HzxLocPFI_acrossTrials=nan(29,24,4,4);
    totalDur_HzxLocPFI_acrossTrials=nan(29,24,4,4);
    
    
    for ippant = 1:29
        
        for itrial=1:24
            
            if PFI_only(ippant).Trial(itrial).Goodtrial==1 %not rejected
                usedata=PFI_only(ippant).Trial(itrial);
                
                
                %location is easy, store first,
                Freq_LocPFI_acrossTrials(ippant,itrial,1) = usedata.PFI_num1target_unadjusted_TL;
                Freq_LocPFI_acrossTrials(ippant,itrial,2) = usedata.PFI_num1target_unadjusted_TR;
                Freq_LocPFI_acrossTrials(ippant,itrial,3) = usedata.PFI_num1target_unadjusted_BL;
                Freq_LocPFI_acrossTrials(ippant,itrial,4) = usedata.PFI_num1target_unadjusted_BR;
                
                mDurperLocPFI_acrossTrials(ippant,itrial,1)=usedata.meanPFIduration_perdisap_1target_TL;
                mDurperLocPFI_acrossTrials(ippant,itrial,2)=usedata.meanPFIduration_perdisap_1target_TR;
                mDurperLocPFI_acrossTrials(ippant,itrial,3)=usedata.meanPFIduration_perdisap_1target_BL;
                mDurperLocPFI_acrossTrials(ippant,itrial,4)=usedata.meanPFIduration_perdisap_1target_BR;
                
                
                totalDurperLocPFI_acrossTrials(ippant,itrial, 1)=usedata.PFItotalduration_disap_1target_TL;
                totalDurperLocPFI_acrossTrials(ippant,itrial, 2)=usedata.PFItotalduration_disap_1target_TR;
                totalDurperLocPFI_acrossTrials(ippant,itrial, 3)=usedata.PFItotalduration_disap_1target_BL;
                totalDurperLocPFI_acrossTrials(ippant,itrial, 4)=usedata.PFItotalduration_disap_1target_BR;
                
                %Now determine correct frequency for each location, before storage.
                %rank locations.
                locswas = [usedata.TL_Freq, usedata.TR_Freq, usedata.BL_Freq, usedata.BR_Freq];
                
                %8Hz first, find the index in Location order, for
                %appropriate frequency.
                for ifreq=1:4
                    switch ifreq
                        case 1
                            ind= find(locswas==8);
                        case 2
                            ind= find(locswas==13);
                        case 3
                            ind= find(locswas==15);
                        case 4
                            ind= find(locswas==18);
                    end
                    
                    Freq_FreqPFI_acrossTrials(ippant,itrial,ifreq) = Freq_LocPFI_acrossTrials(ippant,itrial,ind);
                    mDurperFreqPFI_acrossTrials(ippant,itrial,ifreq)=mDurperLocPFI_acrossTrials(ippant,itrial,ind);
                    totalDurperFreqPFI_acrossTrials(ippant,itrial, ifreq)=totalDurperLocPFI_acrossTrials(ippant,itrial, ind);
                    
                    
                    
                    %now the combined
                    Freq_HzxLocPFI_acrossTrials(ippant,itrial,ifreq, ind) = Freq_LocPFI_acrossTrials(ippant,itrial,ind);
                    mDurper_HzxLocPFI_acrossTrials(ippant,itrial,ifreq, ind)=mDurperLocPFI_acrossTrials(ippant,itrial,ind);
                    totalDur_HzxLocPFI_acrossTrials(ippant,itrial,ifreq, ind)=totalDurperLocPFI_acrossTrials(ippant,itrial,ind);
    
                    
                    
                end
            end
        end
    end
    
    save('PFI_data_concat', 'Freq_FreqPFI_acrossTrials',...
        'Freq_LocPFI_acrossTrials',...
        'mDurperFreqPFI_acrossTrials',...
        'mDurperLocPFI_acrossTrials',...
        'totalDurperFreqPFI_acrossTrials',...
        'totalDurperLocPFI_acrossTrials',...
        'Freq_HzxLocPFI_acrossTrials',...
            'mDurper_HzxLocPFI_acrossTrials',...
                'totalDur_HzxLocPFI_acrossTrials',...
                '-append')
    
    
end
%% based on IG's code (reshuffle_analysis)
if job.createShufflePFIdata_pernum==1
    load('MD_AllBP_perppant.mat')
    nshuff=1000;
    allRandomAllPP =zeros(length(allppants),nshuff,5,3600);
    for ippant=1:length(allppants)
        goodPP=allppants(ippant);
        
        cd(basedir)
        
        
        
        BPdata=ppantTrialDatawithDetails(goodPP).AllBPData;
        
        allRandom = zeros(nshuff,5, 3600);
        for countRand=1:nshuff
            %index of four random trials.
            randTrial=randi([1 24],1,4);
            %make sure none are the same:
%             while length(unique(randTrial))<4
%                 randTrial=randi([1 24],1,4);
%             end
            
            %take a location on screen from separate trials, and concatenate.
            plotTable(1,:)=BPdata(randTrial(1),1,:);
            plotTable(2,:)=BPdata(randTrial(2),2,:);
            plotTable(3,:)=BPdata(randTrial(3),3,:);
            plotTable(4,:)=BPdata(randTrial(4),4,:);
            %accumulative sum
            plotTable(5,:)=sum(plotTable(1:4,:),1);
            
            allRandom(countRand,:,:)= plotTable;
            
        end
        
        allRandomAllPP(ippant,:,:,:)=allRandom;
        
        
    end
    
    save('ShuffledData', 'allRandomAllPP', '-V7.3')
    
    
    %%
end
if job.calcPFIdataperNum_shuffled==1 %resaves into ShuffledData
    %%
    cd(basedir)
    load('ShuffledData')
    for ippant = 1:size(allRandomAllPP,1)
        
        
        for itrial=1:nshuff
            
            accumBP = squeeze(allRandomAllPP(ippant,itrial,5,:))'; %combine each locations BP.
            
            
                
            if excludeTransientlength>0;%; %ms
                % then adjust for short transient button presses:
                
%                 figure(1);
%                 plot(accumBP, 'k'); %plot original.
                
                % find all diffs
                sws = find(diff(accumBP));
                
                %how long between BPs?
                diff_sws = diff(sws);
                
                rmv = find(diff_sws<excludeTransientlength);
                % now work through and replace.
                for irm= rmv
                    
                    tremovesw =sws(irm);
                    %replace values for this period:
                    vals = accumBP(1,tremovesw);
                    %replace to = sws
                    trep = sws(irm+1);
                    %replace now:
                    accumBP(1,tremovesw:trep)=vals;
                    
                    
                end
                
%                 hold on;
%                 plot(accumBP, 'r'); %check method.
%   PFI_only(ippant).Trial(itrial).newaccumBPs= accumBP;
            end
            
            
            
            dur0Disap=0;
            dur1Disap=0; %total time PFI
            dur2Disap=0;
            dur3Disap=0; %3 or more
            
            ZeroTargetDisap=0;
            
            OneTargetDisap=0;
            TwoTargetDisap=0;
            ThreeormoreTargetDisap=0;
            
            
            framestamp0_begin=[];
            framestamp0_end=[];
            
            framestamp1_begin=[];
            framestamp1_end=[];
            framestamp2_begin=[];
            framestamp2_end=[];
            framestamp3_begin=[];
            framestamp3_end=[];
            
            
            
            
            
            for timeDur=2:(length(accumBP)-1)
                
                if accumBP(1,timeDur)==1
                    
                    
                    if accumBP(1,timeDur-1)~=1 %single button onset
                        if framestamp1_begin>0
                            framestamp1_begin= [framestamp1_begin, timeDur];
                        else
                            framestamp1_begin= timeDur;
                        end
                        OneTargetDisap=OneTargetDisap+1;
                        
                    end
                    
                    
                    if accumBP(1,timeDur+1)~=1 || timeDur==length(accumBP)-1
                        %singlebutton offset
                        
                        if framestamp1_end>0
                            framestamp1_end= [framestamp1_end, timeDur+1];
                        else
                            framestamp1_end= timeDur+1;
                        end
                    end
                    
                    
                elseif accumBP(1,timeDur)==2
                    
                    
                    
                    if accumBP(1,timeDur-1)~=2 %double button onset
                        if framestamp2_begin>0
                            framestamp2_begin= [framestamp2_begin, timeDur];
                        else
                            framestamp2_begin= timeDur;
                        end
                        TwoTargetDisap=TwoTargetDisap+1;
                        
                    end
                    
                    if accumBP(1,timeDur+1)~=2 || timeDur==length(accumBP)-1
                        %singlebutton offset
                        
                        if framestamp2_end>0
                            framestamp2_end= [framestamp2_end, timeDur+1];
                        else
                            framestamp2_end= timeDur+1;
                        end
                    end
                    
                    
                elseif  accumBP(1,timeDur)==3
                    
                    
                    
                    if accumBP(1,timeDur-1)~=3 %
                        if framestamp3_begin>0
                            framestamp3_begin= [framestamp3_begin, timeDur];
                        else
                            framestamp3_begin= timeDur;
                        end
                        ThreeormoreTargetDisap=ThreeormoreTargetDisap+1;
                        
                    end
                    
                    if accumBP(1,timeDur+1)~=3 || timeDur==length(accumBP)-1
                        %singlebutton offset
                        
                        if framestamp3_end>0
                            framestamp3_end= [framestamp3_end, timeDur+1];
                        else
                            framestamp3_end= timeDur+1;
                        end
                    end
                    
                    
                elseif accumBP(1,timeDur)==0 % now included for new figs.
                    if accumBP(1,timeDur-1) ~=0 %start of a disappearance.
                        if framestamp0_begin>0
                            framestamp0_begin = [framestamp0_begin, timeDur];
                        else
                            framestamp0_begin = timeDur;
                        end
                        ZeroTargetDisap=ZeroTargetDisap+1;
                    end
                    
                    if accumBP(1,timeDur+1)~=0 || timeDur==length(accumBP)-1 %accounting for end of trial.
                        %singlebutton offset
                        
                        if framestamp0_end>0
                            framestamp0_end= [framestamp0_end, timeDur+1];
                        else
                            framestamp0_end= timeDur+1;
                        end
                        
                    end
                    
                    
                end
                
            end
            
            %in case a button was held until the end of a trial, we
            %need to remove these last disapperances
            
            %adjust for trial beginning.
            if accumBP(1,1)==0
            framestamp0_begin = [1 framestamp0_begin];
            
                ZeroTargetDisap = length(framestamp0_begin);
            
            end
            
            if length(framestamp0_begin) > length(framestamp0_end)
                %adjust length
                framestamp0_begin= framestamp0_begin(1:length(framestamp0_end));
                
            end
            
            if length(framestamp1_begin) > length(framestamp1_end)
                %adjust length
                framestamp1_begin= framestamp1_begin(1:length(framestamp1_end));
                OneTargetDisap = length(framestamp1_begin);
            end
            %same for 2 and three target instances.
            if  length(framestamp2_begin) > length(framestamp2_end)
                
                framestamp2_begin= framestamp2_begin(1:length(framestamp2_end));
                TwoTargetDisap = length(framestamp2_begin);
                
            end
            
            if  length(framestamp3_begin) > length(framestamp3_end)
                
                framestamp3_begin= framestamp3_begin(1:length(framestamp3_end));
                ThreeormoreTargetDisap = length(framestamp3_begin);
            end
            
            
            
            
            
            alldurs0 = framestamp0_end - framestamp0_begin;
            alldurs1 = framestamp1_end - framestamp1_begin;
            alldurs2 = framestamp2_end - framestamp2_begin;
            alldurs3= framestamp3_end - framestamp3_begin;
            
            dur0Disap= sum(alldurs0);
            dur1Disap = sum(alldurs1);
            dur2Disap = sum(alldurs2);
            dur3Disap = sum(alldurs3);
            
              for in=1:4
                    switch in
                        case 1
                            ddurs=alldurs0;
                            nPFI=ZeroTargetDisap;
                            frame_end=framestamp0_end;
                            frame_begin=framestamp0_begin;
                        case 2
                            ddurs=alldurs1;
                            nPFI=OneTargetDisap;
                            frame_end=framestamp1_end;
                            frame_begin=framestamp1_begin;
                        case 3
                            ddurs=alldurs2;
                            nPFI=TwoTargetDisap;
                            frame_end=framestamp2_end;
                            frame_begin=framestamp2_begin;
                        case 4
                            ddurs=alldurs3;
                            nPFI=ThreeormoreTargetDisap;
                            frame_end=framestamp3_end;
                            frame_begin=framestamp3_begin;
                      
                    end
                    
                    keep= find(ddurs>=excludeTransientlength);
                    
                    %now adjust outgoings to only keep the PFI that were
                    %not transients.
                    durDisap = sum(ddurs(keep));
                    nPFI=length(keep);
                    alldurs = ddurs(keep);
                    
                    frame_begs = frame_begin(keep);
                    frame_ends = frame_end(keep);
                    %store
                    switch in
                        case 1
                            %mean duration per individual disap (in seconds).               
                            shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_0target = durDisap/nPFI/60;
                            %  total duration per trial (in seconds)
                            shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_0target = durDisap/60;
                            %frequency of PFI per trial (how many times PFI occurs)
                            shuffledPFI(ippant).Trial(itrial).PFI_num0target_unadjusted = nPFI;
                             %record time stamps (in frames) of these events.
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_0target_framestart= frame_begs;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_0target_frameend= frame_ends;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_0target_durs= alldurs;
                             
                        case 2
                            shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_1target = durDisap/nPFI/60;
                            shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_1target = durDisap/60;
                            shuffledPFI(ippant).Trial(itrial).PFI_num1target_unadjusted = nPFI;
                                shuffledPFI(ippant).Trial(itrial).PFI_disap_1target_framestart= frame_begs;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_1target_frameend= frame_ends;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_1target_durs= alldurs;
                        case 3
                            shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_2target = durDisap/nPFI/60;
                            shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_2target = durDisap/60;
                            shuffledPFI(ippant).Trial(itrial).PFI_num2target_unadjusted = nPFI;
                                shuffledPFI(ippant).Trial(itrial).PFI_disap_2target_framestart= frame_begs;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_2target_frameend= frame_ends;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_2target_durs= alldurs;
                        case 4
                            shuffledPFI(ippant).Trial(itrial).meanPFIduration_perdisap_3ormoretarget = durDisap/nPFI/60;
                            shuffledPFI(ippant).Trial(itrial).PFItotalduration_disap_3ormoretarget = durDisap/60;
                            shuffledPFI(ippant).Trial(itrial).PFI_num3ormoretarget_unadjusted = nPFI;
                                shuffledPFI(ippant).Trial(itrial).PFI_disap_3ormoretarget_framestart= frame_begs;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_3ormoretarget_frameend= frame_ends;
                             shuffledPFI(ippant).Trial(itrial).PFI_disap_3ormoretarget_durs= alldurs;
                       
                end
                end
            %collect
            shuffledPFI(ippant).Trial(itrial).Goodtrial=1;
            
            
            
        end
        
    end
    
    save('ShuffledData', 'shuffledPFI', '-append')
    
end

if job.concatPFIacrossPpants_num_shuffled==1
    cd(basedir)
    load('ShuffledData')
    
    Freq_NumPFI_acrossTrials_shuffled = nan(length(allppants),nshuff,4);
    mDurperNumPFI_acrossTrials_shuffled=nan(length(allppants),nshuff,4);
    totalDurperNumPFI_acrossTrials_shuffled=nan(length(allppants),nshuff,4);
    
    for ippant = 1:length(shuffledPFI)
        
        for itrial=1:nshuff
            
            if shuffledPFI(ippant).Trial(itrial).Goodtrial==1 %not rejected
                usedata=shuffledPFI(ippant).Trial(itrial);
                
                %zero targets first
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,1) = usedata.PFI_num0target_unadjusted;
                
                
                tmp=usedata.meanPFIduration_perdisap_0target;
                if isinf(tmp)
                    error('check')
                else
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,1)=usedata.meanPFIduration_perdisap_0target;
                end
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,1)=usedata.PFItotalduration_disap_0target;
                
                
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,2) = usedata.PFI_num1target_unadjusted;
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,2)=usedata.meanPFIduration_perdisap_1target;
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,2)=usedata.PFItotalduration_disap_1target;
                %double
                
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,3) = usedata.PFI_num2target_unadjusted;
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,3)=usedata.meanPFIduration_perdisap_2target;
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,3)=usedata.PFItotalduration_disap_2target;
                
                %'2 or more'
                
                Freq_NumPFI_acrossTrials_shuffled(ippant,itrial,4) = usedata.PFI_num3ormoretarget_unadjusted;
                mDurperNumPFI_acrossTrials_shuffled(ippant,itrial,4)=usedata.meanPFIduration_perdisap_3ormoretarget;
                totalDurperNumPFI_acrossTrials_shuffled(ippant,itrial,4)=usedata.PFItotalduration_disap_3ormoretarget;
                
            end
            
        end
        
%         mean(mDurperNumPFI_acrossTrials_shuffled(ippant,:,:),2)
    end
    
    save('ShuffledData', 'Freq_NumPFI_acrossTrials_shuffled', 'mDurperNumPFI_acrossTrials_shuffled', 'totalDurperNumPFI_acrossTrials_shuffled', '-append')
    
end

if job.concatPFIacrossPpants_slope_nulldistribution==1;
% cd(basedir)
%     load('ShuffledData')
%     
%     Slope_nulldistribution=nan(200,2, 4); %pearsons r and pval, x (freq, mdur and total)
%     
%     
%     for iDVtype=1:3
%             
%         
%         for itrial=1:200
%             
%             % for each shuff, plot the number of relevant data points.
%             
%             thisshuff=zeros(3,length(allppants));
%             
%             
%                 
%             for ippant = 1:length(allppants)% just good ppants
%                 
%                 usedata=shuffledPFI(ippant).Trial(itrial);
%                 
%                 switch iDVtype
%                     case 1 %frequency.
%                 thisshuff(1,ippant)=  usedata.PFI_num1target_unadjusted;
%                 thisshuff(2,ippant)=  usedata.PFI_num2target_unadjusted;
%                 thisshuff(3,ippant)=  usedata.PFI_num3ormoretarget_unadjusted;
%                     case 2 % dur per.
%                         
%                 thisshuff(1,ippant)=  usedata.meanPFIduration_perdisap_1target;
%                 thisshuff(2,ippant)=  usedata.meanPFIduration_perdisap_2target;
%                 thisshuff(3,ippant)=  usedata.meanPFIduration_perdisap_3ormoretarget;                        
%                         
%                     case 3 %total dur.
%                         thisshuff(1,ippant)=   usedata.PFItotalduration_disap_1target;
%                 thisshuff(2,ippant)=  usedata.PFItotalduration_disap_2target;
%                 thisshuff(3,ippant)=    usedata.PFItotalduration_disap_3ormoretarget;
%                         
%                         
%                 end
%                 
%             end
%                
%             % now calculate slope for this shuffle.
%             xvect= [ones(1,length(allppants)) ; ones(1,length(allppants))*2;ones(1,length(allppants))*3];
%             scatter(xvect,thisshuff);
%             shg
%             shg
%         end
%         
%     end
%     
%     save('ShuffledData', 'Freq_NumPFI_acrossTrials_shuffled', 'mDurperNumPFI_acrossTrials_shuffled', 'totalDurperNumPFI_acrossTrials_shuffled', '-append')
%     
    
    
    
    
    
end

%% Plot Behavioural output.

if job.plotBehaviouraldata ==1
    %%
    cd(basedir)
    load('PFI_data_concat')
    cd('Figures')
    cd('Behavioural Bar Plots')
    icounter=1;
    clf
    plotZero=0;
    for iIV=3%1:2
        switch iIV
            case 1 %look at frequency of target flicker
                p1= Freq_FreqPFI_acrossTrials;
                p2= mDurperFreqPFI_acrossTrials;
                p3= totalDurperFreqPFI_acrossTrials;
                xlabelis = 'Flicker frequency of absent target';
                xticks = [{'8Hz'} {'13Hz'}, {'15Hz'}, {'18Hz'}];
                
            case 2 %look at mean duration of disappearance per trial.
                
                p1= Freq_LocPFI_acrossTrials;
                p2= mDurperLocPFI_acrossTrials;
                p3= totalDurperLocPFI_acrossTrials;
                xlabelis = 'Location of absent target';
%                 xticks = [{'TL'} {'TR'}, {'BL'}, {'BR'}];
               xticks = [{'Left'} {'Right'}];
            case 3
                p1= Freq_NumPFI_acrossTrials;
                p2= mDurperNumPFI_acrossTrials;
                p3= totalDurperNumPFI_acrossTrials;
                xlabelis = '# of invisible targets';
                if plotZero==1
                xticks = [{'0'} {'1'}, {'2'}, {'3 or 4'}];
                else
                xticks = [{'1'}, {'2'}, {'3 or 4'}];
                end
        end

        for iDV=1:3
            %take mean within, then across ppants.
            %within
            
            subplot(2,3,icounter)
            switch iDV
                
                case 1
                    datatoplot=p1;
                    yis={['PFI per trial']};
               
                case 2
                    datatoplot=p2;
                    yis={['PFI duration [s]']};
                case 3
                    
                    datatoplot=p3;
                    yis= 'Total duration [s]';
                    
            end
            
            if plotZero~=1
                datatoplot=datatoplot(:,:,2:4);
            end
            %only look at good ppants.
            datanow1=squeeze(datatoplot(allppants,:,:));
            datanow=squeeze(nanmean(datanow1,2));
            
            % calculate error bars according to 
            %confidence interval for within subj designs. (Cousineau,
            %2005)
            % y = x - mXsub + mXGroup,
            x = datanow;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. RMfactors)
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1,size(x,2)) + repmat(mXgroup, size(x));
            
            % mean to plot
            mData = squeeze(nanmean(NEWdata,1));
            
            %compute stErr %which version?
            stErr = std(NEWdata)/sqrt(size(datanow,1));
            
            %
            if iIV==1 ||iIV==3
            bh=bar(mData);
            bh(1).FaceColor= [.5 .5 .9];
            bh(1).FaceColor= ['b'];
            hold on
            xis=1:length(mData);
            
            %                 title(['PFI data by ' xlabelis], 'fontsize', fontsize)
            elseif iIV==2 %stack Left vs Right and top and bottom.
              %%
                mData = [mData(1), mData(3); mData(2), mData(4)];
                stErr = [stErr(1), stErr(3); stErr(2), stErr(4)];
                bh=bar(mData);
            bh(1).FaceColor= ['m'];            
            bh(2).FaceColor= [.2 .2 .4];
            bh(2).FaceColor= ['b'];
            hold on
              offset=.15;
            xis=[1-offset, 1+offset; 2-offset, 2+offset];
            
            %%
            end
                errorbar(xis,mData,stErr, 'LineStyle', ['none'],'Color', 'k','linewidth', 2 )
                
                
            axis('tight')
           
                    ylabel([yis]);
                  xlabel([xlabelis]);
                  yt=get(gca, 'ylim');
                  set(gca, 'ytick', [0:round((ceil(yt(2)))/10):ceil(yt(2))]); 

            
%             axis tight

            
            if iIV==2
                xlim([.5 2.5])
                legend('Top', 'Bottom ')
            end
            
            
            set(gca, 'xticklabel', xticks)
            set(gca, 'fontsize', fontsize-2)
           icounter=icounter+1; 
        end
        printfilename = ['PFI data by ' xlabelis];
        set(gcf, 'color', 'w')
%         print('-dpng', printfilename)
%         print('-dpng', 'PFI data together, location interaction')
    end
    shg
end

if job.plotBehaviouraldata_num_with_shuffled==1
    %slim version of the above.
    cd(basedir)
    load('PFI_data_concat')
    load('ShuffledData')
    cd('Figures')
    cd('Behavioural Bar Plots')
    %%
    fontsize=30;
    clf
    plot4onlast=0; % change to plot 0-4 on the last panel.
    plot4onall=0;
    
    % new option to plot together or separate:
    %var 1, 2 ,3
    plotrealshufforcombined = 3;
    figure(1); clf
    for iIV=2
        switch iIV
            case 1 %look at frequency of target flicker
                p1= Freq_FreqPFI_acrossTrials;
                p2= mDurperFreqPFI_acrossTrials;
                p3= totalDurperFreqPFI_acrossTrials;
                xlabelis = 'Flicker frequency of absent target';
                
                xticks = [{'8Hz'} {'13Hz'}, {'15Hz'}, {'18Hz'}];
                
          
            case 2 %compare shuffled to number absent
                p1= Freq_NumPFI_acrossTrials;
                p2= mDurperNumPFI_acrossTrials;
                p3= totalDurperNumPFI_acrossTrials;
                xlabelis = 'nPFI';
                xticks = [{'1'} {'2'}, {'3 or 4'}];
                
                % allocate shuffled data.
                s1=Freq_NumPFI_acrossTrials_shuffled;
                s2= mDurperNumPFI_acrossTrials_shuffled;
                s3= totalDurperNumPFI_acrossTrials_shuffled;
                
                
            case 3 %compare location
                
                p1= Freq_LocPFI_acrossTrials;
                p2= mDurperLocPFI_acrossTrials;
                p3= totalDurperLocPFI_acrossTrials;
                xlabelis = 'Location of absent target';
                xticks = [{'TL'} {'TR'}, {'BL'}, {'BR'}];
        end
        clf
        for iDV=1:3
            %take mean within, then across ppants.
            %within
            
            subplot(2,3,iDV)
            switch iDV
                case 1
                    datatoplot=p1;
                    shuffledtoplot=s1;
                    yis='# of PFI events per trial';
                case 2
                    
                    datatoplot=p2;
                    shuffledtoplot=s2;
                    yis='PFI duration per event [s]';

                case 3
                                        datatoplot=p3;
                    shuffledtoplot=s3;
                    yis= 'Total duration [s]';
                    if plot4onlast==1 
                    xticks = [{'0'},{'1'} {'2'}, {'3 or 4'}];
                    end
            end
            if plot4onall==1
                 xticks = [{'0'},{'1'} {'2'}, {'3 or 4'}];
            end
                
            %%
            %only look at good ppants.
            datanowreal=squeeze(nanmean(datatoplot(allppants,:,:),2));
            datanowshuffled=squeeze(nanmean(shuffledtoplot,2));
            
            
            %across
            mDataReal=squeeze(nanmean(datanowreal,1));
            mDataShuff=squeeze(nanmean(datanowshuffled,1));
            %rearrange for plots
            
            if plot4onlast==1 && iDV==3
            mData = [mDataReal(1,1), mDataShuff(1,1); mDataReal(1,2), mDataShuff(1,2); mDataReal(1,3), mDataShuff(1,3);mDataReal(1,4), mDataShuff(1,4)];
            elseif plot4onall==1
                mData = [mDataReal(1,1), mDataShuff(1,1); mDataReal(1,2), mDataShuff(1,2); mDataReal(1,3), mDataShuff(1,3);mDataReal(1,4), mDataShuff(1,4)];
            else
                mData = [mDataReal(1,2), mDataShuff(1,2); mDataReal(1,3), mDataShuff(1,3); mDataReal(1,4), mDataShuff(1,4)];
            end
            
            
            %standard err m  %note different calculations for data and
            %shuffled.
            for id=1:2
                switch id
                    case 1
                        x=datanowreal;
                    case 2
                        x=shuffledtoplot;
                        
                end
                
                if id==1
                mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. RMfactors)
            mXgroup = mean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1,size(x,2)) + repmat(mXgroup, size(x));
            
            std_tmp = nanstd(NEWdata);
            %careful here, as we need to divide by the number of data
            %points
            stErr=zeros(size(std_tmp));
            for idv=1:length(std_tmp);
                %numberobs
                n_isnan = length(find(isnan(NEWdata(:,idv))));
                %correct for stErr.
                stErr(idv)= std_tmp(idv)/ sqrt(size(NEWdata,1)-n_isnan);
            end
            
                else
                    
                    %loop over all shuffles.
                    allstE = zeros(size(x,2), size(x,3));
                    for ishuff=1:size(x,2)
                        newx = squeeze(x(:,ishuff,:));
                       mXppant =squeeze( nanmean(newx,2)); %mean across conditions we are comparing (within ppant ie. RMfactors)
            mXgroup = mean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = newx - repmat(mXppant, 1,size(newx,2)) + repmat(mXgroup, size(newx));
            
            std_tmp = nanstd(NEWdata);
            %careful here, as we need to divide by the number of data
            %points
            stErr=zeros(size(std_tmp));
            for idv=1:length(std_tmp);
                %numberobs
                n_isnan = length(find(isnan(NEWdata(:,idv))));
                %correct for stErr.
                stErr(idv)= std_tmp(idv)/ sqrt(size(NEWdata,1)-n_isnan);
            end
%             stErr = nanstd(NEWdata)/sqrt(size(newx,1));
            
            
            allstE(ishuff, :) = stErr; 
                    end
                
                    %take mean across shuffles.
                    stErr = nanmean(allstE,1);
                end
                    
            switch id
                case 1
                    stErrReal = stErr;
                case 2
                    stErrShuffled= stErr;
            end
            end
            
            
            
            
            if plot4onlast==1 && iDV==3
            stErr=[stErrReal(1), stErrShuffled(1);stErrReal(2), stErrShuffled(2);stErrReal(3), stErrShuffled(3) ;stErrReal(4), stErrShuffled(4) ];
            
            elseif plot4onall==1
                stErr=[stErrReal(1), stErrShuffled(1);stErrReal(2), stErrShuffled(2);stErrReal(3), stErrShuffled(3) ;stErrReal(4), stErrShuffled(4) ];
            
            else
                stErr=[stErrReal(2), stErrShuffled(2);stErrReal(3), stErrShuffled(3) ;stErrReal(4), stErrShuffled(4) ];
                
            end
            
            % adjust x axis as necessary
            offset=.15;
            xis = [1,2,3,4];
            
            %if we have both data types at once, and plotting 4 bars
            if (plot4onlast==1 && iDV==3) || plot4onall==1
            xis=[1-offset, 1+offset; 2-offset, 2+offset; 3-offset,3+offset; 4-offset,4+offset];
            end
            %else remaining types
            if plot4onall==1 && plotrealshufforcombined<3
                xis=1:4;
            elseif plot4onall~=1 && plotrealshufforcombined==3
                 xis=[1-offset, 1+offset; 2-offset, 2+offset; 3-offset,3+offset];
            end
            %
            if plot4onall~=1 && plotrealshufforcombined<3
               xis=[1:3]; 
            end
            
            
                
                       
            
            
            
            %
            if plotrealshufforcombined<3
             mData=mData(:,plotrealshufforcombined);
            stErr = stErr(:,plotrealshufforcombined);
            end
            
           
            bh=bar(mData);
            hold on            
            errorbar(xis,mData,stErr,'linestyle', ['none'], 'linewidth', 2 , 'Color', 'k')
            switch plotrealshufforcombined
                case 1
                    bh.FaceColor=['b'];
                    bh.BarWidth = .2;
                case 2
                    bh.FaceColor= [.5 .5 .5];           
                    bh.BarWidth = .2;
            
                case 3
            bh(1).FaceColor= ['b'];           
            bh(2).FaceColor= [.5 .5 .5];
            if iDV==3
                
            legend 'Observed' 'Shuffled'
            end
            
            end
            
            %%
            
            axis('tight')            
            if plot4onall==1
            xlim([.5 4.5])
            else
                xlim([.5 3.5])
            end
            
            xlabel(xlabelis)
            ylabel(yis)
             yt=get(gca, 'ylim');
             ylim([0 yt(2)*1.1])
             
             if round((ceil(yt(2)))/10) >1
                  set(gca, 'ytick', [0:round((ceil(yt(2)))/10):ceil(yt(2))]); 
             else
                  set(gca, 'ytick', [0:1:ceil(yt(2))]); 
             end

            set(gca, 'xticklabel', xticks)
            set(gca, 'fontsize', fontsize-2)
            

            
            
            end
    end
    %%
     set(gcf, 'color', 'w')
        shg
        %%
       
        if plot4onlast==1 || plot4onall==1
        printfilename = ['PFI data by ' xlabelis ' shuffled'];
        else
            printfilename = ['PFI data by ' xlabelis ' shuffled, no 0 PFI'];
        end
            
        print('-dpng', printfilename)
end
    
    
%
if job.compareSlopesofShuffledvsObserveddata==1
    
   %start with slope across  3?
     cd(basedir)
    load('PFI_data_concat')
    load('ShuffledData')
    cd('Figures')
    cd('Behavioural Bar Plots')
%     clf
   testnPFIinclZERO=0;
   
                p1= Freq_NumPFI_acrossTrials;
                p2= mDurperNumPFI_acrossTrials;
                p3= totalDurperNumPFI_acrossTrials;
                xlabelis = '# of invisible targets';
                xticks = [{'1'} {'2'}, {'3 or 4'}];
                
                % allocate shuffled data.
                s1=Freq_NumPFI_acrossTrials_shuffled;
                s2= mDurperNumPFI_acrossTrials_shuffled;
                s3= totalDurperNumPFI_acrossTrials_shuffled;
                
         
        pvalswere=[];
        figure(1);
        for iDV=1:3
            %take mean within, then across ppants.
            %within
            
            subplot(2,3,iDV)
            switch iDV
                case 1
                    datatoplot=p1;
                    shuffledtoplot=s1;
                    yis='PFI per trial';
                case 2
                    
                    datatoplot=p2;
                    shuffledtoplot=s2;
                    yis='PFI duration [s]';

                case 3
                                        datatoplot=p3;
                    shuffledtoplot=s3;
                    yis= 'Total duration [s]';
                    if testnPFIinclZERO==1 
                    xticks = [{'0'},{'1'} {'2'}, {'3 or 4'}];
                    end
            end
            if testnPFIinclZERO==1
                 xticks = [{'0'},{'1'} {'2'}, {'3 or 4'}];
                 usedata=1:4;
                 xrange=usedata;
            else
                usedata=2:4;
                xrange=1:3;
            end
            %%
            %only look at good ppants.
            pl=[];
            for iregres=1:2
                if iregres==1
             datanowreal=squeeze(nanmean(shuffledtoplot(:,:,usedata),2));
                    Bcolor=[.5 .5 .5];
                    dcolor= [.5 .5 .5];
                else
                   
            datanowreal=squeeze(nanmean(datatoplot(allppants,:,usedata),2));
            Bcolor='b';
            dcolor='b';
                end
                    
            
            %% observed slope = 
            scAx = 1:size(datanowreal,2);
            if iregres==1
                scAx=scAx+.15;
            end
            scAx = repmat(scAx, [size(datanowreal,1),1]);
            
            %reshape into vector for plot
            scAxm = reshape(scAx, 1, size(scAx,1)*size(scAx,2));
            datanowSc = reshape(datanowreal, 1, size(scAxm,1)*size(scAxm,2));
        
            
            %% find index for NaN, so to not use in polyfit.
            exclpoly=find(isnan(datanowSc));
            tmp=1:length(scAxm); 
            tmp(exclpoly)=[];
            
            %remove nans from data for polyfit.
            scAxm=scAxm(tmp);
            datanowSc=datanowSc(tmp);
            %%
     
        
%         take linear regression and plot coeff. and slope
        [p, S] = polyfit(scAxm, datanowSc, 1); %linear fit
%         p1 is the slope, p2 the intercept
        
        if iregres==2 % store the observed slope
        ObservedSlope=p(1);
        end
        
        f1=polyval(p,scAxm);
        hold on; 
        
        %%%%%%%%% plotting
               
%         sc=    scatter(scAxm,datanowSc);
%         sc.MarkerFaceColor=dcolor;
%         sc.MarkerEdgeColor=dcolor;
%         
%         xlabel(xlabelis);
%         xlim([xrange(1)-.5, xrange(end)+.5])
%         
%         ylabel(yis);
%         pl(iregres)=plot(scAxm, f1, 'color',Bcolor, 'linew', 6);
%         
%         if iregres==2 && iDV==3
% %         legend(pl, ['\beta' ' = ' sprintf('%.2f',p(1))])
%         legend([pl(2), pl(1)], {'Observed data', 'Shuffled data'})
%         
%         end
%         
%         % now calculate all slopes from shuffled data:
%         set(gca, 'fontsize', 15)
%         set(gca, 'xtick', xrange, 'xticklabels', xticks)
            end
        %%
        subplot(2,3,iDV+3)
        allslopes = zeros(1,size(shuffledtoplot,2));
        for ishuff=1:size(shuffledtoplot,2)
        
            datanow=squeeze(shuffledtoplot(:,ishuff,usedata));
            
            %% observed slope = 
            scAx = 1:size(datanow,2);
            scAx = repmat(scAx, [size(datanow,1),1]);
            
            %reshape into vector for plot
            scAxm = reshape(scAx, 1, size(scAx,1)*size(scAx,2));
            datanowSc = reshape(datanow, 1, size(scAxm,1)*size(scAxm,2));
        
        %we do need to remove Nans to use polyfit however:
        nanid=find(isnan(datanowSc));
        scAxm(nanid)=[];
        datanowSc(nanid)=[];
            
        %take linear regression and plot coeff. and slope
        [p, S] = polyfit(scAxm, datanowSc, 1); %linear fit
        %p1 is the slope, p2 the intercept
        
%         f1=polyval(p,scAxm);
%         hold on; plot(scAxm, f1)
        allslopes(ishuff)=p(1);
        end
        %plot all.
        shg
        %%
        hb=histogram(sort(allslopes),25);        
        hb.FaceColor=[.5 .5 .5];
         
        % fit CDF
        %note that if we have positive and negative numbers, this will
        %error.
        if any(hb.Data>0) && any(hb.Data<0)
            % shift all values so CDF is correct
            hb4cdf= hb.Data+1000;
            
        else 
            hb4cdf=hb.Data;
        end
        cdf= cumsum(hb4cdf)/ sum(hb4cdf);
        %the X values (actual CV) corresponding to .05
        
        [~,cv05uncorr] = (min(abs(cdf-.975)));
        
        
          [~, c2] = min(abs(hb.Data-ObservedSlope)); %closest value. 
                        pvalis= 1-cdf(c2);
        
        % height of plot?
        yt=get(gca, 'ylim');
        %plot
        hold on; 
        p05=plot([hb.Data(cv05uncorr) hb.Data(cv05uncorr)], [0, yt(2)*.7], ['k:'], 'linew', 3);
                
        text(hb.Data(cv05uncorr), yt(2)*.75, '97.5%', 'fontsize', 15, 'color', 'k');
                
        po=plot([ObservedSlope,ObservedSlope], [0, yt(2)*.7], ['-b'], 'linew', 3);
        


        if iDV==3
%         legend([po, p05], [{'Observed \bf\beta'},{ '95% CV'}])
        legend([po], [{'Observed'}])
        end



% xlabel(['\beta values from shuffled data']);
xlabel(['slope of the linear fit']);
           ylabel('count')
        
           %adjust xlim.
mX=min(hb.Data); maxX= max([ObservedSlope,hb.Data(cv05uncorr)]);
 xlim([ mX maxX + (maxX - mX)*.2])

           
           
           set(gca, 'fontsize', 23)
           pvalswere(iDV) = pvalis;
        end
        %%
        set(gcf, 'color', 'w')
        
        disp(['p values are ' num2str(pvalswere )])
    %%
%  print('-dpng', 'linearregression and shuffle')
end
    
    


%% Now the stats for the PFI data
if job.RMANOVA_PFI_byNum==1
    cd(basedir)
    %load relevant data
    load('PFI_data')
    
    for itest = 1:9
        clearvars 'DV*' 'subjects' 'factorarray'
        switch itest
            case 1
                data= Freq_NumPFI_acrossTrials;
                DVis='Instances of PFI';
                IV = 'Number_absent';
            case 2
                data= mDurperNumPFI_acrossTrials;
                DVis='Mean duration per PFI';
                IV = 'Number_absent';
            case 3
                data= totalDurperNumPFI_acrossTrials;
                DVis = 'Total PFI per trial';
                IV = 'Number_absent';
            case 4
                data= Freq_FreqPFI_acrossTrials;
                DVis='Instances of PFI';
                IV = 'Frequency_flicker';
            case 5
                
                data= mDurperFreqPFI_acrossTrials;
                DVis='Mean duration per PFI';
                IV = 'Frequency_flicker';
            case 6
                data= totalDurperFreqPFI_acrossTrials;
                DVis = 'Total PFI per trial';
                IV = 'Frequency_flicker';
            case 7
                data= Freq_LocPFI_acrossTrials;
                DVis='Instances of PFI';
                IV = 'Location_absent';
            case 8
                
                data= mDurperLocPFI_acrossTrials;
                DVis='Mean duration per PFI';
                IV = 'Location_absent';
            case 9
                data= totalDurperLocPFI_acrossTrials;
                DVis = 'Total PFI per trial';
                IV = 'Location_absent';
        end
        
        
        %% take mean within ppants,
        testd= squeeze(nanmean(data(goodppants,:,:),2));
        %reshape to single column for JvB rmanova
        DV = reshape(testd, [size(testd,1)*size(testd,2),1]);
        %JvB rmanova can't handle nan, so find and replace with zero.
        nanis = find(isnan(DV));
        DV(nanis)=0;
        
        %create identifers
        factorarray = zeros(size(DV,1),1);
        
        %first denote ppants
        nppants = 1:size(testd,1);
        % place in column
        subjects=repmat(nppants', [size(testd,2), 1]);
        
        %factor array for conds (can be 3 or 4 depending on IV).
        for icond=1:size(testd,2)
            %%
            conds= zeros(1, length(nppants))' + icond;
            %place correctly
            rows=nppants + (length(nppants)*(icond-1));
            factorarray(rows,1)=conds;
        end
        %%
        
        %     rmanova_return = rmanova(data,factorarray,subjects [,varnames] [,btw_ss_col])
        rmanova_results = rmanova(DV, factorarray, subjects, { num2str(IV)});
        
        switch itest
            case 1
                anova_FreqbyNumPFI=rmanova_results;
            case 2
                anova_mDurbyNumPFI=rmanova_results;
            case 3
                anova_totalDurbyNumPFI=rmanova_results;
            case 4
                anova_FreqbyFreqPFI=rmanova_results;
            case 5
                anova_mDurbyFreqPFI=rmanova_results;
            case 6
                anova_totalDurbyFreqPFI=rmanova_results;
            case 7
                anova_FreqbyLocPFI=rmanova_results;
            case 8
                anova_mDurbyLocPFI=rmanova_results;
            case 9
                anova_totalDurbyLocPFI=rmanova_results;
        end
    end
    save('PFI_data_rmANOVAresults',...
        'anova_FreqbyFreqPFI',...
        'anova_FreqbyLocPFI',...
        'anova_FreqbyNumPFI',...
        'anova_mDurbyFreqPFI',...
        'anova_mDurbyLocPFI',...
        'anova_mDurbyNumPFI',...
        'anova_totalDurbyFreqPFI',...
        'anova_totalDurbyLocPFI',...
        'anova_totalDurbyNumPFI')
    
    
end



if job.RMANOVA_PFI_byNum_vs_shuffled==1
    cd(basedir)
    %load relevant data
    load('PFI_data')
    load('Shuffleddata')
    load('MD_AllBP_perppant.mat', 'goodppants')
    for itest = 1:3
        clearvars 'DV*' 'subjects' 'factorarray'
        switch itest
            case 1
                data1= Freq_NumPFI_acrossTrials;
                data2=Freq_NumPFI_acrossTrials_shuffled;
                DVis='Instances of PFI';
                IV = 'Number_absent';
            case 2
                data1= mDurperNumPFI_acrossTrials;
                data2= mDurperNumPFI_acrossTrials_shuffled;
                DVis='Mean duration per PFI';
                IV = 'Number_absent';
            case 3
                data1= totalDurperNumPFI_acrossTrials;
                data2= totalDurperNumPFI_acrossTrials_shuffled;
                DVis = 'Total PFI per trial';
                IV = 'Number_absent';
            
        end
        
        
        %% take mean within ppants,
        testd_real= squeeze(nanmean(data1(goodppants,:,:),2));
        testd_shuff= squeeze(nanmean(data2,2)); %only have shuffled for goodppants.
        
        %reshape to single column for JvB rmanova
        DV1 = reshape(testd_real, [size(testd_real,1)*size(testd_real,2),1]);
        DV2 = reshape(testd_shuff, [size(testd_shuff,1)*size(testd_shuff,2),1]);
        DV=[DV1;DV2];
        
        %JvB rmanova can't handle nan, so find and replace with zero.
        nanis = find(isnan(DV));
        DV(nanis)=0;
        
        %create identifers
        factorarray = zeros(size(DV,1),2);
        
        %first denote ppants
        nppants = 1:size(testd_real,1);
        % place in column
        subjects=repmat(nppants', [size(testd_real,2)*2, 1]);
        
        %factor array for conds . First column for number absent.
        
        for icond=1:size(testd_real,2)
            %%
            conds= zeros(1, length(nppants))' + icond;
            %place correctly
            rows=nppants + (length(nppants)*(icond-1));
            factorarray(rows,1)=conds;
        end
        %%
        %second column for datatype
        factorarray(1:78,2)=ones(78,1);
        
        %extend (rough hack)
        factorarray(79:end,1)=factorarray(1:78,1);
        factorarray(79:end,2)=[ones(78,1)]*2;
        %%
        
        %     rmanova_return = rmanova(data,factorarray,subjects [,varnames] [,btw_ss_col])
        rmanova_results = rmanova(DV, factorarray, subjects, { num2str(IV) 'datatype'}, 2);
        
        switch itest
            case 1
                anova_FreqbyNumPFI_vsshuffled=rmanova_results;
            case 2
                anova_mDurbyNumPFI_vsshuffled=rmanova_results;
            case 3
                anova_totalDurbyNumPFI_vsshuffled=rmanova_results;
           
        end
    end
    save('PFI_data_rmANOVAresults',...        
        'anova_FreqbyNumPFI_vsshuffled',...     
        'anova_mDurbyNumPFI_vsshuffled',...
        'anova_totalDurbyNumPFI_vsshuffled', '-append')
%     
%     
end





%% Supplementary analyses.



if job.plotPFIdata_asfunctionofTrialNumber==1
    load('PFI_data')
    for datatype=1:3
        switch datatype
            case 1
                datais= Freq_NumPFI_acrossTrials;
                titleis= 'Instances of PFI per trial, across participants';
                yvec = [0 15];
                ylab = 'Number per trial (#)';
            case 2
                datais= mDurperNumPFI_acrossTrials;
                titleis= 'mean duration of each PFI instance, each trial, across participants';
                ylab = 'Individual PFI duration per trial (s)';
                
            case 3
                datais= totalDurperNumPFI_acrossTrials;
                titleis= 'Total Duration of PFI, each trial, across participants';
                ylab = 'Total PFI duration per trial (s)';
        end
        
        % squeeze per num of targets missing.
        clf
        
        for numt=1:3
            switch numt
                case 1
                    col= 'b';
                    caseis = ' for single targets';
                case 2
                    col= 'k';
                    
                    caseis = ' for two targets';
                case 3
                    col='r';
                    caseis = ' for three (or more) targets';
            end
            %%
            
            plotme = squeeze(datais(goodppants,:,numt));
            errplotme = nanstd(plotme)/sqrt(size(plotme,1));
            %               plot([itrial itrial], nanmean(plotme), ['r' '*'])
            %               bar(itrial, nanmean(plotme))
            hold on
            errorbar(1:24, nanmean(plotme), errplotme , [ col '-*']);
            
        end
        %%
        xt=set(gca, 'xtick', [1:24]);
        xlabel('Trial Number', 'fontsize', fontsize)
        title(titleis, 'fontsize', fontsize);
        ylim(yvec)
        ylabel(ylab, 'fontsize', fontsize)
        legend 'for single target' 'for two targets' 'for three (or more) targets')
        set(gca, 'fontsize', fontsize)
        %%
        print('-dpng', [titleis ])
        %%
    end
end
%%
if job.concatPpantDurationsbyPFI_num==1
    cd(basedir)
    load('PFI_data')
    PFIdurations_struct=[];
    counter=1;
    allsingle=[];
    alldouble=[];
    alltriple=[];
    for ippant = allppants
        
        ppantstruct=[];
        for itrial = 1:24
            %calculate durations and save
            if itrial==1
                ppantstruct.singleDisaps = PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs;
                ppantstruct.doubleDisaps =  PFI_only(ippant).Trial(itrial).PFI_disap_2target_durs;
                ppantstruct.tripleDisaps = PFI_only(ippant).Trial(itrial).PFI_disap_3ormoretarget_durs;
                
            else
                ppantstruct.singleDisaps = [ppantstruct.singleDisaps, PFI_only(ippant).Trial(itrial).PFI_disap_1target_durs];
                ppantstruct.doubleDisaps = [ppantstruct.doubleDisaps, PFI_only(ippant).Trial(itrial).PFI_disap_2target_durs];
                ppantstruct.tripleDisaps = [ppantstruct.tripleDisaps, PFI_only(ippant).Trial(itrial).PFI_disap_3ormoretarget_durs];
                
            end
            
            allsingle = [ allsingle , ppantstruct.singleDisaps];
            alldouble = [ alldouble , ppantstruct.doubleDisaps];
            alltriple = [ alltriple , ppantstruct.tripleDisaps];
            
        end
        PFIdurations_struct(counter).singleDisaps=ppantstruct.singleDisaps;
        PFIdurations_struct(counter).doubleDisaps=ppantstruct.doubleDisaps;
        PFIdurations_struct(counter).tripleDisaps=ppantstruct.tripleDisaps;
        counter=counter+1;
    end
    PFIdurations_allSingle=allsingle;
    PFIdurations_allDouble=alldouble;
    PFIdurations_allTriple=alltriple;
    %%
    save('PFI_data.mat', 'PFIdurations_struct', 'PFIdurations_allSingle','PFIdurations_allDouble', 'PFIdurations_allTriple' , '-append')
end
%%
if job.plothistofDisapdurationbynumbergone==1
    %plot histogram of disapperances by numer of targets PFI.
    cd(basedir)
    load('PFI_data.mat')
    %%
    figure();
    fontsize=10;
    for itarg=1:3
%         subplot(1,3,itarg)
hold on
        switch itarg
            case 1
                usedata=PFIdurations_allSingle/60;
                col='b';
                printname='nPFI 1 durs';
                
            case 2
                usedata=PFIdurations_allDouble/60;
                col='k';
                printname='nPFI 2 durs';
            case 3
                usedata=PFIdurations_allTriple/60;
                col='r';
                printname='nPFI 3 or more durs';
                
        end
        
        %%
        %round to nearest 100ms
        roundusedata = round(usedata,2);
%         clf
binedges=[0:.1:10];
        histogram(roundusedata,binedges,  'Facecolor', [col]);
%         title({[printname];['Across n=' num2str(length(PFIdurations_struct)) ' ppants']}, 'fontsize', fontsize)
title('all PFI durations');
        xlabel('Seconds', 'fontsize', fontsize)
        ylabel('Count', 'fontsize', fontsize)
%         legend(['max duration ' num2str(max(usedata)) 's'])
        set(gca, 'fontsize', fontsize)
        %%
        xlim([0 5])
        ylim([0 3000])
        %%
%         axis tight
        %
%         print('-dpng', printname)
    end
    legend('nPFI = 1', 'nPFI = 2', 'nPFI = 3 or more')
    set(gcf,'color', 'w'); box on;
    set(gca, 'fontsize', 20)
end
%%
if job.plothistofDisapdurationbynumbergone_perppants ==1
    load('PFI_data.mat')
    %%
      figure();
    clf    
    counter=1;
    for ippant= 1:length(allppants)
       
        subplot(6,4,ippant)
        
        
    for itarg=1:3
%         subplot(1,3,itarg)
hold on
        switch itarg
            case 1
                usedata=PFIdurations_struct(ippant).singleDisaps /60;
                col='b';
                printname='nPFI 1 durs';
                
            case 2
                usedata=PFIdurations_struct(ippant).doubleDisaps /60;
                col='k';
                printname='nPFI 2 durs';
            case 3
                usedata=PFIdurations_struct(ippant).tripleDisaps /60;
                col='r';
                printname='nPFI 3 or more durs';
                
        end
        
        %%
        %round to nearest 100ms
        roundusedata = round(usedata,2);
%         clf
binedges=[0:.1:10];
        histogram(roundusedata,binedges,  'Facecolor', [col]);
%         title({[printname];['Across n=' num2str(length(PFIdurations_struct)) ' ppants']}, 'fontsize', fontsize)
title('all PFI durations');
        xlabel('Seconds', 'fontsize', fontsize)
        ylabel('Count', 'fontsize', fontsize)
%         legend(['max duration ' num2str(max(usedata)) 's'])
        set(gca, 'fontsize', fontsize)
        %%
%         axis tight
        ylim([0 40]); xlim([0 5])
        %
%         print('-dpng', printname)
    end
%     legend('nPFI = 1', 'nPFI = 2', 'nPFI = 3 or more')
    set(gcf,'color', 'w'); box on;
    set(gca, 'fontsize', 12)
    title(['Participant ' num2str(ippant) ' data']);
        
    end
end
if job.exportDatatoExcelfolonDesktop==1;
    %%
    cd(basedir)
    %load relevant data
    load('PFI_data_concat')
    
    
    load('Shuffleddata')
    load('MD_AllBP_perppant.mat')
    
    % we have a number of different datatypes to append.
    % ppants x Hz
    % ppants x Loc
    % ppants x Num    
    % ppants x Num shuff. 
    %
    % all have 3 DV types (frequency, duartion, mduration).
    %%
    TableforExport=table();
    thisentry=1; %initialize counter.
    %
    for itype = 1:3
        switch itype
            case 1
                %first case is Hz data
                searchfor ='*FreqPFI*';
            case 2
                searchfor ='*LocPFI*';
            case 3
                searchfor ='*NumPFI_acrossTrials*';            
                
        end
        %index relevant variables in workspace
        varsstore = whos([searchfor], 'variables');
                
                %cycle through and add to table (must have same number of
                %rows (ie ppants).
                for indvar= 1:length(varsstore)
                    tmpdata = eval(varsstore(indvar).name);
                    mtmp = squeeze(nanmean(tmpdata,2));
                    if size(mtmp,1)~=length(allppants);
                        mtmp= mtmp(allppants,:);
                    end
                    %append horizontally to table.
                    TableforExport=[TableforExport, table(mtmp)];
                    TableforExport.Properties.VariableNames(thisentry)={varsstore(indvar).name};
                    thisentry=thisentry+1;
                end
    
    %
    end
    %%
    try cd('/Users/matthewdavidson/Desktop')
    catch
        cd('/Users/MattDavidson/Desktop')
    end
        
    cd('stats for JASP')
    writetable(TableforExport, 'BEHAVIOURAL_PFI.csv')
    
end
if job.exportDatatoExcelfolonDesktop_jamovi ==1
    %% also for JAMOVI - need to code separately and carefully for LME analyses.
    cd(basedir)
    %load relevant data
    load('PFI_data_concat')        
    load('Shuffleddata')
    load('MD_AllBP_perppant.mat')
    
    
    outputtoDESKTOPorSAVEinmatlab=2;
    
    ppantID=[1:length(allppants)];
    
    %cycle through and add to table (must have same number of
    %rows (ie ppants).
    %
    for itype = 3%:4
        switch itype
            %                         case 1
            %                             %first case is Hz data
            %                             searchfor ='*FreqPFI*';
            %                             nameis= 'Hz';
            %                         case 2
            %                             searchfor ='*_LocPFI*';
            %                             nameis= 'Loc';
            case 3
                searchfor ='*NumPFI_acrossTrials*';
                nameis= 'nPFI';
                
            case 4
                searchfor ='*HzxLocPFI*';
                nameis= 'HzxLoc';
        end
        
        TableforExport=table();
        %index relevant variables in workspace
        varsstore = whos([searchfor], 'variables');
        for indvar= 1:length(varsstore)
            
            %For each DV type we need a new table.
            if itype~=3  %then we t need a new table for each DV
                TableforExport=table();
            end
            %
            tmpdata = eval(varsstore(indvar).name);
            %mean over trials.
            mtmp = squeeze(nanmean(tmpdata,2));
            if size(mtmp,1)~=length(allppants);
                
                if itype<4
                    mtmp= mtmp(allppants,:);
                else
                    mtmp= mtmp(allppants,:,:);
                end
            end
            
            
            %create indexes for table entry
            if itype==3
                nreps = size(mtmp,2);
                pIDs = repmat(ppantID', [nreps, 1]);

                %add each levels of this variable
                for ilevel = 1:nreps
                    
                    if mod(indvar,2)==0
                        
                        realrows = ppantID+ ((ilevel-1)*length(allppants))+88;
                    else
                        realrows = ppantID+ ((ilevel-1)*length(allppants));
                    end
                    
                    
                    TableforExport(realrows,1)=table(ppantID');
                
                    TableforExport(realrows,2)=table(repmat(ilevel,[length(realrows),1]));
                    %real or shuffled data?
                    if strcmp(varsstore(indvar).name(end-7:end), 'shuffled')
                        % if shuffled.
                        TableforExport(realrows,3)=table(2);
                    else
                        TableforExport(realrows,3)=table(1);
                    end
                    TableforExport(realrows,4)=table(mtmp(:,ilevel));
                end
                TableforExport.Properties.VariableNames(1)={'subjectID'};
                TableforExport.Properties.VariableNames(2)={nameis};
                TableforExport.Properties.VariableNames(3)={'realvsshuff'};
                TableforExport.Properties.VariableNames(4)={'DV'};
            else
                
                
                %we have an extra dimension in the hzxLoc case
                nreps = size(mtmp,2)*size(mtmp,3);
                
                pIDs = repmat(ppantID', [nreps, 1]);
                TableforExport(1:length(pIDs),1)=table(pIDs);
                TableforExport.Properties.VariableNames(1)={'subjectID'};
                %add each levels of this variable
                
                for ilevel1 = 1:size(mtmp,2) %for first dimension
                    
                    for ilevel2 = 1:size(mtmp,3) %for first dimension
                        
                        realrows = ppantID+ ((ilevel2-1)*length(allppants))+ ((ilevel1-1)*length(allppants)*size(mtmp,3));
                        %include indexing for later analysis:
                        TableforExport(realrows,2)=table(repmat(ilevel1,[length(realrows),1]));
                        TableforExport(realrows,3)=table(repmat(ilevel2,[length(realrows),1]));
                        %and store values
                        TableforExport(realrows,4)=table(mtmp(:,ilevel1,ilevel2));
                    end
                end
                TableforExport.Properties.VariableNames(2)={'Hz'};
                TableforExport.Properties.VariableNames(3)={'Loc'};
                TableforExport.Properties.VariableNames(4)={'DV'};
                
                
                
                
                
                
            end
            
            
            % only save if finished.
            if itype==3 && mod(indvar,2)==0 || itype==4 
                
                if outputtoDESKTOPorSAVEinmatlab==1
                    %output to matlab as friendly file for other
                    %packages.
                    try cd('/Users/matthewdavidson/Desktop')
                    catch
                        cd('/Users/MattDavidson/Desktop')
                    end
                    
                    cd('stats for JAMOVI')
                    
                    switch itype
                       case 3
                           
                           if indvar==2 %save after appending shuffled data to observed, per DV type
                               
                               
                               
                                   writetable(TableforExport, ['FreqPFI_nPFIwShuff_tbl.csv']) 
                               
                           elseif indvar ==4                               
                               writetable(TableforExport, ['mDurPFI_nPFIwShuff_tbl.csv']) 
                           elseif indvar==6
                               
                              writetable(TableforExport, ['tDurPFI_nPFIwShuff_tbl.csv']) 
                           end
                                    
                       case 4
                           
                            if indvar==1                                
                                writetable(TableforExport, ['FreqPFI_HzxLoc_tbl.csv'])
                            elseif indvar== 2
                                                            
                                writetable(TableforExport, ['mDurPFI_HzxLoc_tbl.csv'])
                            elseif indvar== 3
                            
                               writetable(TableforExport, ['tDurPFI_HzxLoc_tbl.csv'])                             
                            end
                    end
                    
                else %keep working in matlab.
                   switch itype
                       case 3
                           
                           if indvar==2 %save after appending shuffled data to observed, per DV type
                               
                               FreqPFI_nPFIwShuff_tbl=TableforExport;
                               try save('BEHAVIOURAL_PFI_forLME', 'FreqPFI_nPFIwShuff_tbl', '-append')
                               catch
                                   save('BEHAVIOURAL_PFI_forLME', 'FreqPFI_nPFIwShuff_tbl')
                               end
                           elseif indvar ==4
                               mDurPFI_nPFIwShuff_tbl=TableforExport;
                               save('BEHAVIOURAL_PFI_forLME', 'mDurPFI_nPFIwShuff_tbl', '-append')
                               
                           elseif indvar==6
                               tDurPFI_nPFIwShuff_tbl=TableforExport;
                               save('BEHAVIOURAL_PFI_forLME', 'tDurPFI_nPFIwShuff_tbl', '-append')
                           end
                                    
                       case 4
                           
                            if indvar==1
                                FreqPFI_HzxLoc_tbl=TableforExport;
                                save('BEHAVIOURAL_PFI_forLME', 'FreqPFI_HzxLoc_tbl', '-append')

                            elseif indvar== 2
                            
                                mDurPFI_HzxLoc_tbl=TableforExport;
                                save('BEHAVIOURAL_PFI_forLME', 'mDurPFI_HzxLoc_tbl', '-append')
                            elseif indvar== 3
                            
                                tDurPFI_HzxLoc_tbl=TableforExport;
                                save('BEHAVIOURAL_PFI_forLME', 'tDurPFI_HzxLoc_tbl', '-append')
                                                            
                            end
                    end
                    
                end
            end
        end
        
    end
end
           
    
    
    if job.LMEonHzxLOC==1
         cd(basedir)
         
         
         load('BEHAVIOURAL_PFI_forLME.mat');
    tresults=table();
    icount=1;
         for id=1:3
             switch id
                 case 1
                     datatable=FreqPFI_HzxLoc_tbl;
                     DVtype='FreqPFI';                     
                 case 2
                     datatable=mDurPFI_HzxLoc_tbl;
                     DVtype='mDurPFI';
                 case 3
                     datatable=tDurPFI_HzxLoc_tbl;                     
                     DVtype='toralDurPFI';
             end
             
             for idrop = 1:3 % compare interaction, hz and loc.
                 switch idrop
                     case 1 %test interaction.
                         UNterm1='DV ~  Hz*Loc + (1|subjectID)';
                         UNterm2='DV ~  Hz*Loc - Hz:Loc + (1|subjectID)';
                         dropis= 'FE Int';
                     case 2 %test hz
                         UNterm1='DV ~  Hz*Loc + (1|subjectID)';
                         UNterm2='DV ~  Hz*Loc - Hz +(1|subjectID)';
                         dropis= 'FE Hz';
                     case 3 %test Loc
                         UNterm1='DV ~  Hz*Loc + (1|subjectID)';
                         UNterm2='DV ~  Hz*Loc-Loc+ (1|subjectID)';
                         dropis= 'FE Loc';
                 end
                 
                 %make sure to change the correct cols to nominal.
                datatable.Hz = nominal(datatable.Hz);
                datatable.Loc =nominal(datatable.Loc);
                datatable.subjectID=nominal(datatable.subjectID);
                
    %creating models:
    lmeUN =fitlme(datatable, UNterm1);
    lmeUN2 =fitlme(datatable, UNterm2);
    
    tmp=compare(lmeUN2, lmeUN);
    
    %store chi-sq; 
    tresults(icount,1) = table({DVtype});
    tresults(icount,2) = table({dropis});
    tresults(icount,3) = table(tmp.LRStat(2));
    tresults(icount,4) = table(tmp.deltaDF(2));
    tresults(icount,5) = table(tmp.pValue(2));
    
    
                icount=icount+1;
                
                
    %alternate ways of calculating:  
%     [h,pValue,stat] =lratiotest(lmeUN.LogLikelihood,lmeUN2.LogLikelihood,abs(lmeUN.DFE-lmeUN2.DFE))

%     [h,pValue,stat] =LogLikelihoodReport(lmeUN,lmeRe);

    
             end
         end
         tresults.Properties.VariableNames(1)= {'DV'};
         tresults.Properties.VariableNames(2)= {'FixedEffect'};
         tresults.Properties.VariableNames(3)= {'teststat'};
         tresults.Properties.VariableNames(4)= {'DoF'};
         tresults.Properties.VariableNames(5)= {'pval'};
         %display results
         disp(tresults)
        
         Finalresults_HzxLoc=tresults;
         save('BEHAVIOURAL_PFI_forLME.mat', 'Finalresults_HzxLoc', '-append')
    end
    
    
    if job.LMEnPFIwShuff==1
        cd(basedir)
        
        addpath('/Users/MattDavidson/Documents/MATLAB/Toolbox/LME resources')
        load('BEHAVIOURAL_PFI_forLME.mat');
        %%
        tresults=table();
        % note to drop the '0' targets PFI, restrict to these rows:
        nPFIonly=1;
        pfirows= [23:88, 111:176];
        
        icount=1;
        for idatatype=1:2
            if idatatype==1
                if nPFIonly==1
                    useindx = 23:88;
                else
                useindx=1:88; %real data in table.
                end
                datatypewas='Observed';
            else
                if nPFIonly==1
                    useindx= 111:176;
                else
                useindx=89:176; %shuffled data in table.
                end
                datatypewas='Shuffled';
            end
            
        for id=1:3
            switch id
                case 1
                    datatable=FreqPFI_nPFIwShuff_tbl(useindx,:);
                    DVtype='PFI per trial';
                case 2
                    datatable=mDurPFI_nPFIwShuff_tbl(useindx,:);
                    DVtype='PFI duration';
                case 3
                    datatable=tDurPFI_nPFIwShuff_tbl(useindx,:);
                    DVtype='Total duration';
            end
            
            
            
            
            for idrop = 4%1:4 % compare interaction, hz and loc.
                switch idrop
                    case 1 %test interaction.
                        UNterm1='DV ~ realvsshuff*nPFI + (1|subjectID)';
                        UNterm2='DV ~ realvsshuff*nPFI -realvsshuff:nPFI+ (1|subjectID)';
                        dropis= 'FE Int';
                    case 2 %test nPFI
                        UNterm1='DV ~   realvsshuff*nPFI +  (1|subjectID)';
                        UNterm2='DV ~  realvsshuff*nPFI - nPFI + (1|subjectID)';
                        dropis= 'FE nPFI';
                    case 3 %test real vs shuff
                        UNterm1='DV ~   realvsshuff*nPFI + (1|subjectID)';
                        UNterm2='DV ~  realvsshuff*nPFI - realvsshuff +(1|subjectID)';
                        dropis= 'FE realvsshuff';
                    case 4 %test intercepts
                        UNterm1='DV ~   nPFI + (1|subjectID)';
                        UNterm2='DV ~  1+(1|subjectID)';
                        dropis= 'FE nPFI.';
                        
                end
                
                %%
                %make sure to change the correct cols to nominal.
                datatable.subjectID = nominal(datatable.subjectID);
                datatable.realvsshuff =categorical(datatable.realvsshuff);
                datatable.nPFI =ordinal(datatable.nPFI);
                
            
                %creating models:
                lmeUN =fitlme(datatable, UNterm1);%, 'fitMethod', 'REML');
                lmeRe =fitlme(datatable, UNterm2);%,'fitMethod', 'REML';
                
                tmp=compare(lmeRe, lmeUN);
                
%                     [h,pValue,stat] =LogLikelihoodReport(lmeUN,lmeRe);
                
                %%
                %store chi-sq;
                
                tresults(icount,1) = table({DVtype});
                tresults(icount,2) = table({dropis});
                tresults(icount,3) = table(tmp.LRStat(2));
                tresults(icount,4) = table(tmp.deltaDF(2));
                tresults(icount,5) = table(tmp.pValue(2));
                tresults(icount,6) = table({datatypewas});
                
                icount=icount+1;
                
                
                %alternate ways of calculating:
                %     [h,pValue,stat] =lratiotest(lmeUN.LogLikelihood,lmeUN2.LogLikelihood,abs(lmeUN.DFE-lmeUN2.DFE))
                
                
                
            end
        end
        
        
    end
        tresults.Properties.VariableNames(1)= {'DV'};
        tresults.Properties.VariableNames(2)= {'Factor'};
        tresults.Properties.VariableNames(3)= {'teststat'};
        tresults.Properties.VariableNames(4)= {'DoF'};
        tresults.Properties.VariableNames(5)= {'pval'};
        tresults.Properties.VariableNames(6)= {'Modeled'};
        %display results
        disp(tresults)
        Finalresults_nPFIxrealvsshuff=tresults;
        %%
         save('BEHAVIOURAL_PFI_forLME.mat', 'Finalresults_HzxLoc', '-append')
        
    end
    
    if job.plotlikelihoodofPFIovertime==1;
       %%
        cd(basedir)
        load('PFI_data')
       %%
       ppantbytrialdur=zeros(length(allppants), 3600);
       
       icounter=1;
       for ippant= allppants;
           
           PFIbytrialdur=[];
           
           for itrial=1:24
               
               allPFI = sum(PFI_only(ippant).Trial(itrial).allBPs,1)  ;
               PFIbytrialdur=[PFIbytrialdur; allPFI];
               %correct to proportion
%                PFIbytrialdur=mean(PFIbytrialdur,1);
               %sanitycheck
              
               
           end 
          ppantbytrialdur(icounter,:)=nanmean(PFIbytrialdur,1); 
          icounter=icounter+1;
       end
       %%
       figure(1);
       
      set(gcf, 'color', 'w')      
       plot([1:3600]./60, nanmean(ppantbytrialdur,1), 'linew', 2)
       ylabel('mean buttons pressed')
       xlabel('Trial time')
       set(gca, 'fontsize', 25)
       axis tight
       %%
      
        
    end