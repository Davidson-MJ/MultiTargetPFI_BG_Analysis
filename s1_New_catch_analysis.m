
%script 1
%Reshape B'l data into more managable chunks.
% MD Jan 2017
clear all
dbstop if error
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/mDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
end 
basefol=pwd;
load('Raw_Behavioral_AllPP.mat')
%%
cd('newplots-MD')
%%
load('MD_AllBP_perppant.mat')
ppantDetails=[];%ppantTrialDatawithDetails;
% cd ../
%
fontsize=15;
%note, current analysis removes ppants without PFI.
PFIppants=[1:3,7:9, 11:29];

%after rejection based on individual catch analysis:
allppants=1:29;
badppants = [8, 15,28, 4 , 7,5 6]; % catch fail
allppants(badppants)=[];


job.updateCatchdetails=0; %
job.epocharoundCatchonset=0;
job.epocharoundCatchoffset=0;
job.calcSignificanceoftraceafterOnset=0;

job.sortCatchtracesbynumTargets=0; %onset and offsets. appended to 'catch performance'
%$

job.calcNEWShuffledCatchOnsetandOffsetBPProb=0; %new ver MD # 04/18.
job.calcNEWShuffledCatchOFFSETBPProb=0; %new ver MD # 04/18. % using not catch offset times, but times when BPs are pressed.


%plotting
job.plotPpantProbandSig=0; %imagesc (OLD output). Ranks participants.
%
job.plotCatchtracesbynumTargets_eachppant=0;  % Ppant responses. also plots overall RTs to catch onset/offset.
job.plotCatchtracesbynumTargets=0; % not really by num targets, plots vs shuffled likelihood.
%
job.plotmissedcatchesbynumremoved=0;
job.Plotmissedcatches_bylocationperppant=0;

%catch output
job.calcMissedcatchesbyTargperppant=1;
job.calcmissedcatches_bylocationperppant=1;




%%
% goodppants = [1:7, 9:14, 16:27,29];
%set up xaxis based on epoch length.
%start with a larger window, but trim the output when plotting.
%note that the 'catch start' time recorded, actually begins the ramp, so
%total target absence is catchstart+1.5seconds.
window = [-4 4].*60; %8 second window.
onset=abs(window(1));
%xaxis spacing for plots using this size window.
xtONSET=[window(1):1:window(2)]/60;
xtOFFSET=[-1*(window(2)):1:abs(window(1))]/60;
% ppantDetails=[];
%
% %ppants 5:24 actually have those targets with a "0" disappear instead of a
% %"1", so change to accomodate.
% %starting at 4th ppant, up to 24th?


% error('check window (changed to -2:8s) for plotting below.')

if job.updateCatchdetails==1
    for i=(24*4)+1:(24*24)
        % adjust bottom Left
        if sum(AllTrialsAllPPBottomLeft(i,6:9))>0 %if it is a trial with atleast one recorded catch, switch the logical index
            for icol=6:9
                if AllTrialsAllPPBottomLeft(i, icol)==0
                    AllTrialsAllPPBottomLeft(i, icol)=1; %actually a catch
                else %is 1
                    AllTrialsAllPPBottomLeft(i, icol)=0; %not a catch
                end
                
                
                if AllTrialsAllPPTopLeft(i, icol)==0
                    AllTrialsAllPPTopLeft(i, icol)=1;
                else
                    AllTrialsAllPPTopLeft(i, icol)=0;
                end
                
                
                if AllTrialsAllPPBottomRight(i, icol)==0
                    AllTrialsAllPPBottomRight(i, icol)=1;
                else
                    AllTrialsAllPPBottomRight(i, icol)=0;
                end
                
                if AllTrialsAllPPTopRight(i, icol)==0
                    AllTrialsAllPPTopRight(i, icol)=1;
                else
                    AllTrialsAllPPTopRight(i, icol)=0;
                end
            end
        end
    end
    %%
    
    for ippant = 1:29
        trials = (1:24)+(24*(ippant-1));
        ppantDATA = zeros(24,4,3600) ;%trial, loc, frames
        %concat data per BP location.
        ppantDATA(:, 1,:) = AllTrialsAllPPTopLeft(trials, 15:3614);
        ppantDATA(:, 2,:) = AllTrialsAllPPTopRight(trials, 15:3614); % Check !
        ppantDATA(:, 3,:) = AllTrialsAllPPBottomLeft(trials, 15:3614); % check!
        ppantDATA(:, 4,:) = AllTrialsAllPPBottomRight(trials, 15:3614);
        
        for itrial=1:24
            %
            
            ppantDetails(ippant).TrialDetails(itrial).TL_Catch= AllTrialsAllPPTopLeft(trials(itrial),6);
            ppantDetails(ippant).TrialDetails(itrial).TR_Catch= AllTrialsAllPPTopRight(trials(itrial),7);
            ppantDetails(ippant).TrialDetails(itrial).BL_Catch= AllTrialsAllPPBottomLeft(trials(itrial),8);
            ppantDetails(ippant).TrialDetails(itrial).BR_Catch= AllTrialsAllPPBottomRight(trials(itrial),9);
            
            ppantDetails(ippant).TrialDetails(itrial).TL_Freq= AllTrialsAllPPTopLeft(trials(itrial),2);
            ppantDetails(ippant).TrialDetails(itrial).TR_Freq= AllTrialsAllPPTopRight(trials(itrial),3);
            ppantDetails(ippant).TrialDetails(itrial).BL_Freq= AllTrialsAllPPBottomLeft(trials(itrial),4);
            ppantDetails(ippant).TrialDetails(itrial).BR_Freq= AllTrialsAllPPBottomRight(trials(itrial),5);
            
            %timing is same for all, or should be.
            %adjusted to account for NO RAMP at beginning, and catch end =
            %end ramp.
            catch_startreal =  AllTrialsAllPPTopLeft(trials(itrial),13)+(1.5*60); %added ramp duration (completely gone)
            catch_endreal=  AllTrialsAllPPTopLeft(trials(itrial),14)- 1.5*60 ; %end at ramp offset begin (no longer completely gone).
            
            catchdur = catch_endreal - catch_startreal;
            ppantDetails(ippant).TrialDetails(itrial).Catchstart_frames=catch_startreal;
            ppantDetails(ippant).TrialDetails(itrial).Catchend_frames= catch_endreal;
            ppantDetails(ippant).TrialDetails(itrial).totalCatchdur=catchdur ;
            ppantDetails(ippant).TrialDetails(itrial).TotalCatchTargetsRemoved= sum(AllTrialsAllPPTopLeft(trials(itrial), 6:9));
            
            ppantDetails(ippant).TrialDetails(itrial).BPData = squeeze(ppantDATA(itrial,:,:));
        end
        ppantDetails(ippant).AllBPData=ppantDATA;
        
    end
    %%
    ppantTrialDatawithDetails=ppantDetails;
    save('MD_AllBP_perppant', 'ppantTrialDatawithDetails')
    clearvars ppantDetails ppantDATA
end
%create catch Epoch for BPs
if job.epocharoundCatchonset==1
    %%NOTE THAT THERE WAS NO RAMP!!, at the end of the disap ramp the catch abruptly was removed, and
    %returned at the end of the ramp.
    for ippant = 1:29
        for itrial=1:24
            Details= ppantTrialDatawithDetails(ippant).TrialDetails(itrial);
            tmpData= Details.BPData;
            Catch_any=[];
            catchstartTrial =Details.Catchstart_frames;
            if Details.TotalCatchTargetsRemoved>0
                %collect catch BP if it exists.
                
                if Details.TL_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(1,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                if Details.TR_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(2,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                if Details.BL_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(3,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                if Details.BR_Catch==1 %epoch around catch
                    Catch_any= [Catch_any; tmpData(4,catchstartTrial-(abs(window(1))):catchstartTrial+window(2))];
                end
                
                %now we have epoched -2:6s around 'catch start',
                %collect all BPs (even if there was no catch) during the longer
                %catch epoch.
                %             tmpAll = tmpData(:, catchstartTrial-(abs(window(1))):catchstartTrial+window(2));
                % ntrain_all= zeros(4,size(Catch_any,2));%
                % ntrain_all(:,catchstartEpoch:Details.totalCatchdur+catchstartEpoch)= tmpAll(:,catchstartTrial:Details.totalCatchdur+catchstartTrial);
                
                %save the long form Epochs, for catch trace BPs
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPs_longepoch = Catch_any; %may not contain 4 rows, and different catch lengths.
                
                %also save JUST the BPs during physical removal (removing pre and post windows), to analyse exact
                %catch accuracy (duration of physical absence recorded).
                ntrain_slim= zeros(size(Catch_any));
                
                physicalabsence = Details.totalCatchdur; %exclude off and on ramp, as not completely absent target.
                start_physicalabsence = abs(window(1)) ; %start in frames
                %record button press for the catch period, excluding the pre
                %catch epoched window,
                %ramp.
                catchlocs = find([Details.TL_Catch, Details.TR_Catch, Details.BL_Catch, Details.BR_Catch]) ;
                catchexactonly = tmpData(catchlocs,catchstartTrial:catchstartTrial+Details.totalCatchdur);
                
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPs_exact= catchexactonly;
                %             ppantTrialDatawithDetails(ippant).TrialDetails(itrial).AllBPs_duringcatchexact= ntrain_all;
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).truephysicalabsence_duration= physicalabsence;
            end
        end
    end
    %%
    
    save('MD_AllBP_perppant', 'ppantTrialDatawithDetails', '-append')
    %% now stack all to see if it matches Irene's plot.
    stackALL=[];
    %%
    for ippant = 1:29
        catchProb=[];
        
        for itrial=1:24
            
            %%skip trials that required 4 button presses.
            tmpCatchBPs=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPs_longepoch;
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TotalCatchTargetsRemoved<4
                
                %only assess catch timing/Prob, if button wasn't already
                %pressed at catch onset.
                
                onset=abs(window(1));

%                                 keeptrials=[];
%                 for iBP=1:size(tmpCatchBPs,1)
%                     %for each trial, only retain for analysis of prob if BP
%                     %not existing
%                     tmpt=tmpCatchBPs(iBP,:);
%                     if sum(tmpt(1,onset-12:onset))<1 %12 frames is .2 seconds.
%                         keeptrials=[keeptrials; tmpt];
%                     end
%                     
%                 end
                keeptrials = tmpCatchBPs;                
                catchProb=[catchProb; mean(keeptrials,1)]; %average
            end
            
        end
        %for save
        ppantTrialDatawithDetails(ippant).eachTrialCatchProb_combOnetoThreetargets = catchProb;
        ppantTrialDatawithDetails(ippant).mTrialCatchProb_combOnetoThreetargets = mean(catchProb,1);
        
        %for plot/calc of pvalues
        
        % average every 250ms, to reduce comp time and redundant info.
        % 250ms = 60frames/4 = 15 frames.
        newvec = [1:15:length(Catch_any)];
        outgoing=zeros(size(catchProb,1), length(newvec));
        %for each catch trace
        for itrace= 1:size(catchProb,1)
            tmptr=squeeze(catchProb(itrace,:));
            %average within bins.
            for ibin = 1:length(newvec)-1
                outgoing(itrace,ibin) = mean(tmptr(newvec(ibin):newvec(ibin+1)));
            end
        end
        %separate trials for significance calculations
        stackALL(ippant).data= catchProb;
        stackALL_binned(ippant).data= outgoing;
        %mean for plots, per ppant.
        stackALLm(ippant,:)= mean(catchProb,1);
        stackALLm_binned(ippant,:)= mean(outgoing,1);
    end
    catchBP_ProbtraceperPpant_onetothreetargets = stackALL;
    catchBP_ProbtraceperPpant_onetothreetargets_binned = stackALL_binned;
    catchBP_ProbtraceperPpant_onetothreetargets_mean=stackALLm;
    catchBP_ProbtraceperPpant_onetothreetargets_meanbinned=stackALLm_binned;
    
    save('Catchperformance', 'catchBP_ProbtraceperPpant_onetothreetargets', ...
        'catchBP_ProbtraceperPpant_onetothreetargets_binned',...
        'catchBP_ProbtraceperPpant_onetothreetargets_mean',...
        'catchBP_ProbtraceperPpant_onetothreetargets_meanbinned');
end
%%
if job.epocharoundCatchoffset==1
    %%NOTE THAT THERE WAS NO RAMP!!, at the end of the disap ramp the catch abruptly was removed, and
    %returned at the end of the ramp.
    for ippant = 1:29
        for itrial=1:24
            Details= ppantTrialDatawithDetails(ippant).TrialDetails(itrial);
            tmpData= Details.BPData;
            Catch_any=[];
            %now catch end frames!
            catchstartTrial =Details.Catchend_frames;
            
            %since we are interested in 5 before, and 1 second after,
            %change window below.
            
            if Details.TotalCatchTargetsRemoved>0
                %collect catch BP if it exists.
                
                
                catchlocs = find([Details.TL_Catch, Details.TR_Catch, Details.BL_Catch, Details.BR_Catch]) ;
                
                
                if catchstartTrial+abs(window(2)) > length(tmpData)% too long to epoch.
                    nantrain = nan(length(catchlocs), length(xtOFFSET));
                %whats the difference?
                tmpdata = tmpData(catchlocs, catchstartTrial-abs(window(1)):end);
                nantrain(:,1:length(tmpdata)) = tmpdata;
                Catch_any = nantrain;                   
                else
                
                Catch_any = tmpData(catchlocs, catchstartTrial-abs(window(1)):catchstartTrial+abs(window(1)));
                end
                
                %now we have epoched -2:6s around 'catch start',
                %collect all BPs (even if there was no catch) during the longer
                %catch epoch.
                %             tmpAll = tmpData(:, catchstartTrial-(abs(window(1))):catchstartTrial+window(2));
                % ntrain_all= zeros(4,size(Catch_any,2));%
                % ntrain_all(:,catchstartEpoch:Details.totalCatchdur+catchstartEpoch)= tmpAll(:,catchstartTrial:Details.totalCatchdur+catchstartTrial);
                
                %save the long form Epochs, for catch trace BPs
                ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPsOFFSET_longepoch = Catch_any; %may not contain 4 rows, and different catch lengths.
                
                %             %also save JUST the BPs during physical removal (removing pre and post windows), to analyse exact
                %             %catch accuracy (duration of physical absence recorded).
                %             ntrain_slim= zeros(size(Catch_any));
                %
                %             physicalabsence = Details.totalCatchdur; %exclude off and on ramp, as not completely absent target.
                %             start_physicalabsence = abs(window(1)) ; %start in frames
                %             %record button press for the catch period, excluding the pre
                %             %catch epoched window,
                %             %ramp.
                %             catchexactonly = Catch_any(:,[start_physicalabsence:start_physicalabsence+ physicalabsence]);
                %
                %             ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPs_exact= catchexactonly;
                %             %             ppantTrialDatawithDetails(ippant).TrialDetails(itrial).AllBPs_duringcatchexact= ntrain_all;
                %             ppantTrialDatawithDetails(ippant).TrialDetails(itrial).truephysicalabsence_duration= physicalabsence;
            end
        end
    end
    %%
    
    save('MD_AllBP_perppant', 'ppantTrialDatawithDetails', '-append')
    %% now stack all to see if it matches Irene's plot.
    stackALL=[];
    %%
    for ippant = 1:29
        catchProb=[];
        
        for itrial=1:24
            
            %%skip trials that required 4 button presses.
            tmpCatchBPs=ppantTrialDatawithDetails(ippant).TrialDetails(itrial).CatchBPsOFFSET_longepoch;
            if ppantTrialDatawithDetails(ippant).TrialDetails(itrial).TotalCatchTargetsRemoved<4
                
                
                
                %%%% below not necessary for catch OFFSET
                %only assess catch timing/Prob, if button wasn't already
                %pressed at catch onset.
                
                %                 onset=abs(window(1));
                %                 keeptrials=[];
                %                 for iBP=1:size(tmpCatchBPs,1)
                %                     %for each trial, only retain for analysis of prob if BP
                %                     %not existing
                %                     tmpt=tmpCatchBPs(iBP,:);
                %                     if sum(tmpt(1,onset-12:onset))<1 %12 frames is .2 seconds.
                %                         keeptrials=[keeptrials; tmpt];
                %                     end
                %
                %                 end
                
                catchProb=[catchProb; mean(tmpCatchBPs,1)]; %average
            end
            
        end
        %for save
        ppantTrialDatawithDetails(ippant).eachTrialCatchProbOFFSET_combOnetoThreetargets= catchProb;
        ppantTrialDatawithDetails(ippant).mTrialCatchProbOFFSET_combOnetoThreetargets = mean(catchProb,1);
        
        %for plot/calc of pvalues
        
        % average every 250ms, to reduce comp time and redundant info.
        % 250ms = 60frames/4 = 15 frames.
        newvec = [1:15:length(Catch_any)];
        outgoing=zeros(size(catchProb,1), length(newvec));
        %for each catch trace
        for itrace= 1:size(catchProb,1)
            tmptr=squeeze(catchProb(itrace,:));
            %average within bins.
            for ibin = 1:length(newvec)-1
                outgoing(itrace,ibin) = mean(tmptr(newvec(ibin):newvec(ibin+1)));
            end
        end
        %separate trials for significance calculations
        stackALL(ippant).data= catchProb;
        stackALL_binned(ippant).data= outgoing;
        %mean for plots, per ppant.
        stackALLm(ippant,:)= mean(catchProb,1);
        stackALLm_binned(ippant,:)= mean(outgoing,1);
    end
    catchBP_ProbtraceperPpantOFFSET_onetothreetargets = stackALL;
    catchBP_ProbtraceperPpantOFFSET_onetothreetargets_binned = stackALL_binned;
    catchBP_ProbtraceperPpantOFFSET_onetothreetargets_mean=stackALLm;
    catchBP_ProbtraceperPpantOFFSET_onetothreetargets_meanbinned=stackALLm_binned;
    
    save('Catchperformance', 'catchBP_ProbtraceperPpantOFFSET_onetothreetargets', ...
        'catchBP_ProbtraceperPpantOFFSET_onetothreetargets_binned',...
        'catchBP_ProbtraceperPpantOFFSET_onetothreetargets_mean',...
        'catchBP_ProbtraceperPpantOFFSET_onetothreetargets_meanbinned', '-append');
end
%%


if job.calcSignificanceoftraceafterOnset==1
    %% Calculate significance of this probability trace, correct for multiple comparisons
    load('Catchperformance')
    ppantCatchsigvalues=[];
%     [~,timezero]=min(abs(xt-0));
    for ippant=1:29
        tmp=catchBP_ProbtraceperPpant_onetothreetargets_binned(ippant).data;
        p=[];
        for itime=1:length(tmp)
            tdata= tmp(:,itime);
            %compare not to zero, but to the likelihood of BP at timepoint
            %zero.
            %             [~,tmpp]=ttest(tdata,tmp(:,timezero));
            [~,tmpp]=ttest(tdata,0);
            if isnan(tmpp)
                tmpp=1;
            end
            p(itime)=tmpp;
        end
        %for save
        ppantTrialDatawithDetails(ippant).mTrialCatchProb_combOnetoThreetargets_sig = log10(p);
        
        %% for plot
        ppantCatchsigvalues(ippant,:)=log10(p);
    end
    
    %%
    catchBP_ProbtraceperPpant_onetothreetargets_sig = ppantCatchsigvalues;
    save('Catchperformance', 'catchBP_ProbtraceperPpant_onetothreetargets_sig', '-append')
end
%% PLOT
if job.plotPpantProbandSig==1
    %     load('CatchPerformance')
    load('Catchperformance.mat')
    clf
    figure(1)
%       st=suptitle(['Catch Report Accuracy:']);
%     st.FontWeight = 'Bold';
%     st.FontSize = fontsize*2;
    colormap('jet')
    
    usebinned=0;
    if usebinned==1
        dataonsetis= catchBP_ProbtraceperPpant_onetothreetargets_meanbinned;
        dataoffsetis= catchBP_ProbtraceperPpantOFFSET_onetothreetargets_meanbinned;
    else
        dataonsetis = catchBP_ProbtraceperPpant_onetothreetargets_mean;
        dataoffsetis = catchBP_ProbtraceperPpantOFFSET_onetothreetargets_mean;
    end
    for ip = 1:2
        switch ip
            case 1
                datais= dataonsetis;
            case 2
                datais = dataoffsetis;
        end
        
        
        %keepdata for only ppants who experienced PFI
        catchBP_ProbtraceperPpant_onetothreetargets=squeeze(datais(PFIppants,:));
        
        %
        y=1:size(catchBP_ProbtraceperPpant_onetothreetargets,1);
        
        %sort according to mean prob
        tmp=mean(catchBP_ProbtraceperPpant_onetothreetargets,2);
        [srtp,ind]=sort(tmp,'descend');
        
        
        %plotting
        subplot(1,2,ip)
        if ip==2
            
            imagesc(xtOFFSET,y,catchBP_ProbtraceperPpant_onetothreetargets(ind,:))
            xlim([-2 2])
            eventis = 'return';
        else
            imagesc(xtONSET,y,catchBP_ProbtraceperPpant_onetothreetargets(ind,:))
            xlim([-2 2])
            eventis = 'removal';
        end
        %
        colb=colorbar;
        ylabel(colb,'P(Button Press)', 'fontsize', 1.5*fontsize)
        
        
        ylabel('Participant Rank', 'fontsize', 15)
        xlabel({['Time [sec] from physical target ' eventis]}, 'fontsize', 15)
        set(gca,'ytick', 1:29)
        set(gca,'fontsize', 20)
        hold on
        rmp=plot([0 0], ylim, ['w' ':'], 'linewidth', 5);
        
                title(['Target ' eventis], 'fontsize', fontsize*2)
    end
    %%
%     st=suptitle(['Catch Report Accuracy:']);
%     st.FontWeight = 'Bold';
%     st.FontSize = fontsize*2;
    cd(basefol)
    cd('newplots-MD')
    cd('Figures')
    cd('Catch analysis')
    %%
    print('-dpng', 'Catch TargetDisap and Reapp prob per ppant ranked')
    %%
    
    % %% no longer plotting significance this way.
    %
    %
    %     clf
    %     colormap('hot')
    %     colormap(flipud(colormap))
    %     catchBP_ProbtraceperPpant_onetothreetargets_sig= squeeze(catchBP_ProbtraceperPpant_onetothreetargets_sig(PFIppants,:,:));
    %
    %     y=1:size(catchBP_ProbtraceperPpant_onetothreetargets,1);
    %     imagesc(xt,y, catchBP_ProbtraceperPpant_onetothreetargets_sig(ind,:))
    %     % multiple correction
    %     %only interested in period immediately after catch
    %     % up to 3.5s
    %     [~,endlook] = min(abs(xt-3.5));
    %     q=fdr(squeeze(catchBP_ProbtraceperPpant_onetothreetargets_sig(:,onset:endlook)),.05);
    %     %%
    %     if q==0
    %         q=.05;
    %     end
    %     set(gca,'ytick', 1:29)
    %
    %     hold on;
    %     rmp=plot([0 0], ylim, ['w' ':'], 'linewidth', 5);
    %
    %     %
    %     c=colorbar;
    %
    %     caxis([min(min(catchBP_ProbtraceperPpant_onetothreetargets_sig)) log10(q)])
    %     caxis([log10(.001) log10(q)])
    %     %
    %     title({['Catch Report Accuracy:'];['Significantly increased likelihood of correct button press, after FDR correction']})
    %     ylabel('Participant Rank', 'fontsize', 15)
    %     xlabel({['Time [sec] from physical target removal'];['(Catch events)']}, 'fontsize', 15)
    %     ylabel(c,{['log_1_0 \itp \rmvalues']}, 'fontsize', 2*15)
    %     xlim([-1 3.5])
    %     set(gca,'fontsize', 20)
    %     shg
    %
    %     %%
    %     %
    %     print('-dpng', 'Catch TargetDisap significance per ppant-ranked');
    %
end
if job.sortCatchtracesbynumTargets==1
    cd(basefol)
    cd('newplots-MD')
    load('MD_AllBP_perppant.mat')
    
    for onsetoroffset = 1:2
        dataOUT=[];
        for ippant = 1:29
            ppantCatch1=[];
            ppantCatch2=[];
            ppantCatch3=[];
            ppantCatch4=[];
            
            for itrial=1:24
                
                Details= ppantTrialDatawithDetails(ippant).TrialDetails(itrial);
                
                Catch_any=[]; %reset variable per trial. how many to keep will change.
                
                numtargetsneeded=Details.TotalCatchTargetsRemoved;
                
                
                
                
                if onsetoroffset==1
                    Catch_any = Details.CatchBPs_longepoch; %10 second epoch of all required channels.
                else
                    Catch_any = Details.CatchBPsOFFSET_longepoch; %10 second epoch of all required channels.
                end
                
                %filter out those with/without button press at t=0.
                if onsetoroffset==1
                    onset=abs(window(1));
                    keeptrials=[];
                    for iBP=1:numtargetsneeded
                        %for each trial, only retain for analysis of prob if BP
                        %   not existing
                        tmpt=Catch_any(iBP,:);
%                         if sum(tmpt(1,onset-12:onset))<1 %if no buttons already pressed.
                            keeptrials=[keeptrials; tmpt];
%                         end
                        
                    end
                else
%                     keeptrials=Catch_any;
                    %at the moment no filtering.
                    onset=abs(window(2));
                    keeptrials=[];
                    for iBP=1:numtargetsneeded
%                         %for each trial, only retain for analysis of prob if BP
                        %   not existing
                        tmpt=Catch_any(iBP,:);
%                         if mean(mean(tmpt(:,onset-12:onset,1))) ==1 %12 frames is 0.2seconds
                            keeptrials=[keeptrials; tmpt];
%                         end
                        
                    end
                end
                %take mean over retained trials, having filtered.
                Catch_any=mean(keeptrials,1);
                
                switch numtargetsneeded
                    case 1
                        ppantCatch1=[ppantCatch1; Catch_any];
                    case 2
                        
                        ppantCatch2=[ppantCatch2; Catch_any];
                    case 3
                        ppantCatch3=[ppantCatch3; Catch_any];
                    case 4
                        ppantCatch4=[ppantCatch4; Catch_any];
                end
                
                
            end  
        
                
                
                %save per ppant.
                
                dataOUT(ippant).ProbTraceforOnetargetcatch=ppantCatch1;
                dataOUT(ippant).mProbTraceforOnetargetcatch=mean(ppantCatch1,1);
                
                dataOUT(ippant).ProbTraceforTwotargetcatch=ppantCatch2;
                dataOUT(ippant).mProbTraceforTwotargetcatch=mean(ppantCatch2,1);
                
                
                dataOUT(ippant).ProbTraceforThreetargetcatch=ppantCatch3;
                dataOUT(ippant).mProbTraceforThreetargetcatch=mean(ppantCatch3,1);
                
                dataOUT(ippant).ProbTraceforFourtargetcatch=ppantCatch4;
                dataOUT(ippant).mProbTraceforFourtargetcatch=mean(ppantCatch4,1);
                
           end     
                
                if onsetoroffset==1
                    catchBP_ProbtraceperPpant_by_numtargetsabsent = dataOUT;
                else
                    catchBP_ProbtraceperPpantOFFSET_by_numtargetsabsent = dataOUT;
                end
            
        end
    
    
    %%
    %also take average
    
     clf
    for onsetoroffset=1:2
        %reset outgoing variable
        stack=[];
        
        if onsetoroffset==1
            datais=catchBP_ProbtraceperPpant_by_numtargetsabsent;
            colis='g';
        else
            datais=catchBP_ProbtraceperPpantOFFSET_by_numtargetsabsent;
            colis='r';
        end
        
        for ippant=1:29
            %collect data for each participant, based on number of catches
            %removed.
            tmp=datais(ippant).mProbTraceforOnetargetcatch;
            try stack(ippant,1,:)=tmp;
            catch
            end
            
            
            tmp=datais(ippant).mProbTraceforTwotargetcatch;
            try stack(ippant,2,:)=tmp;
            catch
            end
            
            tmp=datais(ippant).mProbTraceforThreetargetcatch;
            try stack(ippant,3,:)=tmp;
            catch
            end
            
            %collect anyway, but will likely be ignored.
            
            tmp4=datais(ippant).mProbTraceforFourtargetcatch;
            if tmp4>0
                stack(ippant,4,:)=tmp4; %if we have data for 4 target catch disappearance.
            end
            
            
            
        end
        
        
        if onsetoroffset==1
            mean_catchBP_Probtraceperppant_bynum= stack;
            save('Catchperformance','mean_catchBP_Probtraceperppant_bynum', '-append')
        else
            mean_catchBP_ProbtraceperppantOFFSET_bynum= stack;
            save('Catchperformance','mean_catchBP_ProbtraceperppantOFFSET_bynum', '-append')
        end
        % long form widow was used to create traces.
        %%
    end
    
    
    cd(basefol)
    cd('newplots-MD')
    %%
    
    save('Catchperformance', 'catchBP_ProbtraceperPpant_by_numtargetsabsent',...
        'catchBP_ProbtraceperPpantOFFSET_by_numtargetsabsent',...
        '-append')
    
end
%%

if job.calcNEWShuffledCatchOnsetandOffsetBPProb==1
    %for each participant, calculate likelihood of BP around onset/offset of catch
    % in adjacent trials, to estimate RT.
    
    %%
    cd(basefol)
    cd('newplots-MD')
    load('MD_AllBP_perppant.mat')
    %%
    
    
    
    nshuff=200; 
    
    acrossall_mShuffledCatch_MD=zeros(length(ppantTrialDatawithDetails), 2, nshuff, length(xtOFFSET));
    for ippant = 1:length(ppantTrialDatawithDetails)
        % calc 200 shuffled trials SETs, each composed of a single
        % selection from any (1:24) of the trials. present trial included.
        
        pDetails =ppantTrialDatawithDetails(ippant);
        
        %List trial index by number of targets removed.
        numCatchppant=[];
        for itrial=1:24
            numCatchppant(itrial,:)=pDetails.TrialDetails(itrial).TotalCatchTargetsRemoved;
        end
        
        %for each type of catch trial (loc targets removed). extract a random
        %BP for the same time window, same location(s), any trial.
        
        ppantShuffled=nan(200,2,length(xtOFFSET)); %[nshuffled, disap/reap, nsamps]
       
        for ishuff=1:nshuff
        
        shuffdata=[];
        for itrial = 1:24
            %how many catch removed this trial?
            catchnum = numCatchppant(itrial);
            
            
            %find locations to compare
            locs = [pDetails.TrialDetails(itrial).TL_Catch, pDetails.TrialDetails(itrial).TR_Catch, pDetails.TrialDetails(itrial).BL_Catch, pDetails.TrialDetails(itrial).BR_Catch];
            locs=find(locs);
            
            
            %find trials to compare,
            availTrials = 1:24;
            %             availTrials= find(numCatchppant==catchnum);
            %             availTrials(availTrials==itrial)=[]; %remove the present trial
            
            %find real onset and offset of this catch
            catchtimes = [pDetails.TrialDetails(itrial).Catchstart_frames, pDetails.TrialDetails(itrial).Catchend_frames];
            
                % for each trial onset and offset, select the same time
                % window at random.
                for icatch=1:2 %for disap and reap
                    
                    tmp_time=catchtimes(icatch); %this is the catch start to collect about.
                    
                    
                    %select other trial
                    trialis=availTrials(randi(length(availTrials)));
                    BPis = squeeze(pDetails.AllBPData(trialis,locs,:));
                    %window about same catch time
                    if size(BPis,1)>size(BPis,2)
                        BPis=BPis';
                    end
                    
                    if (tmp_time + window(2))> size(BPis,2)
                              
                        nantrain = nan(size(BPis,1), length(xtOFFSET));
                            %whats the difference?
                        tmpdata = BPis(:, tmp_time-window(2):end);
                        nantrain(:,1:length(tmpdata)) = tmpdata;
                        storeBP= nantrain; 
                    else
                        storeBP = BPis(:, tmp_time-(window(2)):tmp_time+(window(2)));
                    end
                    
                    %store BP of this window
                    
                    shuffData(itrial,icatch,:) = nanmean(storeBP,1);
                    
                end
                
                end
                %                 if searchcount<200
                %store mean across trials for this participant, this
                %shuffle.
                ppantShuffled(ishuff,:,:)=squeeze(nanmean(shuffData,1));
                %                 end
                
                
        end
        %             mShuffledCatchProb = squeeze(nanmean(ppantShuffled,3)); %take mean over shuffles
        %             mShuffledCatchProb = squeeze(nanmean( ppantShuffled,1)); %or
        % take mean over trials.
        %will plot 99%CI of the shuffled mean across trials.
        acrossall_mShuffledCatch_MD(ippant,:,:,:,:) = permute(ppantShuffled,[2 1 3]);
    end
    dimsare = [{'ppants'}, {'Disap/Reap'},{'nshufftrials'},{'nsamps'}];
    cd(basefol)
    cd('newplots-MD')
    save('ShuffledCatch', 'acrossall_mShuffledCatch_MD', 'dimsare', '-append')
end

if job.calcNEWShuffledCatchOFFSETBPProb==1
    %for each participant, calculate likelihood of BP around onset/offset of catch
    % in adjacent trials, to estimate RT.
    
    %%
    cd(basefol)
    cd('newplots-MD')
    load('MD_AllBP_perppant.mat')
    %%
    %also load real catch performance, which we need to estimate the offset
    %shuffled likelihood (new MD#07/18)
    load('Catchperformance', 'mean_catchBP_ProbtraceperppantOFFSET_bynum');
    
    offsetdata = squeeze(nanmean(mean_catchBP_ProbtraceperppantOFFSET_bynum(:,1:3,:),2)); % just 1-3 targs.
    %note that time point zero in offset data is given by:
    offtzero = dsearchn(xtOFFSET', [0]);
    
   %now participant BP at zero = 
   ppantCatchOFFSETlikelihoodATZERO =offsetdata(:,offtzero);
            
   
   minClustlength = 30; % Fs= 60, this is the duration of BP to retain in our offset estimation.
    
    nshuff=200; 
    
    acrossall_mShuffledCatchOffset_MD=zeros(length(ppantTrialDatawithDetails), nshuff, length(xtOFFSET));
    for ippant = 1:length(ppantTrialDatawithDetails)
        % calc 200 shuffled trials SETs, each composed of a single
        % selection from any (1:24) of the trials. present trial included.
               
        pDetails =ppantTrialDatawithDetails(ippant);
        
        %List trial index by number of targets removed.
        numCatchppant=[];
        for itrial=1:24
            numCatchppant(itrial,:)=pDetails.TrialDetails(itrial).TotalCatchTargetsRemoved;
        end
        
        %for each type of catch trial (loc targets removed). extract a random
        %BP for the same time window, same location(s), any trial.
        
        ppantShuffled=nan(200,length(xtOFFSET)); %[nshuffled, disap/reap, nsamps]
       
        
        pCatchLHOOD = ppantCatchOFFSETlikelihoodATZERO(ippant);
        
        shuffdata=[];
        
        for ishuff=1:nshuff
        
        
        for itrial = 1:24
            %how many catch removed this trial?
            catchnum = numCatchppant(itrial);
            
            
            %find locations to compare
            locs = [pDetails.TrialDetails(itrial).TL_Catch, pDetails.TrialDetails(itrial).TR_Catch, pDetails.TrialDetails(itrial).BL_Catch, pDetails.TrialDetails(itrial).BR_Catch];
            locs=find(locs);
            
            
            %find trials to compare,
            availTrials = 1:24;
            %             availTrials= find(numCatchppant==catchnum);
            %             availTrials(availTrials==itrial)=[]; %remove the present trial
            
            
            
            %select other trial
            trialis=availTrials(randi(length(availTrials)));
            BPis = squeeze(pDetails.AllBPData(trialis,locs,:));
            
            %now find any period in BP, where the Pressed buttons equates
            %to the BP at catch offset, to estimate the release speed.
                   
            normBP = nansum(BPis,1)./size(BPis,1);
            %adapted to +- 10%
            if size(BPis,1)>size(BPis,2)
                BPis=BPis';
            end
            
            %find upper and lower bounds per ppant.
            timesBPressed1 = find((normBP > pCatchLHOOD-.2)) ;
            timesBPressed2= find((normBP < pCatchLHOOD+.2)) ;
            %determine clusters of appropriate 'bp'
            timesBPressed = find(ismember(timesBPressed2,timesBPressed1));
             
            if length(timesBPressed)>2
            %find clusters
              vect1 = diff(timesBPressed);
                v1 = (vect1(:)==1);
                d = diff(v1);
                
                
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                
                %any large clusters?
                cLength = clusterSTandEND(:,2)-clusterSTandEND(:,1);
                
                usecls = find(cLength>minClustlength);
                %for these clusters, take any random point within, and use
                %as centre for this shuffle.
                
                storeBP_tmp = nan(size(clusterSTandEND,1), length(xtOFFSET));
                for icl=usecls'
                    
                    
                    STC=timesBPressed(clusterSTandEND(icl,1));
                    ENDC=timesBPressed(clusterSTandEND(icl,2)+1);
                    
                    %select random point between.
                    tmp_time = randi([STC ENDC]);
                    
                    
                    %                 [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                     
                    try storeBP_tmp = BPis(:, tmp_time-(window(2)):tmp_time+abs(window(1)));
                    catch
                    end
                end
                
            else storeBP_tmp= nan(1,length(xtOFFSET));
            end
                
            
                
                
            
            
            %store BP of this window
            
            shuffdata(itrial,:) = nanmean(storeBP_tmp,1);
            end
        
                
                
                %                 if searchcount<200
                %store mean across trials for this participant, this
                %shuffle.
                ppantShuffled(ishuff,:)=squeeze(nanmean(shuffdata,1));
                                
                
                
        end
        %             mShuffledCatchProb = squeeze(nanmean(ppantShuffled,3)); %take mean over shuffles
        %             mShuffledCatchProb = squeeze(nanmean( ppantShuffled,1)); %or
        % take mean over trials.
        %will plot 99%CI of the shuffled mean across trials.
        acrossall_mShuffledCatchOffset_MD(ippant,:,:) = (ppantShuffled);
    end
    
    cd(basefol)
    cd('newplots-MD')
    save('ShuffledCatch', 'acrossall_mShuffledCatchOffset_MD', '-append')
end


if job.plotCatchtracesbynumTargets_eachppant==1;
 cd(basefol)
    cd('newplots-MD')
    load('Catchperformance.mat')
    load('ShuffledCatch.mat')
    %collect onset trace across ppants for each type
    %%
    clf
    %collect RTs for later.
    onsetRTs=[];offsetRTs=[];
    for onsetoroffset=1%:2
        %reset outgoing variable

        figure(11);    
%         clf
        
        if onsetoroffset==1
            
            dtAll= squeeze(nanmean(mean_catchBP_Probtraceperppant_bynum(:,1:3,:),2));% just use 1-3 targets
            x=xtONSET;
       eventis='onset';
            
            col = ['r'];
        else
            dtAll= squeeze(nanmean(mean_catchBP_ProbtraceperppantOFFSET_bynum(:,1:3,:),2));
            x=xtOFFSET;
            col = ['b'];
            eventis='offset';
        end
        %%
%         ippant=2; is a clear example



        for ippant=2%:29
%             subplot(5,6,ippant)
    
    RT=[];
    %% plot results of shuffled
    
    if onsetoroffset==1
    acrossall_mShuffledCatch= squeeze(acrossall_mShuffledCatch_MD(:,1,:,:)); %has an extra dimension for onsets/offsets.
    else
        acrossall_mShuffledCatch= acrossall_mShuffledCatchOffset_MD; %new versio of offset shuffled likelihood (separate job - see above).
    end
    
        
        shuffdata= squeeze(acrossall_mShuffledCatch(ippant,:,:));

        %remove the nan trials, they mess up averaging.
        rmv = find(isnan(shuffdata(:,1)));
        shuffdata(rmv,:)=[];
        %%
%         mShuff = squeeze(nanmean(shuffdata,1)); %m across trials (200 shuffs per timepoint).
        
        
        Shuffupperbound=[];
        Shufflowerbound=[];
        % for each time point, calculate width of distribution, based on CI of mean:
        for ip = 1:length(shuffdata)
          %%
            shd= shuffdata(:,ip)';
%with logit transform
%             shd(shd==0)=.0001;
%             %note not normally distributed!
%             %try logit transform
%             lp = log(shd)-log(1-shd);
%                 
%             %now rearrange so that we can place correct scores.
%               [lpn,id]=sort(lp);%convert to zscore
%               %sort first.
%             tz= zscore((lpn));
%                         
%             %find closest to zscores to 99% CI (2.57), 95% CI= 1.96;
% %             AB=dsearchn(tz, [-2.57 2.57]');
%             AB=dsearchn(tz', [-1.96 1.96]');              
%             %this is upper and lower bound.
%             CI= [lpn(AB(1)), lpn(AB(2))];
%             
%             %now convert back from logit.
% %CI upper bound (95%)
% Shuffupperbound(1,ip) = exp(CI(2)) ./ (exp(CI(2))+1);
% Shufflowerbound(1,ip) = exp(CI(1)) ./ (exp(CI(1))+1);
%           
%without
shd=sort(shd);
tz=zscore(shd);
AB=dsearchn(tz', [-1.96, 1.96]');
  CI= [shd(AB(1)), shd(AB(2))];
Shuffupperbound(1,ip) = CI(2);
Shufflowerbound(1,ip) = CI(1);

        end
        
     
        %% plot shuffled
  shf=plot(x, squeeze(nanmean(shuffdata,1)), 'color', [.1 .1 .1], 'linew', 3);
            %add patch
            %Calculate the error bars
            uE=Shuffupperbound;
            lE=Shufflowerbound;
            
            %Make the patch
            yP=[lE,fliplr(uE)];
            xP=[x,fliplr(x)];
            
            %remove nans otherwise patch won't work
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            hold on
            colp=[.1 .1 .1];
            H.patch=patch(xP,yP,1,'facecolor',colp,...
                'edgecolor','none',...
                'facealpha',.15);
            
            
            %Make pretty edges around the patch.
            H.edge(1)=plot(x,lE,'-','color',colp);
            H.edge(2)=plot(x,uE,'-','color',colp);
             xlim([-1 4])
        
        %%
        % now repeat the process for observed data.
        hold on    

            dt = squeeze(dtAll(ippant,:));
            if onsetoroffset==1
            pl=plot(xtONSET, dt, 'color', [col], 'linew', 3);

            else
                pl=plot(xtOFFSET, dt, 'color', col, 'linew', 3);
            end

            hold on
            
            
            %store real catch trace across avail trials for comparison of
            %significant departure from shuffled.
%             anovadata(2,1:size(stack,1),:)= stack;
            p=[];
            
% %             %% perform anovas per timepoint.
% %             subs = repmat(1:22 ,1, 4)';
% %             trace= [ones(1,22)'; 2*ones(1,22)'; 3*ones(1,22)'; 4*ones(1,22)'];
%%           
%%
% %         %plot first sig point after onset on graff
% %         %fdr correct

% %         
            if onsetoroffset ==1
                %find first point data lower is greater than shuff upper
%             sigpoint=find((dt(120:end))>Shuffupper(120:end), 1 ); %finds first.
            sigpoints=find(dt>Shuffupperbound); %finds first.
            sigpoint = min(sigpoints(sigpoints>window(2)));
            placement = (sigpoint)/60 -abs(window(1)/60);
           
            try onsetRTs(ippant)=placement;
            catch
                onsetRTs(ippant)=nan;
            end
            placeY=0.95;
            colis= [0 .5 0];
            
            else %for offsets, we want the first ns point.
%             sigpoint=find(dt(480:end)<Shuffupper(480:end), 1 ); %finds first.
% Shufflower = dt - stErrShuff/2;            
sigpoints=find(dt<Shufflowerbound);
            
            %first after offset
            sigpoint = min(sigpoints(sigpoints>window(2))); 
            
            placement = (sigpoint)/60 -abs(window(1)/60) ;
            try offsetRTs(ippant)=placement;
                
                if placement<.02
                    
                    offsetRTs(ippant) =nan;
                end
            catch
                offsetRTs(ippant)=nan;
            end
            placeY=0.95;
            colis = 'r';
            end
            
            % now plot
            
%             try sigt=plot(placement, placeY, '*','color', colis, 'markersize', 5, 'linewidth', 3);
                try plot([placement placement], [0 1], '--','color', 'b', 'linewidth', 1);
                placementT=1;
            catch
                placementT=0;
                end
        
            
                

% title(['Participant ' num2str(ippant)]);
            ylim([ 0 1])
%         xlabel('onset')
%%
hold on
plot([ 0 0], ylim, ['k-'], 'linew', 1)
        ylabel('Button press likelihood')
        set(gca, 'fontsize', 15)
        shg
        end
        %%
    legend([pl shf], 'Observed', 'Baseline')
    %%
%     st=suptitle('Reaction Time to Catch events');
%         st.FontSize=30;
%         st.FontWeight='bold';
set(gcf, 'color', 'w')    
cd(basefol)
    cd('newplots-MD')
    cd('Figures')
    cd('Catch analysis')
    %%
    
    end
    
    
    % if plotting single ppant:
      xlabel(['Time from PMD ' eventis ' [sec]'], 'fontsize', 15)
        
        ylabel('Button press likelihood', 'fontsize', 15)

        set(gca, 'fontsize', 1.5*fontsize)
        
    %%
    print('-dpng', 'Timecourse of each CatchDisap BP, by ppant (catch combined')
    
    
    %% also plot bar
    figure(1);
    clf
    goodppants=allppants;
% goodppants=1:29
    plotbar=[nanmean(onsetRTs(goodppants)); nanmean(offsetRTs(goodppants))];
    plotst = [nanstd(onsetRTs(goodppants)); nanstd(offsetRTs(goodppants))];
    bh=bar(plotbar);
    bh.FaceColor= ['r'];
%     bh.FaceColor= [0 .5 0];
    %
    hold on
    bh = bar(2, plotbar(2), 'r');
    shg
    %%
    hold on
    
    % adjust within subject error bars.
    tmpX = [onsetRTs(goodppants); offsetRTs(goodppants)];
    
    %within subj mean
    mX = nanmean(tmpX,1);
    %overallm
    mmX=nanmean(mX); 
    
    %adjust. 
    newX= tmpX - repmat(mX,[size(tmpX,1),1]) + repmat(mmX, [size(tmpX,1),size(tmpX,2)]);
    
    plotst = nanstd(newX,0,2)/sqrt(size(newX,2));
    %
    
    e=errorbar(plotbar, plotst, 'color', 'k'); 
    e.LineWidth=2;
    e.LineStyle='none';
    %%
    set(gca, 'fontsize', 15*1.5')%, 'xtickmark', {'Onset' 'End'})
    set(gca, 'XTickLabel', {'Onset' 'Offset'})
    ylabel('reaction time [sec]')
    set(gcf, 'color', 'w')
    ylim([0 1.2])
%     axis square
    shg
    %%
    print('-dpng', 'Catch RT across participants');
end


if job.plotCatchtracesbynumTargets==1
    cd(basefol)
    cd('newplots-MD')
    load('Catchperformance.mat')
    load('ShuffledCatch.mat')
    
    colis='r';
    %%% %%% for these plots, we can use CI or SEM for error bars
    useSEMorCI=2;
    
    
    
    %collect onset trace across ppants for each type
    %%
    
    
    
   goodppants=allppants;

   %    goodppants=badppants(1:3);%

    RT=[];
    % plot results of shuffled
    figure(1);
    clf
    for onsetoroffset = 1%:2
        
        anovadata = zeros(length(goodppants),4, length(xtOFFSET)); %for comparing each time point.
        
        if onsetoroffset==1
            stack=mean_catchBP_Probtraceperppant_bynum;
             acrossall_mShuffledCatch= squeeze(acrossall_mShuffledCatch_MD(:,1,:,:)); %has an extra dimension for onsets/offsets.
                  
            x = xtONSET;
                 eventis= 'onset' ;
        else
            stack=mean_catchBP_ProbtraceperppantOFFSET_bynum;
            acrossall_mShuffledCatch= acrossall_mShuffledCatchOffset_MD; %new versio of offset shuffled likelihood (separate job - see above).
                 
            x = xtOFFSET;
                 eventis= 'offset' ;
        end
        
        for idatatype=1:2
            switch idatatype
                case 1 %plot shuffled first
        %mean across trials, within participants.
        tmpA= squeeze(acrossall_mShuffledCatch(goodppants,:,:));
        
        %mean across shuffs( 26 ppants).
        tmp= squeeze(nanmean(tmpA,2));
        shuffdata=tmp;
        
        
        mShuff = squeeze(nanmedian(tmp,1)); %m across ppants.
      
   colsh=[.1 .1 .1];
   
                case 2
                       
            dt= squeeze(nanmean(stack(goodppants,1:3,:),2));
            shuffdata = dt;
            
            mShuff = squeeze(nanmedian(dt,1));
            colsh='r';
            end
                    
           % calculate CI         
        Shuffupperbound=[];
        Shufflowerbound=[];
        stErrShuff=[];
        % for each time point, calculate width of distribution:
        for ip = 1:length(shuffdata)
                     
            %%
            shd= shuffdata(:,ip)';
            shd(shd==0)=.0001;
            %note not normally distributed!
            %try logit transform
            lp = log(shd)-log(1-shd);
                
            %now rearrange so that we can place correct scores.
              [lpn,id]=sort(lp);%convert to zscore
              %sort first.
            tz= zscore((lpn));
                        
            %find closest to zscores to 99% CI (2.57), 95% CI= 1.96;
%             AB=dsearchn(tz, [-2.57 2.57]');
            AB=dsearchn(tz', [-1.96 1.96]');              
            %this is upper and lower bound.
            CI= [lpn(AB(1)), lpn(AB(2))];
            
            %now convert back from logit.
%CI upper bound (95%)
Shuffupperbound(1,ip) = exp(CI(2)) ./ (exp(CI(2))+1);
Shufflowerbound(1,ip) = exp(CI(1)) ./ (exp(CI(1))+1);
          

        end
        
     
        %%

        plot(x,mShuff, 'color', colsh, 'linew', 3);
       % if plotting bad ppants.
%         plot(x,dt', 'color', 'b', 'linew', 3);
        hold on
            %add patch
            %Calculate the error bars
            uE=Shuffupperbound;
            lE=Shufflowerbound;
            
            %Make the patch
            yP=[lE,fliplr(uE)];
            xP=[x,fliplr(x)];
            
            %remove nans otherwise patch won't work
            xP(isnan(yP))=[];
            yP(isnan(yP))=[];
            
            
            H.patch=patch(xP,yP,1,'facecolor',colsh,...
                'edgecolor','none',...
                'facealpha',.15);
            
            
            %Make pretty edges around the patch.
            H.edge(1)=plot(x,lE,'-','color',colsh);
            H.edge(2)=plot(x,uE,'-','color',colsh);
             xlim([-1 4])
%         
        %%
        % now repeat the process for observed data.
        hold on 
     
        end
        
        ylim([0 1])
        cr=plot([0 0 ], ylim, ['k' '-'], 'linewidth', 1);
%         plot(xlim, [0 0 ], ['k', '-'], 'linewidth', 2);
        xlabel(['Time from PMD ' eventis ' [sec]'], 'fontsize', 15)
        
        ylabel('Button press likelihood', 'fontsize', 15)

        set(gca, 'fontsize', 1.5*fontsize)
        %     ylim([0 1])
    end
    set(gcf, 'color', 'w')
    %%
    cd(basefol)
    cd('newplots-MD')
    cd('Figures')
    cd('Catch analysis')
    %%
    set(gcf,'color', 'w')
    
    print('-dpng', ['Timecourse of each CatchDisap BP, across all with CI'])
end
%%
if job.calcMissedcatchesbyTargperppant==1
   cd(basefol)
    cd('newplots-MD')
    %%
    load('MD_AllBP_perppant.mat')
    pD=ppantTrialDatawithDetails;
    catchStruct=[{'ppant'}, {'trial'}, {'catchtype'}, {'failed'}, {'Locations'}, {'RejectTrial'}];
%     req_frames=60;
    catchData_numbergone=[];
    catchData_location=[];
    
    %manually identified 'false' switches to ignore
    %ppant trial % These are when no fourth button could be pressed, on
    %account of three other buttons already being pressed.
    
    falseswitches = [12 17; ...
        16 19;...
        16 22; ...
        17 15; ...
        19 17; ...
        21 11; ...
        22 17; ...
        26 5];
        
    
    for ippant=1:29
        counter1=1;
        counter2=1;
        counter3=1;
        counter4=1;
        loc1=1;
        loc2=1;loc3=1;loc4=1;
        
        for itrial=1:24
            catchtype=pD(ippant).TrialDetails(itrial).TotalCatchTargetsRemoved;
            chData=sum(pD(ippant).TrialDetails(itrial).CatchBPs_exact,2);
            
            totalabsence = pD(ippant).TrialDetails(itrial).totalCatchdur; 
            
            % for both, adjust for RTs of 1 second.
            totalabsence = totalabsence-60;

            
            
            missedcatches=length(find(chData< totalabsence*.50));
            missedindx = find(chData<totalabsence*.5);
            %catch needed, find locations we missed.
            
            catchneeded = [pD(ippant).TrialDetails(itrial).TL_Catch, pD(ippant).TrialDetails(itrial).TR_Catch,pD(ippant).TrialDetails(itrial).BL_Catch,pD(ippant).TrialDetails(itrial).BR_Catch];
            %index of which buttonp presses were catch related
            bpindex=find(catchneeded>0);
            if find(chData<totalabsence*.5)>0 % if we had a miss.
                
                locationsmissed=bpindex(missedindx);
            else
                locationsmissed=0;
            end
            
            
            
            %record for later analysis in matrix.
            switch catchtype
                case 1
                    catchData_numbergone(ippant,counter1, catchtype) = missedcatches;
                    counter1=counter1+1;
                    
                    %Also include an index for whether to reject the trial or not.
                    
                    if missedcatches>0
                        rejectTrial=1; %failed single target catch.
                    else
                        rejectTrial=0;
                    end
                    
                case 2
                    catchData_numbergone(ippant,counter2, catchtype) = missedcatches;
                    counter2=counter2+1;
                    
                    
                    if missedcatches>1
                        rejectTrial=1; %failed dual target catch, but if they reported one, good enough.
                    else
                        rejectTrial=0;
                    end
                case 3
                    catchData_numbergone(ippant,counter3, catchtype) = missedcatches;
                    counter3=counter3+1;
                    
                    if missedcatches>2
                        rejectTrial=1;
                    else
                        rejectTrial=0;
                    end
                    
                case 4
                    catchData_numbergone(ippant,counter4, catchtype) = missedcatches;
                    counter4=counter4+1;
                    
                    if missedcatches>3
                        rejectTrial=1;
                    else
                        rejectTrial=0;
                    end
                    
                case 0
                    missedcatches=0;
                    locationsmissed=0;
                    rejectTrial=0;
            end
            
            
            % input data to structure.
            %override if this trial is a false switch (detected manually);
            if any(falseswitches(:,1)==ippant)
                
                tid=find(falseswitches(:,1)==ippant);
                %check if this trial needs rejecting.
                for ich=1:length(tid)
                    
                    if falseswitches(tid(ich),2)==itrial;
                        rejectTrial=nan;
                    end
                end
            end
            catchStruct=[catchStruct; {num2str(ippant)}, {num2str(itrial)}, {num2str(catchtype)}, {num2str(missedcatches)}, {num2str(locationsmissed)}, {num2str(rejectTrial)}];
            
            
            
        end
        %%
    end
    %%
    req_frames= ['.5*totalabsence'];
    
    save('Catchperformance', 'catchStruct', 'catchData_numbergone', 'req_frames', '-append')
end
if job.plotmissedcatchesbynumremoved==1
    %% unpack
    %mean per ppant by targ type
    for ipl=2%1:3
        switch ipl
            case 1 %use 1:4 targs
                catchDatapl = squeeze(mean(catchData_numbergone(goodppants,:,:),2));  %average over trials
                bartoplot= squeeze(mean(catchDatapl,1)); %w/in ppants.
                errplot=std(catchDatapl)/sqrt(size(catchDatapl,1));
                printname= 'CatchTargets missed 1 to 4';
            case 2 %use 1:3 no 4.
                catchDatapl = squeeze(mean(catchData_numbergone(goodppants,:,1:3),2));
                bartoplot= squeeze(mean(catchDatapl,1)); %w/in ppants.
                errplot=std(catchDatapl)/sqrt(size(catchDatapl,1));
                printname='CatchTargets missed 1 to 3';
            case 3 %use combined 3 and 4.
                catchDatapl = squeeze(mean(catchData_numbergone(goodppants,:,:),2));
                bartoplot1= squeeze(mean(catchDatapl,1)); %w/in ppants.
                errplot1=std(catchDatapl)/sqrt(size(catchDatapl,1));
                
                bartoplot=bartoplot1(1:2); %need to calculate separately for third column.
                errplot=errplot1(1:2);
                %vector of all tri
                bartoplot(3)=mean(bartoplot1(3:4));
                errplot(3)=mean(errplot1(3:4));
                printname= 'Catch Targets missed 1 to comb';
        end
        
        clf
        bh= bar(bartoplot, 'k');
        hold on
        errorbar(bartoplot, errplot, 'linestyle', 'none', 'color', 'r', 'linewidth', 3)
        %         ylim([0 3])
        xlabel({['Numer of targets removed'];['(Catch Events)']}, 'fontsize', 15);
        if ipl==3
            set(gca,'Xticklabel', [{'1'}, {'2'}, {'>2'}])
        end
        
        title({['Catch targets missed by participants across all trials, by type of event'];['with standard error of the mean across participants']}, 'fontsize', 15)
        ylabel('Single targets missed', 'fontsize', 15)
        ylim([0 .75])
        set(gca, 'fontsize', 15)
        %%
        cd(basefol)
        cd('newplots-MD')
        cd('figures')
        cd('Catch analysis')
        %%
        print('-dpng', printname)
        %
    end
    
    %Total number of removed targets = 60 per ppant (excluding 4 catch trials).
    %total missed per ppant:
    allmissed1= squeeze(sum(catchData_numbergone(goodppants,:,:),2));
    %from 1:3 types
    allmissed2=sum(allmissed1(:,1:3),2);
    totalacrossppants= sum(allmissed2);
    percentmissed = totalacrossppants/((6+12+18)*26);
    %%
    
end
%%
%createstructure per ppant, to see if location effect
if job.calcmissedcatches_bylocationperppant==1
    s_loc=[];
    for ippant = 1:29
        indx = [(2:25) + (24*(ippant-1))];
        
        
        %% all missed catches by location
        allmissed=[];
        for r=1:24
            if (str2num(catchStruct{indx(r),3})) ~= 4 %if not a 4 target removal.
                allmissed= [allmissed, str2num(catchStruct{indx(r),5})]; %if a legit trial, store the location of missed catches.
            end
        end
        
        s_loc(ippant).allmissed= allmissed; %allmissed should now be 18 trials in length( removing the 6, 4 catch trials)
        s_loc(ippant).countmissed_TL = length(find(allmissed==1));
        s_loc(ippant).countmissed_TR = length(find(allmissed==2));
        s_loc(ippant).countmissed_BL = length(find(allmissed==3));
        s_loc(ippant).countmissed_BR = length(find(allmissed==4));
    end
    
    %%
    %prepare for ANOVA.
    for ippant=1:29
        datab(ippant,1)=s_loc(ippant).countmissed_TL;
        datab(ippant,2)=s_loc(ippant).countmissed_TR;
        datab(ippant,3)=s_loc(ippant).countmissed_BL;
        datab(ippant,4)=s_loc(ippant).countmissed_BR;
    end
    catchData_bylocation = datab;
    catchStruct_bylocation=s_loc;
%     if strcmp(pwd, '/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy/newplots-MD')~=1
%         cd('newplots-MD')
%     end
    %%
    save('CatchPerformance', 'catchData_bylocation', 'catchStruct_bylocation', '-append');
end

if job.Plotmissedcatches_bylocationperppant==1
    %%%%%%%
    %%%%%%%%%%%%%%%%%% Unfinished
    
    
    %%
    bartoplot= squeeze(mean(catchData_bylocation,1)); %w/in ppants.
    errplot=std(catchData_bylocation)/sqrt(size(catchData_bylocation,1));
    printname= 'All CatchTargets missed by location';
    %%
    
    clf
    bh= bar(bartoplot, 'k');
    hold on
    errorbar(bartoplot, errplot, 'linestyle', 'none', 'color', 'r', 'linewidth', 3)
    %         ylim([0 3])
    xlabel({['Location of targets removed'];['(Catch Events)']}, 'fontsize', 15);
    
    set(gca,'Xticklabel', [{'TL'}, {'TR'}, {'BL'}, {'BR'}])
    
    
    title({['Catch targets missed by participants across all trials, by location'];['with standard error of the mean across participants (n=' num2str(length(goodppants)) ')']}, 'fontsize', 15)
    ylabel('Single targets missed', 'fontsize', 15)
    %         ylim([0 2])
    set(gca, 'fontsize', 15)
    %%
    cd(basefol)
    cd('newplots-MD')
    cd('Figures')
    cd('catch analysis')
    %%
    print('-dpng', printname)
    %
    
    %% plot freq of disappearances per location.
    
end
