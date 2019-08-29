% s3_Cc_takeRESSchunked

%take two for RESS, focusing on smaller windows to avoid transients.
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/MatthewDavidson/Desktop/FromIrenes Mac/Data Experiment copy')
end
%%
basefol=pwd;

clearvars -except basefol allppants
dbstop if error
%based on/replaces Mkat_MT_taper from Irene.



%%

job.sortTrialindicesbyFreqandLocation=0; %store relevant catch/target hz data.

job.EpochperppantCATCH=0; %Epoch each participants Catch windows
%also stores the BP trace for next job.

job.ppantCatch_topotime=1; %used to append multichannel data for spatial correlation.
job.concatTOPOTIMEacrossppants=1; %concat the above
job.plotTOPOtimeacrossppants=1;

%also plot topotime ver of catch.


job.erpimageCATCHusingSNR=0; % saves at ppant level.


job.concaterpdataacrossppants=0;

job.erpimageacrossppants=0;
job.BPandSSVEPtimecourseacrossppants=0;


% job.calcRESSforcatchtrials=1; % also performs a baseline removal.
% job.concatacrossppants=1;
% job.plotacrossppants=1;




%based on/replaces Mkat_MT_taper from Irene.
cd(basefol)
cd('newplots-MD')
%% load data to determine physical catch timing.
load('MD_AllBP_perppant')
load('Catchperformance') ; %to identify failed catches. (no BP).
%%

%remaining participants after behavioral data analysis/exclusion
allppants=1:29;
badppants = [8, 15,28, 4 , 7,5 6];
allppants(badppants)=[];

% window=[-2 4];
window=[-3 3];

srate=250;

rmvbaseEPOCH=1; %0 for no removal,  1 for trial baseline, 2 for pooled freq baseline (compare across conds).
baserem1 = [-2.99 -2]; % baseline to remove in secs for ONSET of catch
baserem2 = baserem1;%[2 2.99]; % baseline to remove in secs for offset of catch
epochdur = sum(abs(window))*srate;
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;


if job.sortTrialindicesbyFreqandLocation==1
    
    if job.sortTrialindicesbyFreqandLocation==1
        %% This stores the 'like' target location x Hz combinations, for RESS analysis across like trials,
        % which needs spatial stability.
        %%
        cd(basefol)
        cd ./
        load('Raw_Behavioral_AllPP.mat')
        
        %%
        for ifol=allppants
            cd(basefol)
            %change into ppant folders
            cd(num2str(ifol))
            
            for Loc = 1:4 %for each location, (TL,TR,BL,BR) , store trial indices
                %by Hz, (8, 13,15,18)
                % 24 rows per participant.
                corr_rows = [1:24] + (24*(ifol-1));
                
                targindx=zeros(4,6);
                catchindx=zeros(4,6);
                
                for ihz=1:4
                    switch ihz
                        case 1
                            searchhz = 8;
                        case 2
                            searchhz = 13;
                        case 3
                            searchhz = 15;
                        case 4
                            searchhz = 18;
                    end
                    %the columns for target location start at 2. Find when Hz
                    %is in each location
                    ftrials=find(AllTrialsAllPPBottomLeft(corr_rows,(Loc+1))== searchhz);
                    targindx(ihz,:)= ftrials;
                    
                    
                    
                    % also whenn this
                    % hzxloc was removed.
                    for icatch=1:length(ftrials)
                        trialid = ftrials(icatch);
                        if AllTrialsAllPPBottomLeft(corr_rows(trialid),(Loc+5))==1
                            
                            catchindx(ihz,icatch) = trialid;
                            
                        end
                    end
                    
                    
                    
                    
                end
                
                switch Loc
                    case 1
                        TopLeftTrialindexbyHz =targindx;
                        TopLeftCatchindexbyHz=catchindx;
                    case 2
                        TopRightTrialindexbyHz =targindx;
                        TopRightCatchindexbyHz=catchindx;
                    case 3
                        BottomLeftTrialindexbyHz =targindx;
                        BottomLeftCatchindexbyHz=catchindx;
                    case 4
                        BottomRightTrialindexbyHz =targindx;
                        BottomRightCatchindexbyHz=catchindx;
                        
                end
                
            end
            save('TrialIndicesbyLocationandHz', 'TopLeftTrialindexbyHz', ...
                'TopRightTrialindexbyHz', 'BottomLeftTrialindexbyHz',...
                'BottomRightTrialindexbyHz','TopLeftCatchindexbyHz', ...
                'TopRightCatchindexbyHz', 'BottomLeftCatchindexbyHz',...
                'BottomRightCatchindexbyHz')
        end
    end
    
    
end
if job.EpochperppantCATCH==1
    %%
    cd(basefol)
    cd('newplots-MD')
    load('Catchperformance.mat', 'catchStruct') % note that for some trials no catch occured (and screen froze!).
    %%
    for ifol = allppants
        cd(basefol)
        cd(num2str(ifol))
        %load new rereferenced data.
        %         load('ppantRESSwholetrial')
        load(['P' num2str(ifol) '_autopreprocd.mat'])
        %there are separate types of catch we will store:
        %.1=all catches at on/offset (regardless of response),
        %.2=failed catches at on/offset (target absence no response),
        %.3= BP centred catches (changes timing).
        
        %.1
        ppant_SNREEG_catchramponset=nan(24,64,epochdur+1);
        ppant_SNREEG_catchrampoffset = nan(24,64,epochdur+1);
        
        %.3
        ppant_SNREEG_disapBPwithincatch = nan(24,64,epochdur+1);
        ppant_SNREEG_reapBPaftercatch = nan(24,64,epochdur+1);
        
        %.2
        ppant_SNREEG_Disapfailedcatches=nan(24,64,epochdur+1);
        ppant_SNREEG_Reapfailedcatches=nan(24,64,epochdur+1);
        
        
        %how many did this ppant pass/fail
        ppantrows = [1:24] + (24*(ifol-1)) + 1;
        %         tmp=str2num(cell2mat(catchStruct(ppantrows,6)));
        %         goodtrials = 24- sum(tmp);
        %         badtrials = sum(tmp);
        
        
        goodcounter=1;
        badcounter=1;
        
        catchonsetBPs=nan(24,361);
        catchoffsetBPs=nan(24,361);
        withincatchonsetBPs=nan(24,361);
        postcatchoffsetBPs=nan(24,361);
        
        firstBPwithincatch_RTs=nan(24,1);
        firstBPaftercatch_RTs=nan(24,1);
        
        for itrial=1:24
            
            
            ppant_SNREEG_catchramponset_tmp=[];
            ppant_SNREEG_catchrampoffset_tmp=[];
            ppant_SNREEG_disapBPwithincatch_tmp=[];
            ppant_SNREEG_reapBPaftercatch_tmp=[];
            
            
            
            
            
            
            
            trialdata= ppantTrialDatawithDetails(ifol).TrialDetails(itrial);
            
            catchtruestart = trialdata.Catchstart_frames/60;
            catchtrueend =  trialdata.Catchend_frames/60; %perhaps remove the ramp?
            
            
            %which target frequencies were the 'catch' in this trial?
            catchlocs = [trialdata.TL_Catch, trialdata.TR_Catch, trialdata.BL_Catch, trialdata.BR_Catch];
            
            
            
            
            
            % take 64 channel data for this EPOCH.
            tmpEpoch = squeeze(EEG_preprocd(:,:,itrial)); %both 20 and 40 hz.
            %%
            
            %initialize
            tmpCatchONsetOFFset=[];
            if sum(catchlocs)~=0 %skip trial.
            %onset and offset EEG epochs. (not accounting for BP)
            for iwin=1:2
                switch iwin
                    case 1
                        centrep = catchtruestart;
                    case 2
                        centrep = catchtrueend;
                end
                
                windowstart = centrep-abs(window(1));
                %                 windowend = centrep+(window(2));
                
                [~, rampst] = min(abs(tt-windowstart));
                
                rampend=rampst+epochdur; %add epoch length.
                
                
                %also collect BP data for comparative plots (next job).
                
                switch iwin
                    case 1
                        ppant_SNREEG_catchramponset_tmp = tmpEpoch(:,rampst:rampend);
                        
                    case 2
                        ppant_SNREEG_catchrampoffset_tmp = tmpEpoch(:,rampst:rampend);
                        
                end
            end
            
           
                catchlocs= find(catchlocs);
                
                % now epoch around BP , first within catch, then after.
                % Was there a BP for catch? ie catch (pass/fail)?
                
                rowis=itrial + (24*(ifol-1)) + 1;
                badcatch = str2num(catchStruct{rowis, 6}); %6th column for reject or no.
                
                if badcatch~=1 %then collect BP centred data also
                    
                    CatchBPs = nansum(trialdata.CatchBPs_exact,1);
                    
                    %if BPdata during catch removal
                    firstBPwithincatch= min(find(diff(CatchBPs)>0));
                    
                    if firstBPwithincatch>0 %change in BP state (not already pressed)
                        
                        
                        
                        
                        %Epoch around this point (notice catch gone)
                        
                        %                     correct for difference in timing units.
                         % note that the ramp was also present, and 'catchrampstart' which delays total removal of flickering target.
                        BP_disapwithincatch_secs = firstBPwithincatch/60 + catchtruestart;
                        
                        [~, BPdisap1] = min(abs(tt-(BP_disapwithincatch_secs-abs(window(1)))));
                        BPdisap2 = BPdisap1+epochdur;
                        
                        ppant_SNREEG_disapBPwithincatch_tmp= tmpEpoch(:,BPdisap1:BPdisap2);
                        
                        firstBPwithincatch_RTs(itrial)=firstBPwithincatch;
                        
                        
                        %store associated button press Epoch.
                        try withincatchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(BP_disapwithincatch_secs*60-window(2)*60):ceil(BP_disapwithincatch_secs*60+window(2)*60)),1));
                        catch
                            withincatchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(BP_disapwithincatch_secs*60-window(2)*60):ceil(BP_disapwithincatch_secs*60+window(2)*60+1)),1));
                        end
                        
                        
                        
                        
                    else %no response to catch?
                        ppant_SNREEG_disapBPwithincatch_tmp= nan(64,epochdur+1);
                        firstBPwithincatch_RTs(itrial)=nan;
                    end
                    
                    %%%%% also collect and epoch around first BP after catch
                    %%%%% return.
                    
                    %note that the 'long epoch' is 4s before catch start and
                    %4seconds after catch START
                    
                    catchdur=trialdata.totalCatchdur;
                    
                    longcatchEpoch=nansum(trialdata.CatchBPsOFFSET_longepoch,1);
                    
                    %midpoint, true catch offset
                    mp=ceil(length(longcatchEpoch)/2);
                    postOffset=longcatchEpoch(1,(mp):end);
                    
                    %first change in BP, (could also use return to zero?)
                    
                    firstBPaftercatch= min(find(diff(postOffset)<0));
                    
                    %                     firstBPaftercatch= min(find(postOffset~=0));
                    
                    %ajust this time, to account for seconds in whole trial EEG
                    
                    catchreturn_secs = trialdata.Catchend_frames/60 + firstBPaftercatch/60;
                    
                    try %will fail if BPdisap2 is beyond the end of the trial.
                        [~, BPdisap1] = min(abs(tt-(catchreturn_secs-abs(window(1)))));
                        BPdisap2 = BPdisap1+epochdur;
                        
                        ppant_SNREEG_reapBPaftercatch_tmp= tmpEpoch(:,BPdisap1:BPdisap2);
                        
                        firstBPaftercatch_RTs(itrial)=firstBPaftercatch;
                        
                        
                        
                        %store BP aligned EPOCH
                        
                        
                        try postcatchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchreturn_secs*60-window(2)*60):ceil(catchreturn_secs*60+window(2)*60)),1));
                        catch
                            postcatchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchreturn_secs*60-window(2)*60):ceil(catchreturn_secs*60+window(2)*60+1)),1));
                        end
                        
                        
                    catch
                        ppant_SNREEG_reapBPaftercatch_tmp= nan(64,epochdur+1);
                        firstBPaftercatch_RTs(itrial)=nan;
                    end
                    
                    
                    
                    
                    
                    
                    goodcounter=goodcounter+1;
                else %just record the failed catch EEG, no BP window
                    
                    ppant_SNREEG_reapBPaftercatch_tmp= nan(64,epochdur+1);
                    ppant_SNREEG_disapBPwithincatch_tmp= nan(64,epochdur+1);
                    %                 badcounter=badcounter+1;
                end
                
                
                % STORE and save EEG per ppant.
                
                
                ppant_SNREEG_catchramponset(itrial,:,:)=ppant_SNREEG_catchramponset_tmp;
                
                
                ppant_SNREEG_catchrampoffset(itrial,:,:)=ppant_SNREEG_catchrampoffset_tmp;
                
                ppant_SNREEG_disapBPwithincatch(itrial,:,:)= ppant_SNREEG_disapBPwithincatch_tmp;
                
                ppant_SNREEG_reapBPaftercatch(itrial,:,:)= ppant_SNREEG_reapBPaftercatch_tmp;
                
                
                
                % store the BPs associated
                %catch aligned
                try catchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchtruestart*60-window(2)*60):ceil(catchtruestart*60+window(2)*60)),1));
                catch
                    catchonsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchtruestart*60-window(2)*60):ceil(catchtruestart*60+window(2)*60+1)),1));
                end
                
                try catchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchtrueend*60-window(2)*60):ceil(catchtrueend*60+window(2)*60)),1));
                catch
                    catchoffsetBPs(itrial,:)= squeeze(nansum(trialdata.BPData(catchlocs,ceil(catchtrueend*60-window(2)*60):ceil(catchtrueend*60+window(2)*60+1)),1));
                end
                
                
                
                
            end
            
        end
        
        savename = 'ppant_Catch_Epoched';
        
        save(savename, ...
            'ppant_SNREEG_catchrampoffset',...
            'ppant_SNREEG_catchramponset',...
            'ppant_SNREEG_disapBPwithincatch',...
            'ppant_SNREEG_reapBPaftercatch',...
            'window','tt', 'firstBPaftercatch_RTs', 'firstBPwithincatch_RTs', ...
            'catchonsetBPs', 'catchoffsetBPs',...
            'withincatchonsetBPs', 'postcatchoffsetBPs')
        
                disp(['Finished ppant ' num2str(ifol)])
    end
end

if job.ppantCatch_topotime==1
    
   getelocs
    window=[-3 3];
    
    srate=250;
    rmvbase=0;
    epochdur = sum(abs(window))*srate;
    
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
    
    onsetc = ceil(epochdur)/2;
    % peakfreqsare=[20,40]; %hz
    %timing
    tt = 0:1/srate:60;
    
    %%
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.Fpass= [0 50];
    param_spctrm.trialave=0;
    
    param_spcgrm.tapers = [1 1];
    param_spcgrm.Fs= [250];
    param_spcgrm.Fpass= [0 50];
    param_spcgrm.trialave=0;
    movingwin=[1,.15];
    
    kernelw = [-.25 -.25 0 0 1 0  0 -.25 -.25];
    
    %%
    for ifol = allppants
        for hzis=1:2
            switch hzis
                case 1
                    usehz=20;
                case 2
                    usehz=40;
            end
            
            
            
            icounter=1;
            
            cd(basefol)
            cd(num2str(ifol))
            
            load(['ppant_Catch_Epoched'])
            
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).
                    datatouse = ppant_SNREEG_disapBPwithincatch;
                   
                    ctype = 'report catch onset';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = ppant_SNREEG_reapBPaftercatch;
                    
                    
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                %             %rmvbaseline from EEG.
                datais=datatouse;
                bsdata = zeros(size(datais));
                
                for ichan=1:64
                for itrial = 1:size(datais,1)
                    td = detrend(squeeze(datais(itrial,ichan,:)), 'linear');
%                     tdrm= mean(td(1,1:250));
%                     rmb= repmat(tdrm, [1 length(td)]);
                    bsdata(itrial,ichan,:) = td;%-rmb;
                end
                end
                datais=bsdata;
                %%
                snr_grmout=zeros(64,size(datais,1),33);
                for ichan=1:64
                    datacr=squeeze(datais(:,ichan,:));
                [sgrm ,tgrm, fgrm] = mtspecgramc(datacr', movingwin, param_spcgrm);
                %%
                %conv SNR
                snr_sgrm =zeros(size(sgrm));
%                 kernelw=
                for itrial=1:size(sgrm,3)
                    %compute SNR
                    tmps= squeeze(sgrm(:,:,itrial));
                    
                    for itime= 1:size(tmps,1)
                        checkput = conv(log(tmps(itime,:)), kernelw,'same');
                        if ~isreal(checkput)
                            snr_sgrm(itime,:,itrial)= nan(1, size(tmps,2));
                        else
                            snr_sgrm(itime,:,itrial)= conv(log(tmps(itime,:)), kernelw,'same');
                        end
                    end
                end
                
                %store just stim freq.
                [~, hzid]= min(abs(fgrm-usehz));
                %%
                %reduce size.
                snrgrm20=squeeze(snr_sgrm(:,hzid,:));
                tbase= tgrm-3;
                %%
                if rmvbase==1
                    %rmv baseline ?
                    acrSNRsort_rmv=zeros(size(snrgrm20));
                    for itrial=1:size(snrgrm20,2);
                        
                        tmp=snrgrm20(:,itrial)';
                        %which baseline?
                        tidx= dsearchn(tbase', [bsrem]');
                        
                        %rmvbase
                        bs= mean(tmp(1,tidx(1):tidx(2)));
                        bs=repmat(bs, 1, length(tmp));
                        acrSNRsort_rmv(:,itrial)= tmp-bs;
                        %                         acrSNRsort_rmv(itrial,:)= tmp./bs;
                        
                    end
                    snrgrm20=acrSNRsort_rmv;
                end
                
%                 snr_grmout(ichan,:)=squeeze(nanmean(snrgrm20,2));
                snr_grmout(ichan,:,:)=snrgrm20';
                end
                %%
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        
                        Ppant_onsetSNR_allchan=snr_grmout; 
                        
                    case 2
                        
                        
                        Ppant_offsetSNR_allchan=snr_grmout; 
                        
                end
            end
            
            
            %%
            switch hzis
                case 1
                    savename='Catchperformance_withSNR_20,BPaligned';
                case 2
                    savename='Catchperformance_withSNR_40,BPaligned';
            end
            
            
            save(savename,...
                'Ppant_onsetSNR_allchan','Ppant_offsetSNR_allchan','-append');
                
            
        end
    end
end

if job.concatTOPOTIMEacrossppants==1
    for ihz=1:2
        
        if ihz==1
            loadname='Catchperformance_withSNR_20,BPaligned';
        else
            loadname='Catchperformance_withSNR_40,BPaligned';
        end
        storeacrossPpant_onsetSNR_chans=zeros(22,64,33);
        storeacrossPpant_offsetSNR_chans=zeros(22,64,33);
        
        icounter=1;
        
        for ippant = allppants
            cd(basefol)
            cd(num2str(ippant));
            %onset types
            
            load(loadname)

            storeacrossPpant_onsetSNR_chans(icounter,:,:)=squeeze(nanmean(Ppant_onsetSNR_allchan,2));
            storeacrossPpant_offsetSNR_chans(icounter,:,:)=squeeze(nanmean(Ppant_offsetSNR_allchan,2));
            
            icounter=icounter+1;
        end
        
        cd(basefol)
        cd('newplots-MD')
        
        
        switch ihz
            case 1
                savename=['GFX_Catchperformance_withSNR_20,BPaligned_allchan'];
            case 2
                savename=['GFX_Catchperformance_withSNR_40,BPaligned_allchan'];
        end
        
        save(savename, 'storeacrossPpant_onsetSNR_chans',...
            'storeacrossPpant_offsetSNR_chans', 'tgrm');
        
    end
    
end
   
    

%%

if job.plotTOPOtimeacrossppants==1
    clearvars -except job basefol
    getelocs;
    cd(basefol)
    
    cd('newplots-MD')
   
    %% %% %
    job2.plotSpacedTimetopo=0;
    job2.plotMeanTIMEtopo_andtvals=0;
    job2.plotSpatialCorrelation_overtime=1;
    %
%
%
%
%
%
%
%
%
    
    if job2.plotSpacedTimetopo==1
   cd(basefol)
   cd('newplots-MD')
        icount=1;
%     colormap('parula')
  for itimezero=1:2; %onset and offset
    for ihz=1:2

        switch ihz
            case 1
        load('GFX_Catchperformance_withSNR_20,BPaligned_allchan')
        titlep='1st harmonic';
            case 2
                load('GFX_Catchperformance_withSNR_40,BPaligned_allchan')
                titlep='2nd harmonic';
        end
      
            switch itimezero
                case 1
                    dPLOT = storeacrossPpant_onsetSNR_chans;
                    TIMING = [0];
%                     titlep=
                case 2
                    dPLOT = storeacrossPpant_offsetSNR_chans;
                    TIMING = [0];
            end
        
        %
%     TIMING = [  -.5 -.4 -.2 0 .2 .4 .5];
            figure(1)
    
    timeIND = dsearchn([tgrm-3]', TIMING');
    %%
    hold on
    pvals=nan(64,length(timeIND));
    tvals=nan(64,length(timeIND));
    for itime=1:length(timeIND)
        
        subplot(2,2,icount)
        
        tid= timeIND(itime);
%         if itimezero==1
            ipl=itime;
%         else
%                 ipl=itime+length(timeIND);
%         end
           ipl= ipl + (length(TIMING)*(ihz-1));
        % plot sig. 
        timeTOPO=squeeze(dPLOT(:,:,tid));
        
        for ichan=1:64
            [~,pvals(ichan,itime),~, stat]= ttest(timeTOPO(:,ichan), 0, 'tail', 'both');
            tvals(ichan,itime)=stat.tstat;
                        
        end
        q=fdr(pvals(:),.05);
        pmask = (pvals(:,itime)<q);
        
        
%         subplot(2,length(timeIND),ipl)
%         subplot(1,2,ipl)
       tp= topoplot(squeeze(mean(dPLOT(:,:,tid),1)), elocs(1:64), 'pmask', pmask, 'conv', 'on'); 
       %else plot tscores
%        tp= topoplot(tvals, elocs(1:64), 'pmask', pmask, 'conv', 'on'); 
       
       caxis([0 1])
       set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
       
%         topoplot(squeeze(mean(dPLOT(:,:,tid),1)), elocs(1:64), 'emarker2', {find(pmask), 'o', 'w', 2}); caxis([-1 1])
%         if itime==length(timeIND)
            c=colorbar;
            
            ylabel(c, 'log(SNR)');
%         end
%         title([ sprintf('%.2f',(tgrm(tid)-3)) ' s'])
        title(titlep)
        set(gca, 'fontsize', 25)
        icount=icount+1;
    end
        end
c=colormap('jet');        %
c(1,:)=[ 0 0 0];
colormap(c)
%         subplot(3,1,3)
%         colorbar; caxis([-1 1])
        set(gcf, 'color', 'w')
    end
    %%
    cd('figures')
    print('-dpng', 'Offset - spatial correlation')
    end
    if job2.plotMeanTIMEtopo_andtvals==1
        %%
        cd(basefol)
        cd('newplots-MD')
%looking at PFIINC>PFIdec
        %         timewinds = [-3 -1.1]; %'20Hz early'
%         timewinds = [-.75 .25];%'20Hz late'

%         timewinds = [-3 -1.5]; %'40Hz early'
%         timewinds = [-.5 0];%'40Hz late'

%looking at baseline:
% timewinds=[-1 0.24]; %PFI increase 20Hz from baseline
% timewinds=[-1.28 0.1]; %PFI increase 40Hz from baseline
% 
% timewinds=[-.52 1.8]; %PFI decrease 20Hz from baseline
% timewinds=[-.5 .5]; %PFI decrease 40Hz from baseline
        ipl=1;
        clf
        
        %only use for offsets.
        rmvbase=0;
        
        
        
     for ihz=1%1:2%:2;
         cd(basefol)
        cd('newplots-MD')
        switch ihz
             case 1
        load('GFX_Catchperformance_withSNR_20,BPaligned_allchan')
        titlep='1st harmonic';
            case 2
                load('GFX_Catchperformance_withSNR_40,BPaligned_allchan')
                titlep='2nd harmonic';
        end
        dPLOT=[];
                    dPLOT(1,:,:,:) = storeacrossPpant_onsetSNR_chans;
                 
                    
                    if rmvbase==1
                        bsrem=[-3 -1];
                        tbase=tgrm-3;
                        tmp = zeros(size(storeacrossPpant_offsetSNR_chans));
                        tIND = dsearchn(tbase', [bsrem]');
                        for ippant=1:size(storeacrossPpant_offsetSNR_chans,1)
                            for itrial = 1:size(storeacrossPpant_offsetSNR_chans,2)
                                tbase = squeeze(nanmean(storeacrossPpant_offsetSNR_chans(ippant,itrial,tIND(1):tIND(2)),3));
                                rbase = repmat(tbase, [1  size(storeacrossPpant_offsetSNR_chans,3)]);
                                
                                tmp(ippant,itrial,:) =squeeze(storeacrossPpant_offsetSNR_chans(ippant,itrial,:))' - rbase;
                            end
                        end
                        storeacrossPpant_offsetSNR_chans=tmp;
                    end
                    
                    dPLOT(2,:,:,:) = storeacrossPpant_offsetSNR_chans;
        
            
        
        %
%     TIMING = [  -.5 -.4 -.2 0 .2 .4 .5];
    
    timeIND = dsearchn([tgrm-3]', timewinds');
    %
    hold on
    %take mean over timewindow
    dTest= mean(dPLOT(:,:,:,timeIND(1):timeIND(2)),4);
    
    pvals=[];
    tvals=[];
    
               
        for ichan=1:64
            [~,pvals(ichan),~, stat]= ttest(squeeze(dTest(2,:,ichan)));
%             [~,pvals(ichan),~, stat]= ttest(squeeze(dTest(1,:,ichan)), squeeze(dTest(2,:,ichan)));
            tvals(ichan)=stat.tstat;                        
        end
        pmask = (pvals<.05);
        %
        topoplot(tvals, elocs(1:64), 'pmask', pmask)
        shg
        checkchans=find(pmask);
        % perform spatial cluster FDR:
            %% compare the increase in ASH to zero (since a diff already)
        % extract tvalues at spatially coincident p<05 electrodes,
        %(uncorrected)

%         checkchans= [25,62,30,61,29,63,31]; % PFI 20Hz, early.
%         checkchans= [29,30,61:63]; % PFI 20Hz, late.
% checkchans= [4 9 18 23 24 25 26 29 30 47 54 56 57 58 60 62 63];
% checkchans=[10,44,53,58,59,26,60:63,29:31];
        tvalsperchan = tvals(checkchans);
        
         % this the accrued test statistic (total), we have to 
% check against, after shuffling the data (Maris&Oostenveld)
         observedCV = sum(abs(tvalsperchan));
         %
        % now repeat the above process, but create random samples:
       sumTestStatsShuff = zeros(1,2000);
        for irand = 1:2000
            %testing the null that it isn't mismatched - matched at time 2 
            % which creates a diff. so select from either!
            shD= zeros(2,length(tvalsperchan),21);
            for ipartition = 1:2
                for ippant = 1:21
                    for chan=1:length(checkchans)
                    if mod(randi(100),2)==0 %if random even number
                    pdata = dTest(1,randi(21), checkchans(chan)); %select both chans

                    else %
                                            pdata = dTest(2,randi(21), checkchans(chan)); %select both chans

                    end
                
                shD(ipartition,chan,ippant) = pdata;
                    end
                end
            end
        
        %now compute difference between out hypothetical topoplots,
        % and test for sig, checking the accumulated test statistic at our
        % chans of interest
            tvalsperchan = zeros(1,length(checkchans));
            if length(checkchans>1)
            testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
            else
                testdata = (shD(1,:,:)) - (shD(2,:,:));
                testdata=testdata';
            end
            
            for itest = [1:length(checkchans)] %test each channel
               
                if length(checkchans)>1
                [~, p, ~,stat]= ttest(testdata(itest,:));
                else
                    [~, p, ~,stat]= ttest(testdata);
                end
                
                tvalsperchan(1,itest) = stat.tstat;
            end
            
            sumTestStatsShuff(1,irand) = sum(abs(tvalsperchan));
        end %repeat nshuff times
        %
        subplot(211)
        %plot histogram:
        H=histogram(abs(sort(sumTestStatsShuff)));
        % fit CDF
        cdf= cumsum(H.Data)/ sum(H.Data);
        %the X values (actual CV) corresponding to .01
        [~,cv05uncorr] = (min(abs(cdf-.95)));
        [~,cv01uncorr] = (min(abs(cdf-.99)));
        [~,cv001uncorr] = (min(abs(cdf-.999)));
        %THE Q VALUE FOR OUR OBSERVED DATA:
        
        
        hold on
         pCV=plot([observedCV observedCV], ylim, ['r-']);
         p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
         plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
         plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
         legend([pCV p05], {['observed'] ['p005'] })
         % what is the Pvalue? 
         
         %
         if observedCV>H.Data(cv05uncorr)
             %observed pvalue in distribution=
                    [~, c2] = min(abs(H.Data-observedCV)); %closest value.
                    pvalis= 1-cdf(c2);
             title({['Spatial Cluster  Significant!'];['Tvals = ' num2str(observedCV) ', p =' num2str(pvalis)]})
             %display Q
%              checkchans= [25,62,30,61,29,63,31];
             pmaskchecked = zeros(1,64);
             pmaskchecked(checkchans)=1;
         else
             title('Spatial Cluster  ns!')
         end
    
        
        
        
        
        
        
        
        %%
        subplot(2, 1, 2) ; 
% pmaskchecked(9)=[0]; PFI in 20 vs 0
% pmaskchecked(39)=[0]; %PIF in 40 vs 0

% pmaskchecked(45)=[0];pmaskchecked(19)=[0]; %PIF in 40 vs 0

        tp=topoplot(tvals, elocs(1:64), 'pmask', pmaskchecked, 'conv', 'on');
%         topoplot(tvals, elocs(1:64), 'emarker2', {find(pmask), 'o', 'w', 20}); 
        
        c= colorbar;
        ylabel(c, 't-value')
        caxis([-4 4])
        set(c, 'location', 'Eastoutside')
        set(gca, 'fontsize', 25)
        %%
        title({['\it p \rm\bf < .001 for'];[sprintf('%.2f',tgrm(timeIND(1))-3) ' to ' sprintf('%.2f',(tgrm(timeIND(2))-3)) 's']})
        %%
        set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('figures')
       
        shg
     end
     %%
    shg
    set(gcf, 'color', 'w')
     print('-dpng', 'Topo PreBP 20Hz late PFIin>PFIde.png')
    end
    
    
    
    
    
    
         if  job2.plotSpatialCorrelation_overtime==1;
             
             interpcatch=1;    
             
             cd(basefol)
        cd('newplots-MD')
             figure(1)
            clf
               
            selectchans=[1:64]; %use whole head
%             selectchans = [23:32,56:64]; % parieto-occipital only
        
        load('GFX_Catchperformance_withSNR_20,BPaligned_allchan')       
        onsetChans_20 = storeacrossPpant_onsetSNR_chans;
        offsetChans_20 = storeacrossPpant_offsetSNR_chans;
        
        load('GFX_Catchperformance_withSNR_40,BPaligned_allchan')
                 onsetChans_40 = storeacrossPpant_onsetSNR_chans;
               offsetChans_40 = storeacrossPpant_offsetSNR_chans;
                 
               %%
                clf
               leg=[];
               ttestdata= zeros(2, size(onsetChans_20,1),size(onsetChans_20,3));
               for iPFIdir=1:2%1:2 %onset and offset
                  hold on
                   switch iPFIdir
                       case 1
                           d1=onsetChans_20;
                           d2=onsetChans_40;
                           chis= 'catch target invisible';
                   %calculate spatial correlation over time:
                   colis='r';
                   linestyle='-';
                       case 2
                           d1=offsetChans_20;
                           d2=offsetChans_40;
                           chis= 'catch target visible';
                           linestyle='--';
                   end
                   
                   
                   
                   
                   %if restricting the channels for comparison
                   d1=d1(:,selectchans,:);
                   d2=d2(:,selectchans,:);
                   
                  
                   
                   corr_time = zeros(size(d1,1), size(d1,3));
%                    p_time = zeros(size(d1,1), size(d1,3));
                   % per ppant, calculate corr over time
                   for ippant = 1:size(d1,1)
                      for itime= 1:size(d1,3)
                       
%                              scatter(squeeze(d1(ippant,:,itime))', squeeze(d2(ippant,:,itime))');
                          [r,p]= corr(squeeze(d1(ippant,:,itime))', squeeze(d2(ippant,:,itime))');
                       corr_time(ippant,itime)=r;
%                        p_time(ippant,itime)=r;
                      end
                                             
                   end
                   
                    %%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%
                   %%%%% %%%%%%%%% %%%% Interp catch section
                   if interpcatch==1 && iPFIdir==1
                       %may need to interpolate at trial level. see what this looks like.
                       %time vector
                       tvector=tgrm-3;
                       points= 1:length(tvector);
                       %remove bad section from trace.
                       %we used a one second sliding window, so:
                      
                           badsec = dsearchn(tvector', [-1.45 -.36]');
%                       badsec = dsearchn(tvector', [-1.7 -.36]');
                       
                       tmp=points;
                       tmp(badsec(1):badsec(2))=[];
                       %input x vector with points missing.
                       xIN=tmp;
                       plotd=corr_time;
                       outg=[];
                           for ippant=1:size(d1,1)
                               
                                   
                                   reald= squeeze(plotd(ippant,:));
                                   
                                   %remove 'bad' points
                                   realrem = [reald(1, 1:(badsec(1)-1)), reald(1,badsec(2)+1:end)];
                                   
                                   %interp
                                   querypoints = badsec(1):badsec(2);
                                   y=interp1(xIN, realrem, querypoints, 'spline');
                                   
                                   %now replace the data with the interpolated values
                                   
                                   reald2=reald;
                                   reald2(1,badsec(1):badsec(2))=y;
                                   
                                   %plot for comparison.
                                   %       clf;
                                   %       plot(reald); hold on; plot(reald2)
                                   
                                   
                                   %      replace for plotting
                                   
                                   outg(ippant,:)=reald2;
                               
                           end
                           corr_time=outg;
                   end
                   
                   
                   
                 
                   %%%%% %%%%%%%%% %%%%%%%%% %%%%%%%%% %%%%
                   %%%%% %%%% end interp catch section
                     
                   
                   %within subj error bars
                   
                    
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = corr_time;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
                   
                   corr_time=NEWdata;
                   
                   
%                    subplot(1,2,iPFIdir)
                   mP= squeeze(mean(corr_time,1));
                   
                   
                   stP = std(corr_time)/sqrt(size(corr_time,1));
                   
                   
                   st=shadedErrorBar(tgrm-3, mP, stP,[],1);
                   st.mainLine.Color= colis;
                   st.mainLine.LineStyle=linestyle;
                   st.mainLine.LineWidth=3;
                   st.patch.FaceColor=colis;
                   st.edge(1).Color=colis;
                   st.edge(2).Color=colis;
                   
                   leg(iPFIdir)=st.mainLine;
                   
%                    ylim([-.1 .3])
axis tight
                  xlim([-2.5 2.5])
                   
                   
                   
                   ylabel(['1f vs. 2f correlation [\itr\rm ]'])
%                    xlabel(['Time from ' chis ])
                   xlabel(['Time from catch report'])
                   set(gca, 'fontsize', 25)
                   set(gcf, 'color', 'w')
                  
                   ttestdata(iPFIdir,:,:) =x;
                   
                   
               end
               
               %plot sig
               pvals = zeros(1,size(ttestdata,3));
               for itime=1:size(ttestdata,3)
                   
                   [h, pvals(itime) ]= ttest(ttestdata(1,:,itime), ttestdata(2,:,itime)) ;
                                  
               end
               q=fdr(pvals, .05);
               sigspots = find(pvals<q);
               
               
               for itime = 1:length(sigspots)
                   
                   rtime = sigspots(itime);
                   hold on
                   
                   %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                   plot(tgrm(rtime)-3, .45, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
               end
%                legend([leg(1) leg(2)], {'target invisible' 'target visible'})
%%
               ylim([.24 .55])
% axis tight
shg
             %%
             cd('Figures')
             %%
             print('-dpng', 'SpatialCorrelation BG freqs Catch')
         end
    
         
         %and topos
         
    
     
end
   




if job.erpimageCATCHusingSNR==1
    getelocs
    window=[-3 3];
    
    srate=250;
    
    epochdur = sum(abs(window))*srate;
    
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
    
    onsetc = ceil(epochdur)/2;
    % peakfreqsare=[20,40]; %hz
    %timing
    tt = 0:1/srate:60;
    
    %%
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.Fpass= [0 50];
    param_spctrm.trialave=0;
    
    param_spcgrm.tapers = [1 1];
    param_spcgrm.Fs= [250];
    param_spcgrm.Fpass= [0 50];
    param_spcgrm.trialave=0;
    movingwin=[1,.15];
    
    rmvbase=0;
    
    %%
    for ifol = allppants
        
        
        for hzis=1:2
            switch hzis
                case 1
                    usehz=20;
                case 2
                    usehz=40;
            end
            
        
            for alignment=1:2
                
                cd(basefol)
                cd(num2str(ifol))
                
                load(['ppant_Catch_Epoched'])
%                 load('ProposedCATCHonsettiralindextorej')
                for itimezero = 1:2
                    if itimezero==1
                        
                        
                        if alignment==1
                            datatouse = ppant_SNREEG_catchramponset;
                            RTstouse = firstBPwithincatch_RTs;
                            BPstouse = catchonsetBPs;
                            ctype = 'onset';
                        else
                            
                            datatouse = ppant_SNREEG_disapBPwithincatch;
                            RTstouse = firstBPwithincatch_RTs;
                            BPstouse = withincatchonsetBPs;
                            ctype = 'onset, BP aligned';
                        end
                        bsrem = [-3 -1]; %seconds
                        
                    else
                        
                        if alignment==1
                            datatouse = ppant_SNREEG_catchrampoffset;
                            RTstouse = firstBPaftercatch_RTs;
                            BPstouse = catchoffsetBPs;
                            ctype = 'offset';
                        else
                            datatouse = ppant_SNREEG_reapBPaftercatch;
                            RTstouse = firstBPaftercatch_RTs;
                            BPstouse = postcatchoffsetBPs;
                            ctype = 'offset, BP aligned';
                            
                            
                        end
                        
                        
                        bsrem = [1 3]; %seconds
                        
                    end
                    
                    
                    
                    %plot the topo for sanity check: (pre disap)
                    windcheck= [-3 -0.1];
                    tidx=dsearchn(timeid', [windcheck]');
                    
                    %reduce size.
                    topod= squeeze(datatouse(:,:,tidx(1):tidx(2)));
                    %calc statSNR per channel, over trials
                    
                    
                    %%
                    snr20=zeros(64,1);
                    
                    for ichan= 1:64
                        
                        [s,f]=mtspectrumc(squeeze(topod(:,ichan,:))', param_spctrm)    ;
                        
                        %comp SNR;
                        kernelw = [-.25 -.25 0 0 1 0  0 -.25 -.25];
                        
                        sNOW=zeros(size(s));
                        for tmptr=1:size(s,2)
                            
                            sNOW(:,tmptr)=conv(log(s(:,tmptr)), kernelw, 'same');
                        end
                        
                        %store just stim freq.
                        [~, hzid]= min(abs(f-usehz));
                        
                        snr20(ichan)=squeeze(nanmean(sNOW(hzid,:),2));
                    end
                    %%
                    %             [~,usechan]=max(snr20);
                    usechan=30; %Oz for all
                    
                    
                    %%
                    figure(1)
                    clf
                    subplot(4,2,1:2)
                    topoplot(snr20, elocs(1:64), 'emarker2', {usechan, 'o', 'w'});
                    c=colorbar;
                    caxis([0 2])
                    ylabel(c, 'SNR')
                    title({[num2str(usehz) 'Hz SNR '];['-3000:-100ms']})
                    
                    %%
                    
                    %plot BP for comparison
                    subplot(4,2, [3 5])
                    %sort trial index by longest RT first.
                    %sort by duration?
                    
                    %                 [sortedRTs, cid] = sort(RTstouse, 'descend');
                    %
                    %sort by accum Targ gone.
                    %sort by accum PFI in window after onset.
                    if itimezero==1
                        accumPeriod = sum(BPstouse(:,180:end),2);
                    else
                        accumPeriod = sum(BPstouse(:,1:180),2);
                    end
                    
                    
                    
                    
                    [checkBP, cid] = sort(accumPeriod, 'descend');
                    
                    %% MOVE nan to bottom of order (not top)
                    nanind= find(isnan(checkBP));
                    %new start point
                    cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                    %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                    
                    %rearrange.
                    BPstouse=BPstouse(cid,:);
                    datais=squeeze(datatouse(cid,usechan,:));
                    
                    
                    
                    
                    imagesc(-3:1/60:3,1:24,BPstouse)
                    c=colorbar;
                    ylabel(c, 'buttons pressed')
                    set(gca, 'ytick', 1:24, 'yticklabel', cid)
                    ylabel('Trialind')
                    xlabel('time [sec]')
                    hold on
                    plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                    title({['sorted BP data for catch ' ctype ','];['ppant' num2str(ifol)]})
                    set(gca, 'fontsize', 10)
                    %% show mean over time.
                    subplot(4,2,7)
                    plot([-3:1/60:3], nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                    %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                    set(gca, 'fontsize', 15)
                    hold on;
                    ylabel('nanmean BP ')
                    %
                    xlabel('Time secs')
                    axis tight
                    plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                    
                    
                    set(gca, 'fontsize', 10)
                    %%
                    % now plot spectrogram SNR.
                    
                    
                    %
                    %             %rmvbaseline from EEG.
                    bsdata = zeros(size(datais));
                    
                    for itrial = 1:size(datais,1)
                        td = detrend(datais(itrial,:), 'linear');
                        tdrm= mean(td(1,1:250));
                        rmb= repmat(tdrm, [1 length(td)]);
                        bsdata(itrial,:) = td-rmb;
                    end
                    datais=bsdata;
                    %%
                    [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                    %%
                    %conv SNR
                    snr_sgrm =zeros(size(sgrm));
                    
                    kernelw= [ -1/6 -1/6 -1/6 0 0 1 0 0 -1/6 -1/6 -1/6];
                    
                    for itrial=1:size(sgrm,3)
                        %compute SNR
                        tmps= squeeze(sgrm(:,:,itrial));
                        
                        for itime= 1:size(tmps,1)
                            checkput = conv(log(tmps(itime,:)), kernelw,'same');
                            if ~isreal(checkput)
                                snr_sgrm(itime,:,itrial)= nan(1, size(tmps,2));
                            else
                                snr_sgrm(itime,:,itrial)= conv(log(tmps(itime,:)), kernelw,'same');
                            end
                        end
                    end
                    
                    %store just stim freq.
                    [~, hzid]= min(abs(fgrm-usehz));
                    %%
                    %reduce size.
                    snrgrm20=squeeze(snr_sgrm(:,hzid,:))';
                    %
                    
                    tbase = tgrm-3;
                    tidFREEZE = dsearchn(tbase', [-.15 .35]');
                    
                    if rmvbase==1
                        %rmv baseline ?
                        acrSNRsort_rmv=zeros(size(snrgrm20));
                        for itrial=1:size(snrgrm20,1);
                            
                            tmp=snrgrm20(itrial,:);
                            %which baseline?
                            tidx= dsearchn(tbase', [bsrem]');
                            
                            %rmvbase
                            bs= mean(tmp(1,tidx(1):tidx(2)));
                            bs=repmat(bs, 1, length(tmp));
                            acrSNRsort_rmv(itrial,:)= tmp-bs;
                            %                         acrSNRsort_rmv(itrial,:)= tmp./bs;
                            
                        end
                        snrgrm20=acrSNRsort_rmv;
                        
                        
                        
                        
                        
                        
                    end
                    %
                    
                    % smooth across trials
                    sm_snrgrm20=zeros(size(snrgrm20));
                    for itime=1:size(snrgrm20,2)
                        sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),3, 'moving');
                    end
                    
                    subplot(4,2,[4 6]);
                    imagesc(tgrm-3,  1:size(snr_sgrm,3), sm_snrgrm20);
%                     imagesc(tgrm-3,  1:size(snr_sgrm,3), snrgrm20);
                    c=colorbar;
                    ylabel(c, 'SNR')
                    xlabel('Time secs')
                    caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
                    %                     caxis([-10 10])
                    set(gca, 'ytick', 1:24, 'yticklabel', cid)
                    hold on
                    %%
                    plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                    title({['sorted ' num2str(usehz) 'Hz SNR data (smoothed) for catch ' ctype ','];['ppant' num2str(ifol)]})
                    set(gca, 'fontsize', 10)
                    
                    %% show mean over time.
                    subplot(4,2,8)
                    plot(tbase, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)
                    %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                    set(gca, 'fontsize', 15)
                    hold on;
                    ylabel('nanmean SNR (basesub) ')
                    
                    axis tight
                    plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                    %
                    xlabel('Time secs')
                    set(gcf,'color', 'w');
                    
                    %         ylim([0 20])
                    
                    %%
                    
%                     print('-dpng', ['figure_catch_' ctype '_summary_BPandSSVEP_' num2str(usehz) '.png']);
                    
                    %sanitychecks
%                     figure(2)
%                     clf
%                         for ipl=1:24
%                             subplot(6,4,ipl), 
%                             plot(tgrm-3, snrgrm20(ipl,:));
%                             title([num2str(cid(ipl))]);
%                             ylim([-20 20]); hold on;
%                             plot([0 0], ylim, ['k:']);
%                         end
                    %%
                    switch itimezero
                        case 1 %store for across ppant plots:
                            Ppant_onsetBP=BPstouse;
                            Ppant_onsetSNR=snrgrm20; %sorted.
                            Ppant_onsetRTs=RTstouse;
                            Ppant_onsetTOPO=snr20;
                        case 2
                            
                            Ppant_offsetBP=BPstouse;
                            Ppant_offsetSNR=snrgrm20; %always sorted in descending order of PFI.
                            Ppant_offsetRTs=RTstouse;
                            Ppant_offsetTOPO=snr20;
                    end
                end
                
                
                switch hzis
                    case 1
                        savename='Catchperformance_withSNR_20';
                    case 2
                        savename='Catchperformance_withSNR_40';
                end
                
                if alignment==2
                    savename = [savename ',BPaligned'];
                end
                save(savename,...
                    'Ppant_onsetBP','Ppant_offsetBP',...
                    'Ppant_onsetSNR', 'Ppant_offsetSNR', ...
                    'Ppant_onsetRTs','Ppant_offsetRTs',...
                    'Ppant_onsetTOPO', 'Ppant_offsetTOPO', 'tgrm')
                
                
                
                
            end
            
        end
        
        
        
    end
end


if job.concaterpdataacrossppants==1
    %%
    removeFrozentrials=0; % can remove the 'bad' trials per ppant.
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=1:2
        
        if ihz==1
            loadname='Catchperformance_withSNR_20';
        else
            loadname='Catchperformance_withSNR_40';
        end
        for ialignment=1:2
            if ialignment==2
                loadname=[loadname ',BPaligned'];
            end
            storeacrossPpant_onsetBP=[];
            storeacrossPpant_onsetSNR=[];
            storeacrossPpant_onsetRTs=[];
            storeacrossPpant_onsetTOPO=[];
            storeacrossPpant_offsetBP=[];
            storeacrossPpant_offsetSNR=[];
            storeacrossPpant_offsetRTs=[];
            storeacrossPpant_offsetTOPO=[];
            
            icounter=1;
            for ippant = allppants
                cd(basefol)
                cd(num2str(ippant));
                %onset types
                
                load(loadname)
                
                
                %alsoload index of frozen frames:
                keeptrials=[];
                if ihz==2 && removeFrozentrials==1
%                 load('ProposedCATCHonsettrialindextorej')
                keeptrials=1:24;
                % note only for catch ONSETS (physical).
                keeptrials(descendingtrialREJ)=[];                                
                Ppant_onsetSNR=Ppant_onsetSNR(keeptrials,:);
                Ppant_onsetBP=Ppant_onsetBP(keeptrials,:);
                end
                
                
                
                if any(isnan(Ppant_onsetBP(:,1)))
                    lc= min(find(isnan(Ppant_onsetBP(:,1))));
                    lc=lc-1;
                else
                    lc=size(Ppant_onsetBP,1);
                end
                
                % we need to resample the BP and SNR data, to equate across
                % trial types    (ignoring nan)
                % currently at 100 fs.
                %%
                Ppant_onsetBP= resample(Ppant_onsetBP(1:lc,:),100,lc);
                Ppant_onsetSNR= resample(Ppant_onsetSNR(1:lc,:),100,lc);
                
                 
                %offsets too
                if any(isnan(Ppant_offsetBP(:,1)));
                    lc= min(find(isnan(Ppant_offsetBP(:,1))));
                    lc=lc-1;
                else
                    lc=size(Ppant_offsetBP,1);
                end
                
                % we need to resample the BP and SNR data, to equate across
                % trial types    (ignoring nan)
                % currently at 100 fs.
                %%
                Ppant_offsetBP= resample(Ppant_offsetBP(1:lc,:),100,lc);
                Ppant_offsetSNR= resample(Ppant_offsetSNR(1:lc,:),100,lc);
                
               % also apply participant level smoothing.
                if ppantsmoothing==1
                    %smooth across trials, per ppant for SNR
                    
                    
                    % note, need to adjust for edge artefacts in smoothig
                    % process, repmat the end trials so smoothing doesnt
                    % contaminate.
                    
                    
                    for ionoff=1:4
                        
                        switch ionoff
                            case 1
                                pData = Ppant_onsetSNR;
                            case 2
                                pData = Ppant_offsetSNR;
                            case 3
                                pData = Ppant_onsetBP;
                            case 4
                                pData = Ppant_offsetBP;
                        end
                        
                        
                        %new data with expanded width.
                        pdataN = zeros(size(pData,1)+32, size(pData,2));
                        %top of image
                        pdataN(1:16,:) = repmat(pData(1,:), [16,1]);
                        %and bottom;
                        pdataN(end-15:end,:) = repmat(pData(end,:), [16,1]);
                    
                        %fill centre
                        pdataN(17:size(pData,1)+16,:) = pData;
                        
                    
                    %now safe to smooth
                    
                    pdataOUT= zeros(size(pdataN));
                    for itime=1:size(pdataN,2)
                        pdataOUT(:,itime)= smooth(pdataN(:,itime),15);
                    end
                    
                    %now collapse to original size.
                    pdataOUT = pdataOUT(17:end-16,:);
                    
                    switch ionoff
                        case 1
                            Ppant_onsetSNR = pdataOUT;
                        case 2
                            Ppant_offsetSNR = pdataOUT;
                            case 3
                            Ppant_onsetBP = pdataOUT;
                        case 4
                            Ppant_offsetBP = pdataOUT;
                            
                    end
                    
                    end
                
                end
                
                
                
                
                storeacrossPpant_onsetBP(icounter,:,:)=Ppant_onsetBP;
                storeacrossPpant_onsetSNR(icounter,:,:)=Ppant_onsetSNR;
                storeacrossPpant_onsetRTs(icounter,:)=Ppant_onsetRTs;
                storeacrossPpant_onsetTOPO(icounter,:)=Ppant_onsetTOPO;
               
                
                 storeacrossPpant_offsetBP(icounter,:,:)=Ppant_offsetBP;
                storeacrossPpant_offsetSNR(icounter,:,:)=Ppant_offsetSNR;
                storeacrossPpant_offsetRTs(icounter,:)=Ppant_offsetRTs;
                storeacrossPpant_offsetTOPO(icounter,:)=Ppant_offsetTOPO;
                
                
                
                
                
                icounter=icounter+1;
            end
            
            %save appropriately
            cd(basefol)
            cd('newplots-MD')
            
            
            switch ihz
                case 1
                    savename='GFX_Catchperformance_withSNR_20';
                case 2
                    savename='GFX_Catchperformance_withSNR_40';
            end
            
            if ialignment==2
                savename = [savename ',BPaligned'];
            end
            save(savename,...
                'storeacrossPpant_onsetBP','storeacrossPpant_offsetBP',...
                'storeacrossPpant_onsetSNR', 'storeacrossPpant_offsetSNR', ...
                'storeacrossPpant_onsetRTs','storeacrossPpant_offsetRTs',...
                'storeacrossPpant_onsetTOPO', 'storeacrossPpant_offsetTOPO', 'tgrm')
            
            
            
        end
    end
    
    
end


if job.erpimageacrossppants==1
    %% now plot across ppants.
    %reshape manually.
    
    %%
    getelocs
    rmvbase=0;
    for alignment=1%1:2 % 1 for physical change at time zero, 2 for BP noticing catch change at time zero.
        clf
        
        for hzis=1:2
            cd(basefol)
            cd('newplots-MD')
            switch hzis
                case 1
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_20')
                    else
                        load('GFX_Catchperformance_withSNR_20,BPaligned')
                    end
                    
                    usehz=20;
                    
                case 2
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_40')
                    else
                        load('GFX_Catchperformance_withSNR_40,BPaligned')
                    end
                    usehz=40;
                    
            end
            tbase = tgrm-3;
            %%
            for itimezero=1%:2
                
                
                switch itimezero
                    case 1
                        useBP=storeacrossPpant_onsetBP;
                        useSNR=storeacrossPpant_onsetSNR;
                        useRTs=storeacrossPpant_onsetRTs;
                        useTOPO=storeacrossPpant_onsetTOPO;
                        
                        if alignment==1
                            ctype= 'onset';
                            chtype='catch onset';
                        else
                            ctype= 'onset, BPaligned';
                            chtype='reporting catch onset';
                        end
                        bsrem = [-3 -2]; %seconds
                    case 2
                        
                        useBP=storeacrossPpant_offsetBP;
                        useSNR=storeacrossPpant_offsetSNR;
                        useRTs=storeacrossPpant_offsetRTs;
                        useTOPO=storeacrossPpant_offsetTOPO;
                        
                        if alignment==1
                            ctype= 'offset';
                            chtype= 'catch offset';
                        else
                            ctype= 'offset, BPaligned';
                            chtype='reporting catch offset';
                        end
                        
                        bsrem = [2 3]; %seconds
                        
                end
                
                
                
                %                 for ip= 1:length(allppants)
                %                     trialind = [1:24] + 24*(ip-1);
                %                     acrBP(trialind,:,:) = useBP(ip,:,:);
                %                     acrSNR(trialind,:,:)=    useSNR(ip,:,:);
                %                     acrRTs(trialind)=    useRTs(ip,:);
                %
                %                 end
                
                %take mean across ppants for erp images..
                acrBP= squeeze(nanmean(useBP,1));
                
                acrSNR=squeeze(nanmean(useSNR,1));
                acrTOPO=squeeze(nanmean(useTOPO,1));
                
                
                
                %sort across all
                %sort trial index by longest RT first.
                %%
                figure(1)
%                 clf
%% not topoplotting anymore.                
%                 usechan=30;
%                 subplot(2,2,1:2)
%                 topoplot(acrTOPO, elocs(1:64), 'emarker2', {usechan, 'o', 'w'});
%                 c=colorbar;
%                 caxis([0 15])
%                 ylabel(c, {['10log10(SNR)'];['-3000:-100ms']})
%                 hold on
%                 title({[num2str(usehz) ' Hz SSVEP']})
%                 set(gca, 'fontsize', 20)
                %%
                subplot(3,1,1)
                %         [sortedRTs, cid] = sort(acrRTs, 'descend');
                
                %
                imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
                %
                title({['Buttons Pressed'] }, 'fontsize', 25)
                c=colorbar
                ylabel(c, 'Total')
                %                 ylabel('resampled catch trials')
                xlim([-2.5 2.5])
                set(gca, 'fontsize', 25, 'ytick', [])
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                plot([0 0] , ylim, ['k:'], 'linewidth', 2)
                %                 subplot(4,2,7)
                
                ppantMeanBP= squeeze(mean(useBP,2));
                %adjust standard error as per COusineau(2005)
                %confidence interval for within subj designs.
                % y = x - mXsub + mXGroup,
                x = ppantMeanBP;
                
                mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                
                %for each observation, subjtract the subj average, add
                %the group average.
                NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                
                %compute new stErr %which version?
                %             stE = std(NEWdata)/sqrt(size(x,1));
                %                 shadedErrorBar([-3:1/60:3], mean(ppantMeanBP,1),stE)
                %                 set(gca, 'fontsize', 15)
                %                 hold on;
                                ylabel({['Normalized'];['trial count']})
                %                 plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                %%
                xlabel(['Time from ' chtype])
                
                %SNR
                subplot(3, 1, hzis+1);
                
                %reorder SNR
                acrSNRsort=acrSNR;
                
                %%
                %                 %smooth across trials
                sm_snrgrm20=zeros(size(acrSNRsort));
                
                for itime=1:size(acrSNRsort,2)
                    sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
                end
                
                
                %     imagesc(acrSNR);
                imagesc(tgrm-3, 1:size(acrSNRsort,1), acrSNRsort);
                c=colorbar;
                ylabel(c, {['SNR(dB/Hz)'];['baseline corrected']})
                %
                caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])
%                 caxis([-5 5])
                        caxis([0 2])
                title([num2str(usehz) 'Hz dyn-SSVEP SNR'])
                %             set(gca, 'ytick', 1:24, 'yticklabel', round(sortedRTs./60,2))
                hold on
                plot([0 0] , ylim, ['k:'], 'linewidth', 3)
                plot([0 0] , ylim, ['w:'], 'linewidth', 2)
                xlim([-2.5 2.5])
                set(gca, 'fontsize', 25, 'ytick', [])
                %                 subplot(4,2,8)
                %
                %plot across ppant trace
                ppantMeanSNR= squeeze(mean(useSNR,2));
                
                %adjust standard error as per COusineau(2005)
                %confidence interval for within subj designs.
                % y = x - mXsub + mXGroup,
                x = ppantMeanSNR;
                
                mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                
                %for each observation, subjtract the subj average, add
                %the group average.
                NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                
                %             %compute new stErr %which version?
                %             stE = std(NEWdata)/sqrt(size(x,1));
                %
                %                 shadedErrorBar(tgrm-3, mean(ppantMeanSNR,1),stE,'k',[])
                %                 set(gca, 'fontsize', 15)
                %                 hold on;
                %                 %%
                %
                %                 ylabel('mean SNR')
                %                 axis tight
                %                 plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                xlabel(['Time from ' chtype])
                ylabel({['Normalized'];['trial count']})
                set(gca, 'fontsize', 25)
                cd(basefol)
                cd('newplots-MD')
                cd('Figures')
                set(gcf, 'color', 'w')
                shg
                %%
                if itimezero==1 && usehz==40 && alignment==1
                    yt=ylim;
                    xV = [-.5 -.5 .5 .5]; %x coordinates of vertices.
                    yV = [yt(1) yt(2) yt(2) yt(1)]; %y coordinates of vertices.
                    pch=patch(xV,yV,[.5 .5 .5]);
                    pch.FaceAlpha =(.75);
                    shg
                end
                %%
                
            end
        end
        %%
           print('-dpng', ['Catch SNR summary,'  chtype '.png'])
             
    end
end

if job.BPandSSVEPtimecourseacrossppants==1;
    %%
    getelocs
    rmvbase=0;
    checksigON=0;
    
    for alignment=2%:2%1:2 % 1 for physical change at time zero, 2 for BP noticing catch change at time zero.
        figure(1);
        clf
        plcount=1;
        legendprint=[];
        for hzis=1:2
            cd(basefol)
            cd('newplots-MD')
            switch hzis
                case 1
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_20')
                    else
                        load('GFX_Catchperformance_withSNR_20,BPaligned')
                    end
                    lint='-';
                    usehz=20;
                    sigheight= 4.1; %where on figure to place '*' if significant.
                case 2
                    if alignment==1
                        load('GFX_Catchperformance_withSNR_40')
                    else
                        load('GFX_Catchperformance_withSNR_40,BPaligned')
                    end
                    usehz=40;
                    lint=':';
                    sigheight= 3.9; %where on figure to place '*' if significant.
            end
            tbase = tgrm-3;
            %%
            %             clf
            ttestdata=[];
            for itimezero=1:2%1%:2
                
                
                switch itimezero
                    case 1
                        useBP=storeacrossPpant_onsetBP;
                        useSNR=storeacrossPpant_onsetSNR;
                        useRTs=storeacrossPpant_onsetRTs;
                        useTOPO=storeacrossPpant_onsetTOPO;
                        
                        if alignment==1
                            ctype= 'onset';
                            chtype='catch start-end';
                            xlabelis = {['Time from catch onset [s]']};
                        else
                            ctype= 'onset, BPaligned';
                            xlabelis = {['Time from reporting'];['catch onset [s]']};
                            
                        end
                        bsrem = [-3 -2]; %seconds
                        col=[0 .5 0]; %dark green.
                    case 2
                        
                        useBP=storeacrossPpant_offsetBP;
                        useSNR=storeacrossPpant_offsetSNR;
                        useRTs=storeacrossPpant_offsetRTs;
                        useTOPO=storeacrossPpant_offsetTOPO;
                        
                        if alignment==1
                            ctype= 'offset';
                            xlabelis = {['Time from catch offset [s]']};
                            
                        else
                            ctype= 'offset, BPaligned';
                            xlabelis = {['Time from reporting'];['catch offset [s]']};
                            chtype='catch report';
                        end
                        
                        bsrem = [2 3]; %seconds
                        col='r';
                        
                        
                end
                
                
                %%
                figure(1)
                for ip=1:2
                    switch ip
                        case 1
                            used=useBP;
                            col='b';
                            lint='-';
                            ylimsare = [0 2.5];
                            basep=1;
                        case 2
                            used=useSNR;
                            basep=3;
                            if hzis==1
                                col='r';
                                lgc=1;
                            else
                                col=[0, .5, 0];
                                lgc=2;
                            end
                            lint=':';
                            ylimsare= [ 1 1.9];
                    end
                    placeme = basep+ (1*itimezero-1);
                    
                    subplot(2,2,placeme)
                    
                    
                    hold on
                    
                    ppantMeanSNR = squeeze(mean(used,2));
                    %adjust standard error as per COusineau(2005)
                    %confidence interval for within subj designs.
                    % y = x - mXsub + mXGroup,
                    
                    x = ppantMeanSNR;
                    if ip~=1
                    mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
                    mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
                    
                    %for each observation, subjtract the subj average, add
                    %the group average.
                    NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
                    
                    else
                        NEWdata=ppantMeanSNR;
                    end
                    
                    %compute new stErr %which version?
                    stE = std(NEWdata)/sqrt(size(x,1));
                    
                    if ip==1
                        sh=shadedErrorBar(-3:1/60:3, mean(ppantMeanSNR,1),stE,[lint],[1]);
                        ylabel(['Buttons pressed'])
                    else
                        sh=shadedErrorBar(tgrm-3, mean(ppantMeanSNR,1),stE,[lint],[1]);
                        ylabel({['log(SNR)']})
                    end
                    sh.mainLine.LineWidth=3;
                    
                    sh.mainLine.Color = col;
                    sh.patch.FaceColor = col;
                    sh.edge(1).Color = col;
                    sh.edge(2).Color = col;
                    
                    %                 if hzis==2
                    %                     sh.mainLine.Color= 'k';
                    %                 end
                    set(gca, 'fontsize', 15)
                    hold on;
                    %                 %%
                    %
                    
                    
                    
                    
                    %             title({[num2str(usehz) ' Hz SSVEP']})
                    xlabel(xlabelis)
                    %             xlabel('Time from perceptual report')
                    set(gca, 'fontsize', 25)
                    xlim([-3 3])
                    ylim( [ylimsare])
                    
                    
                    
                    
                    
                    set(gcf, 'color', 'w')
                    
                    if ip==2
                        %             ttestdata(plcount,:,:) = ppantMeanSNR;
                        legendprint(lgc)=sh.mainLine;
                    end
                    
                    plcount=plcount+1;
                    
                    
                    %place legend
                    if hzis==2 && ip==2
                        if itimezero==1
                            lg=legend([legendprint(1) legendprint(2)], {'20 Hz', '40 Hz'});
                        end
                        set(lg, 'location', 'NorthEast')
                        
                    end
                    
                end
                
                        %plot patch:
        if itimezero==1 && alignment==1 && hzis==2
            xV = [-.5 -.5 .5 .5]; %x coordinates of vertices.
            yV = [-2 1.5 1.5 -2]; %y coordinates of vertices.
            pch=patch(xV,yV,[.5 .5 .5]);
            pch.FaceAlpha= 1;
        end
                
            end
        
        %%
        if checksigON==1
        %check for sig
        pvals=zeros(1,size(ttestdata,3));
        tvals=zeros(1,size(ttestdata,3));
        for itime = 1:size(ttestdata,3)
            
            try [h,pvals(itime),~,stat]=ttest(ttestdata(1,:,itime), ttestdata(2,:,itime));
                shuffType=1;
            catch
                [h,pvals(itime),~,stat]=ttest(ttestdata(plcount-1,:,itime)); %compares to zero.
                shuffType=2; %whether or not to skip the non-parametric test for sig.
            end
            
            tvals(itime)= stat.tstat;
        end
        %%
        %             q=fdr(pvals,.05);
        
            sigs=find(pvals<.05);
            
            %perform cluster based correction.
            if length(sigs)>2
                % find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                %grab largest
                %                 ignore bad points.
                if hzis==2 &&itimezero==1
                    pvals(11:22)=nan;
                end
                
                % find biggest cluster:
                %finds adjacent time points
                sigs = find(pvals<.05);
                
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
                if itimezero==1 && hzis==2 && alignment==2
                    maxClust=3;
                end
                %%
                for icl=1%:size(clusterSTandEND,1)
                    
                    %start and end are now:
                    STC=sigs(clusterSTandEND(maxClust,1));
                    ENDC=sigs(clusterSTandEND(maxClust,2)+1);
                    checktimes =STC:ENDC;
                    observedCV = sum(abs(tvals(checktimes)));
                    % now shuffle condition labels to see if this cluster is
                    % sig (compared to chance).
                    sumTestStatsShuff = zeros(1,2000);
                    
                    for irand = 1:2000
                        %testing the null that it isn't mismatched - matched at time 2
                        % which creates a diff. so select from either!
                        shD= zeros(2,length(checktimes),size(ttestdata,2));
                        
                        %change shuffle parameters based on test of
                        %interest. (ie between conditions, or temporal
                        %null).
                        if shuffType==1 %null is that no condition differences.
                            for ipartition = 1:2
                                for ippant = 1:21
                                    for itime=1:length(checktimes)
                                        
                                        if mod(randi(100),2)==0 %if random even number
                                            pdata = ttestdata(1,randi(21), checktimes(itime)); %select both chans
                                        else
                                            pdata = ttestdata(2,randi(21), checktimes(itime));
                                        end
                                        
                                        shD(ipartition,itime,ippant) = pdata;
                                    end
                                end
                            end
                        else %null is that there are no temporal coincident sig values.
                            for ipartition = 1:2
                                for ippant = 1:21
                                    for itime=1:length(checktimes)
                                        
                                        %take random timepoint.
                                        pdata = ttestdata(1,ippant, randi(size(ttestdata,3)));
                                        
                                        
                                        shD(ipartition,itime,ippant) = pdata;
                                    end
                                end
                            end
                        end
                        %now compute difference between out hypothetical topoplots,
                        % and test for sig, checking the accumulated test statistic at our
                        % times of interest
                        tvalspertimepoint = zeros(1,length(checktimes));
                        
                        testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                        
                        for itest = 1:length(checktimes) %test each time point
                            
                            [~, p, ~,stat]= ttest(testdata(itest,:));
                            
                            tvalspertimepoint(1,itest) = stat.tstat;
                        end
                        
                        sumTestStatsShuff(1,irand) = sum(abs(tvalspertimepoint));
                    end %repeat nshuff times
                    
                    
                    %is the observed greater than CV?
                    % plot histogram:
                    figure(2);
                    if plcount==2
                        clf
                    end
                    subplot(2,1, plcount-1)
                    H=histogram(abs(sort(sumTestStatsShuff)));
                    % fit CDF
                    cdf= cumsum(H.Data)/ sum(H.Data);
                    %the X values (actual CV) corresponding to .01
                    [~,cv05uncorr] = (min(abs(cdf-.95)));
                    [~,cv01uncorr] = (min(abs(cdf-.99)));
                    [~,cv001uncorr] = (min(abs(cdf-.999)));
                    hold on
                    pCV=plot([observedCV observedCV], ylim, ['r-']);
                    
                    p05=plot([H.Data(cv05uncorr) H.Data(cv05uncorr)], ylim, ['k:']);
                    plot([H.Data(cv01uncorr) H.Data(cv01uncorr)], ylim, ['k:']);
                    plot([H.Data(cv001uncorr) H.Data(cv001uncorr)], ylim, ['k:']);
                    legend([pCV p05], {['observed'] ['p01'] })
                    
                    %%
                    if observedCV>H.Data(cv05uncorr)
                        title(['sum tvals = ' num2str(observedCV)]);
                        %              title('Spatial Cluster  Significant!')
                        for itime=checktimes
                            figure(1);
                            hold on
                            plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                        end
                    end
                    
                    %
                end
            end
        end
        %%
        
        
        hold on
%         plot(xlim, [0 0], ['k:']);
%         plot([0 0] , ylim, ['k:'], 'linewidth', 1)

        
        shg
        
        hold on
        
        %%
        
    end
   %%
    cd('Figures')
    print('-dpng', ['Catch trace Bground SSVEP summary, during ' num2str(chtype) '.png'])
end
end% end