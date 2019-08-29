% s3_D_CreateandApplyPpantRESSfilter
%based on/replaces Mkat_MT_taper from Irene.



clear all
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/MatthewDavidson/Desktop/FromIrenes Mac/Data Experiment copy')
    addpath('/Users/matthewdavidson/Documents/GitHub/SSVEP-PFI_BackgroundTag-master')
end

basefol=pwd;

job.applyRESSperppant_toCATCH=1; 


%sanity check 1 - confirm RESS is not operating on neural noise (vs SNR
%during catch periods)
job.compareSNR_RESSperppant_toCATCH_Notargetspresent=0; %operates on pre and post catch
 
job.concatSanitycheck=0; %across ppant pre vs post catch, stores for JASP on desktop.
% 
job.plotSanitycheckprepostRESSoncatch=0; % compares Pre to post PFI for a single target BP(0->1)
job.performLMEonRESS1 =0;
% %sanity check 2 - confirm RESS is not operating non-uniformly, compare
% pre-RESS SNR during precatch and pre PFI (targets present).
% to post-RESS SNR during precatch and prePFI ( targets present). 
% 
% 
job.construct_and_compareRawtoRESSonBGhz=0;

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

rmvbaseEPOCH=1; %0 for no removal,  1 for trial baseline, 2 for pooled freq baseline (compare across conds).

epochdur = sum(abs(window))*srate;
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;




% which frequencies to analyze?/apply?
        peakfreqsare=[8,13,15,18,20,40];
        peakfreqsare=[16,26,30,36]; % can change to 2*f of TGs
        


%%

if job.applyRESSperppant_toCATCH==1
  
    for ifol = allppants
        
        
        %%
        cd(basefol)
        cd(num2str(ifol))
        
        
        %% load the relevant PFI data.
        load('ppant_Catch_Epoched');
        
        if length(peakfreqsare)==6
            load('RESSfilterdata');
        else
            load('RESSfilterdata_2fTG');
        end
        load('TrialIndicesbyLocationandHz');
        %%
        
        %%%%%% POOL all conditions to increase the quality of the
        %%%%%% covariance matrix.
        dataIN=[];
        
        %%
        
        for id=1:4
            
            switch id
                case 1
                    dataIN=ppant_SNREEG_catchramponset;
                    %                     ress_catchonset=zeros(6,4,24,size(dataIN,3));
                    ress_catchonsetTGs=zeros(6,4,6,size(dataIN,3));
                    ress_catchonsetBGs=zeros(6,4,24,size(dataIN,3));
                case 2
                    dataIN=ppant_SNREEG_catchrampoffset;
                    ress_catchoffsetTGs=zeros(6,4,6,size(dataIN,3));
                    ress_catchoffsetBGs=zeros(6,4,24,size(dataIN,3));
                    %                     ress_catchoffset=zeros(6,4,24,size(dataIN,3));
                case 3
                    dataIN=ppant_SNREEG_disapBPwithincatch;
                    ress_BPcatchonsetTGs=zeros(6,4,6,size(dataIN,3));
                    ress_BPcatchonsetBGs=zeros(6,4,24,size(dataIN,3));
                                        
                case 4
                    dataIN=ppant_SNREEG_reapBPaftercatch;
                    ress_BPcatchoffsetTGs=zeros(6,4,6,size(dataIN,3));
                    ress_BPcatchoffsetBGs=zeros(6,4,24,size(dataIN,3));
                                        
                    
            end
            %since all catch, no durs, 
            durscheck=[];
            searchtrials=1:24;
            
            for ifreq=1:length(peakfreqsare)
                for iloc=1:4
                    % apply the correct spatial filter per trial
                    %%
                    
                    
                    %because we have different configurations on different
                    %trials, we need to apply the correct spatial filter
                    %for each type.
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
                    
                    
                    
                     %reduce size.
                    datast= squeeze(dataIN(trialtypesperTargPresent,:,:));
                    
                    
                    
                       %% check for bad trials (noisy)
%                     %std per trial(average over electrodes)
%                     tmp=[];
%                     if ndims(datast)<3
%                         tmp(1,:,:)= datast;
%                         datast=tmp;
%                     end
%                     datastSD = nanstd(squeeze(nanmean(datast,2)),0,2);
%                     
%                     %remove those with 2.5*std from mean.
%                     trialSD=nanstd(datastSD);
%                     mSD=nanmean(datastSD);
%                     keeptrials=1:size(datast,1);
%                     
%                     %remove the trials with excessive noise.
%                     badtrials=find(datastSD>mSD+2.5*trialSD)';
%                     
                    % also skip the trials which were only transient button
                    % presses. (less than one second).                    
%                     shorttrials = find(durscheck<60);
%                     badtrials = [badtrials, shorttrials];
%                     
%                     % remove these from consideration.
%                     keeptrials(badtrials)=[];
%                     datast=datast(keeptrials,:,:);
%                     
%                     
                    
                    
                    %now we have the correct trials, get the appropriate
                    %filter                    
                    evecs = squeeze(ressEVEC_byHzxLoc(ifreq,iloc,:)); %
                    
                    
                    
                    
                    ress_ts1=zeros(size(datast,1), size(datast,3));
                    for ti=1:size(datast,1)
                        ress_ts1(ti,:) = evecs'*squeeze(datast(ti,:,:));
                    end
                    
                    
                    %%
                    
                    %Save the ress data per type.
                    switch id
                        case 1
                            if ifreq<5
                                    ress_catchonsetTGs(ifreq,iloc,:,:)= ress_ts1;

                            else
                                    ress_catchonsetBGs(ifreq,iloc,:,:)= ress_ts1;

                            end
                        case 2
                            
                            if ifreq<5

                                    ress_catchoffsetTGs(ifreq,iloc,:,:)= ress_ts1;

                            else

                                    ress_catchoffsetBGs(ifreq,iloc,:,:)= ress_ts1;

                            end
                            
                        case 3
                            if ifreq<5
                                    ress_BPcatchonsetTGs(ifreq,iloc,:,:)= ress_ts1;
                                
                            else
                                    ress_BPcatchonsetBGs(ifreq,iloc,:,:)= ress_ts1;
                            
                            end
                            
                        case 4
                             if ifreq<5
                                    ress_BPcatchoffsetTGs(ifreq,iloc,:,:)= ress_ts1;
                            
                             else
                                    ress_BPcatchoffsetBGs(ifreq,iloc,:,:)= ress_ts1;
                            
                            end
                    end
                end
            end
            clearvars dataIN
            
        end
        
        if length(peakfreqsare)==6
            savename='ppant_Catch_Epoched_RESS';
        else
            savename='ppant_Catch_Epoched_RESS_2fTG';
        end
        save(savename, 'ress_catchonsetTGs', 'ress_catchonsetBGs',...
            'ress_catchoffsetTGs', 'ress_catchoffsetBGs',...
            'ress_BPcatchonsetTGs', 'ress_BPcatchonsetBGs',...
            'ress_BPcatchoffsetTGs', 'ress_BPcatchoffsetBGs')
        
    end
    
    
end


%% %% construct participant TARGET present  vs absent (catch)


if job.compareSNR_RESSperppant_toCATCH_Notargetspresent==1
    
    % Compute the average SNR for PFI pre disappearance (actual targets present)
    % Store and compare to post CATCH (actual targets missing - to control
    %for overfitting bias).
     
    epochdur = sum(abs([-3 3]))*srate;
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
 
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.fpass= [0 50];
    param_spctrm.trialave=0;
    
    
    for ifol=allppants
        cd(basefol)
        cd(num2str(ifol))
        %%
        
%         load('ppant_PFI_RESS')
%         load('ppant_Catch_RESS')
        load('ppant_Catch_Epoched')

        % now combining to compare BOTH TG2 and TG2.
        
        load('RESSfilterdata.mat')
        ress1=ressEVEC_byHzxLoc;        
        load('RESSfilterdata_2fTG.mat') %now checking 2nd harmonic also.
        ress2=ressEVEC_byHzxLoc;
        ressEVEC_byHzxLoc=[];
        ressEVEC_byHzxLoc(1:6,:,:)=ress1;
        ressEVEC_byHzxLoc(7:10,:,:)=ress2;
        load('TrialIndicesbyLocationandHz.mat')
        %%
           
    windowsmall=[];
    windowsmall(1,:) = [-2 -.1] ;% short window targets present
    windowsmall(2,:) = [.1 2] ;% short window targets absent
    windowsmall(1,:) = [-3 -.1] ;% short window targets present
    windowsmall(2,:) = [.1 3] ;% short window targets absent
    
        
        
        dtype=2;%
        
        switch dtype
            case 1 %use PFI data (lotsmore trials)
                %
                allTS= ress_catchonset.ressTS;
                markert='x';
                dtypep='CATCH';
                
            case 2 %USE epoched EEG< construct new RESS filters -300:-100ms.
                allTS= ppant_SNREEG_catchramponset;
                markert='x';
                dtypep='CATCH separated';
                
        end
        
        
        peakfreqsare = [8,13,15,18, 20,40, 16,26,30,36];
        
        
        ppantRESSsanitycheck =zeros(2,4, length(peakfreqsare), 205); %window, locations, freqs, hz
        mbar=[];barstore=[];
        %
        for ifreq=1:length(peakfreqsare)
            peakfreq1= peakfreqsare(ifreq);
            %%
            
            %
            if ifreq>6
                usefreq=ifreq-6; %needed for trial indexing
            else
                usefreq=ifreq;
            end
                
            for wind=1:2
                if wind==1
                    col='b';
                    %                 allTS = ressPFI_0_1.ressTS;
                else
                    col='k';
                end
                
                windcheck= windowsmall(wind,:);
                tidx=dsearchn(timeid', [windcheck]');
                
                                
                
                for iloc=1:4
                    
                    %%
                    % Ensure we are only grouping trials withthe same spatial
                    % configuration.
                    
                    if usefreq<5 %&& wind==2%
                        
                        if wind==1
                            switch iloc
                                case 1
                                    trialsthisloc = TopLeftTrialindexbyHz(usefreq,:);
                                case 2
                                    trialsthisloc = TopRightTrialindexbyHz(usefreq,:);
                                case 3
                                    trialsthisloc = BottomLeftTrialindexbyHz(usefreq,:);
                                case 4
                                    trialsthisloc = BottomRightTrialindexbyHz(usefreq,:);
                            end
                            
                        else %  find when missing 
                            switch iloc
                                case 1
                                    trialsthisloc = TopLeftCatchindexbyHz(usefreq,:);
                                case 2
                                    trialsthisloc = TopRightCatchindexbyHz(usefreq,:);
                                case 3
                                    trialsthisloc = BottomLeftCatchindexbyHz(usefreq,:);
                                case 4
                                    trialsthisloc = BottomRightCatchindexbyHz(usefreq,:);
                            end
                            
                            
                        end
                        trialsthisloc(trialsthisloc==0)=[];
                        
                        allTSnow= allTS(trialsthisloc,:,:);
                    else
                        %if 20 or 40hz, no loc separation.
                        allTSnow = allTS;
                        
                    end
                    
                    %reduce for appropriate window.
                    datast=allTSnow(:,:,tidx(1):tidx(2));
                   %multiply by appropriate ress component
                   
                    
                    %now we have the correct trials, get the appropriate
                    %filter                    
                    evecs = squeeze(ressEVEC_byHzxLoc(ifreq,iloc,:)); %
                    
                    
                    
                    ress_ts1=zeros(size(datast,1), size(datast,3));
                    for ti=1:size(datast,1)
                        ress_ts1(ti,:) = evecs'*squeeze(datast(ti,:,:));
                    end
                    
                    
                    %calculate SNR
                    
                    
                    [s,f]=mtspectrumc(ress_ts1', param_spctrm)    ;
                    
                    %TAKE mean
                    sMEAN = squeeze(nanmean(log(s),2));
                    
                    %% comp SNR;
%                     kernelw = [-1/6 -1/6 -1/6   0 0 1  0 0 -1/6  -1/6 -1/6];
                    kernelw = [-1/4 -1/4   0 0 1  0 0   -1/4 -1/4];
                    %                 
%                     %take mean
                    sNOW= conv(sMEAN, kernelw, 'same');
                %
%                 figure(1)
%                 hold on
%                 plot(f,sNOW);   
%                 shg
                    
                    
                    ppantRESSsanitycheck(wind, iloc, ifreq,:)= sNOW;
                    % also collect SNR across trials for SD/variance estimate.
                [~,fis]=min(abs(f-peakfreq1));
                
                snrR=squeeze(nanmean(ppantRESSsanitycheck,2)); %mean over locations.
                mbar(wind,ifreq,:) = nanmean(snrR(fis));
                % 
                end %each location
                
                
                
                
            end
        end
        %%% plot new SNR and save.
        %%        
%         figure(3);

%                 clf
%                 hold on %plot pre and post for comparison.
%                 plot(f, squeeze(snrR(1,:,:)), ['r' markert '-' ])
%                 col='k'
%                 plot(f, squeeze(snrR(2,:,:)), [col markert '-' ])
%                 hold on
%                 xlim([0 45])
%                 
                
                  
%         %%
%         set(gcf,'color', 'w')
%         legend({'Target PRESENT', ['Target ' dtypep]})
%         xlim([ 0 45])
%         ylim([0 50])
%         set(gca, 'fontsize', 20)
%         xlabel('Hz')
%         ylabel('SNR')
%         title('Instances of Catch disappearance only')
%         %     print snr
%         print('-dpng', ['figure_TARG SNR_BPzerod.png'])
        %also print log/linear var
        %     %%
%             figure(2)
%             icounter=1;
%         
%             for ifreq=1:length(peakfreqsare)
%                subplot(6,2,icounter)
%                % pre orpost?
%                if mod(icounter,2)
%                    %even numbers
%                    windp=2;
%                else
%                    windp=1;
%                end
%                hist(mbar(windp, ifreq));
%                hold on
%                hist(log(mbar(windp,ifreq,:)));
%                icounter=icounter+1;
%             end
%         hzUsed=hz;
%
%take mean across locs
if any(isnan(ppantRESSsanitycheck(:)))
    error('')
end
ppantRESSsanitycheck=squeeze(nanmean(ppantRESSsanitycheck,2));
        %%
        save('ppantRESSsanitycheck', 'ppantRESSsanitycheck', 'f', 'windowsmall');
    end
end
%%
if job.concatSanitycheck==1
    cd(basefol);
    peakfreqsare=[8,13,15,18,20,40,16,26,30,36];
    storeall=[];
    snrall=[];
    icounter=1;
    %%
    for ifol=allppants
        
        cd(basefol)
        cd(num2str(ifol));
        load('ppantRESSsanitycheck.mat')
        for ifreq=1:length(peakfreqsare)
            % also collect SNR across trials for SD/variance estimate.
            peakfreq1=peakfreqsare(ifreq);
            [~,fis]=min(abs(f-peakfreq1));
            
            storeall(icounter,:,ifreq) = squeeze(ppantRESSsanitycheck(:,ifreq,fis));
            
            snrall(icounter,:,ifreq,:) =squeeze(ppantRESSsanitycheck(:,ifreq,:,:));
        end
        icounter=icounter+1;
    end
    %%
    hzUsed=f;
    cd(basefol)
    cd('newplots-MD')
    save('RESS-sanitycheck_Catchbased', 'storeall', 'snrall', 'hzUsed')
    %%
    
    %% export for stats analysis
    %
    TableforExport=table();
    thisentry=1; %initialize counter.
    %
    for iwin = 1:2
        for ihz=1:length(peakfreqsare)
            
            %loc data =
            ppantdata = squeeze(storeall(:, iwin, ihz));
            
            %append horizontally to table.
            TableforExport=[TableforExport, table(ppantdata)];
            TableforExport.Properties.VariableNames(thisentry)={['Tidx' num2str(iwin) '_' num2str(peakfreqsare(ihz)) 'Hz']};
            thisentry=thisentry+1;
        end
    end
    
    try cd('/Users/matthewdavidson/Desktop')
    catch
        cd('/Users/mattdavidson/Desktop')
    end
    
    writetable(TableforExport, 'RESSprevspost.csv')
end

if job.plotSanitycheckprepostRESSoncatch==1
    
    
    % load
    %Spectrum and barchart.
    peakfreqsare=[8,13,15,18,20,40,16,26,30,36];
    cd(basefol)
    cd('newplots-MD')
    load('RESS-sanitycheck_Catchbased.mat')
    hfig=figure(1);
    clf

    
%reorder HZ in bar plot.

 storeallRearranged = storeall(:,:,[1, 2, 7,3,4, 5,8,9,10,6]) ;
 
%else the order is first harms, second harms.
storeallRearranged = storeall;

 %rearrange into a single vector. all present then all absent
    tmpD = [[squeeze(storeallRearranged(:,1,:))],[squeeze(storeallRearranged(:,2,:))]];
    %or each hz?
%     tmpD = [[squeeze(storeallRearranged(:,1,:))],[squeeze(storeallRearranged(:,2,:))]];

subplot(2,2,3:4)


    mbar= squeeze(nanmean(storeallRearranged ,1))';
    tmp=[];
    mbar= mbar([1,2,3,4,5,7,8,9,10,6],:);
    %order for stacked bar plots, all first harmonics then all second.
    for ihz = 1:size(mbar,1)
    tmp= [tmp;mbar(ihz,1), mbar(ihz,2)]; 
    end
    %
    %
    bh=bar(tmp);
    bh(2).FaceColor=[.7 .7 .7];
    
    %now means are plotted,
    %%need to rearrange for calc of w/in subj error bars:

    X=storeallRearranged ;
    % participant mean...
    pX=squeeze(mean(nanmean(storeallRearranged ,2),3));

    meanOverall = mean(pX(:));
    %
%adjust outgoing according to Cousineau (2005).
    adjustedST = tmpD - repmat(pX,[1 20])+ repmat(meanOverall, [size(tmpD,1) size(tmpD,2)]);
    
    barstore=squeeze(nanstd(adjustedST))./sqrt(size(adjustedST,1));
    %
    %now rearrange to place on barchart (ie permute dimensions).
    %order for stacked bar plots at the moment barstore is all present then
    %all absent.
    errbars =[];
    errbars(:,1) = barstore(1:10);
    errbars(:,2) = barstore(11:20);
    
    %now match original plot. (errbar is not stacked)
    
    errbars= errbars([1,2,3,4,5,7,8,9,10,6],:);
    %
    %order for stacked bar plots, all first harmonics then all second.
    tmp=[];
    for ihz = 1:size(errbars,1)
    tmp= [tmp;errbars(ihz,1), errbars(ihz,2)]; 
    end
    
    
    numgroups = length(peakfreqsare);
    numbars = 2;%length(peakfreqsare);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    % determine locations for errorbars:
    hold on
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, mbar(:,i), errbars(:,i), 'k', 'linestyle', 'none', 'linew',1);
    end
    %
%     title({['Pre and post ' dtypep ' onset']})
  set(gca,'xticklabel',{'8', '13', '15' ,'18', '20', '16','26','30', '36', '40'}, 'fontsize', 20)   
%     set(gca,'xticklabel',{'8', '13', '15' ,'16', '18', '20', '26', '30', '36','40'}, 'fontsize', 20)
    
    xlabel('Hz')
    ylabel('log(SNR)')
%     title('RESS-SNR')
%      lga=legend('Targets physically present', 'Targets physically absent');
%      set(lga, 'location', 'northwest');
    set(gcf, 'color', 'w')
%     xlim([0 13])
    
     
    % can also plot color code:
    xt=get(gca, 'xtick');
    for ix= 1:length(xt)
        if ix<5
        lg2(1)=plot([xt(ix) xt(ix)], [-1 0], ['r-'], 'linew', 3);
        elseif ix==5
            lg2(2)=plot([xt(ix) xt(ix)], [-1 0], ['r:'], 'linew', 3);
        elseif ix<10
            lg2(3)=plot([xt(ix) xt(ix)], [-1 0], 'color', [0 .5 0], 'linew', 3);
        elseif ix==10
            lg2(4)=plot([xt(ix) xt(ix)], [-1 0], 'color', [0 .5 0], 'linestyle', ':','linew', 3);
        end
    end
    %
    %make new axis for the new legend.     
    
%     ax=axes('Position',get(gca,'Position'),'Visible','Off');
    % plot second legend. first copy the first so we retain it.
%     lg2a=legend(ax, [lg2(:)],{'Targ 1st harm.','Background first harm.','Targ 2nd harm.','Background 2nd harm.'});
%%
    lg2a=legend([lg2(:)],{'Targ 1st harm.','Background first harm.','Targ 2nd harm.','Background 2nd harm.'});
    set(lg2a, 'Location', 'NorthWest')
     %%
         
     % also alternate for paper fig
    
    
    plotmeanver=1; % or 2 to plot hz seprate.
    for isub=1:2
        subplot(2,2,isub)
        switch isub
            case 1
    usefreqs=[1:4,7:10]; %1:4 = 1fTG, 7:10 = 2fTGs
            case 2            
                usefreqs=[5,6]; %1:4 = 1fTG, 7:10 = 2fTGs
        end
%     clf
    
    if plotmeanver==1
        
        mbar=mean(squeeze(mean(storeall(:,:,usefreqs),3)),1);
        
        bh=bar(mbar);
        hold on
        %%
        bar(2,mbar(2),'FaceColor',[.7 .7 .7])
        %%
% now adjust w/subj error bars.
        errbarstore= mean(squeeze(storeall(:,:,usefreqs)),3);
        
        hold on
        % adjust errorbars.
        pM = mean(errbarstore,2);
        mM= mean(pM);
        adj = errbarstore - repmat(pM,1,2) + repmat(mM,size(storeall,1) ,2);
        
        errplot = std(adj)/sqrt(21);
        %%
        errorbar(1:2, mbar, errplot, 'k', 'linestyle', 'none', 'linew', 1)
%         lga=legend({'8 Hz', '13 Hz', '15 Hz' ,'18 Hz'});         
    else
    bar(mbar2)
    hold on
    barstore2= barstore(usefreqs,:)';    
     numgroups = 2;
    numbars = length(usefreqs);%length(peakfreqsare);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    % determine locations for errorbars:
    
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        et=errorbar(x, mbar2(:,i), barstore2(:,i), 'k', 'linestyle', 'none');
    end
    lga=legend({'8 Hz', '13 Hz', '15 Hz' ,'18 Hz','16','26','30','36'}); 
    set(lga, 'location', 'northEast');
    end
    
    set(gca,'xticklabel',{'Targets Present', 'Targets Absent'}, 'fontsize', 20)
    
    ylabel('log(SNR)')
     switch isub
         case 1
     ylim([0 .8]);
     title({['Target Flicker'];['[8,13,15,16,18,26,30,36 Hz]']})
         case 2
             ylim([ 2.8 3.6])
             title({['Background Flicker'];['[20,40 Hz]']})
     end
%      axis tight
%     xlim([.2 3])
    end
    %%
    cd(basefol)
    cd('newplots-MD')
    cd('figures')
    print('-dpng', 'RESS-sancheck1')
end
%%

if job.performLMEonRESS1 ==1
   cd('/Users/matthewdavidson/Desktop/stats for JASP') 
    datatable=readtable('RESSprevspost_longformTARGETS.csv');
    
    %perform LME.
    icount=1;
    for idrop = 1:4 % which terms to check for sig?
        switch idrop
            
            case 1 %test interaction.
                UNterm1='SNR ~ x___Flicker*Targets + (1|SubjectID)'; %note no harmonic. (nonsig)
                UNterm2='SNR ~ x___Flicker*Targets - x___Flicker:Targets + (1|SubjectID)';
                dropis= 'FE Int';
            case 2 %test Target flicker freq.
                UNterm1='SNR ~  x___Flicker*Targets +  (1|SubjectID)';
                UNterm2='SNR ~  x___Flicker*Targets - x___Flicker + (1|SubjectID)';
                dropis= 'FE Flicker';
            case 3 %test real vs shuff
                UNterm1='SNR ~  x___Flicker*Targets +  (1|SubjectID)';
                UNterm2='SNR ~  x___Flicker*Targets - Targets + (1|SubjectID)';
                dropis= 'FE Target presence';
        end
        
        datatable.SubjectID=nominal(datatable.SubjectID);
        datatable.x___Flicker=nominal(datatable.x___Flicker);
        datatable.Targets =nominal(datatable.Targets);
        
            
                %creating models:
                lmeUN =fitlme(datatable, UNterm1);%, 'fitMethod', 'REML');
                lmeRe =fitlme(datatable, UNterm2);%,'fitMethod', 'REML';
                
                tmp=compare(lmeRe, lmeUN);
                    [h,pValue,stat] =LogLikelihoodReport(lmeUN,lmeRe);
                
                %%
                %store chi-sq;
                tresults(icount,1) = table({'SNR'});
                tresults(icount,2) = table({dropis});
                tresults(icount,3) = table(tmp.LRStat(2));
                tresults(icount,4) = table(tmp.deltaDF(2));
                tresults(icount,5) = table(tmp.pValue(2));
                
                
                icount=icount+1;
    end
    
     tresults.Properties.VariableNames(1)= {'DV'};
        tresults.Properties.VariableNames(2)= {'FixedEffect'};
        tresults.Properties.VariableNames(3)= {'teststat'};
        tresults.Properties.VariableNames(4)= {'DoF'};
        tresults.Properties.VariableNames(5)= {'pval'};
        %display results
        disp(tresults)
end


if job.construct_and_compareRawtoRESSonBGhz==1
    %% first job is to save SNR per ppant, in Timewindow 1 for PFI increase 0_1( TG present)
    % and CAtCH 1 (TG PRESENT).
    
    %params for SNR>
    srate=250;
    epochdur = sum(abs([-3 3]))*srate;
    timeid = [0:1/srate:epochdur];
    timeid= timeid-3;
    
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.fpass= [0 50];
    param_spctrm.trialave=0;
    
    tcheck = [-3 -.1];
    timewindow = dsearchn([timeid]', [tcheck(1) tcheck(2)]');
    
    usechan=30;
    acrossALLRawvsRESS= zeros(length(allppants), 2, 2, 2);
%     [nppants, nhz, ndataPROCESSEDtype, ndataEXPtype]= size(acrossALLRawvsRESS);
    
    icounter=1;
    for ifol=allppants
        
        cd(basefol)
        cd(num2str(ifol))
        
        
        %load RESS TS
        load('ppant_Catch_Epoched_RESS')
        load('ppant_Catch_Epoched')
        load('ppant_PFI_Epoched')
        load('ppant_PFI_Epoched_RESS')
        % get SNR for 20 and 40 Hz.
        for id=1:4
            datause=[];
            switch id
                case 1 %raw Catch
                    datause = squeeze(ppant_SNREEG_catchramponset(:,usechan,timewindow(1):timewindow(2)));
                    nPROC=1; %raw
                    nEXP=1; %catch
                case 2 %rawPFI
                    datause = squeeze(ppant_SNREEG_PFI_0_1(:,usechan,timewindow(1):timewindow(2)));
                    nEXP=2;
%                     trialstouse = ress_PFI_0_1_HzxLoc.retainedtrialIndex;
                case 3 %ress catch. locs are repeated. 
                    datause(1,:,:) = squeeze(ress_catchonsetBGs(5,1,:,timewindow(1):timewindow(2)));
                    datause(2,:,:) = squeeze(ress_catchonsetBGs(6,1,:,timewindow(1):timewindow(2)));
                    nPROC=2;
                    nEXP=1;
                case 4 %ress PFI
                    datause(1,:,:) = squeeze(ress_PFI_0_1_HzxLoc(5,1).ressTS(:,timewindow(1):timewindow(2)));
                    datause(2,:,:) = squeeze(ress_PFI_0_1_HzxLoc(6,1).ressTS(:,timewindow(1):timewindow(2)));
                    nEXP=2; 
            end
            
            
            if sum(unique(isnan(datause(:))))>0
                %remove nans from dataset to avoid stuffing spectrum.
                if id<3 % trials is first dimension    
                alltrials = [1:size(datause,1)];
                    ntrials = find(isnan(datause(:,1)));
                    alltrials(ntrials)=[];
                    datause= datause(alltrials,:);
                else % trials are second dimension.
                    alltrials = [1:size(datause,2)];
                    ntrials = find(isnan(datause(1,:,1)));
                    alltrials(ntrials)=[];
                    datause= datause(:,alltrials,:);
                end
                    
                    
            elseif sum(unique(isinf(datause(:))))>0
            error('pause');
            end
            
            %take SNR mtspectrum of whole EPOCH, 2 xSNR to extract, since
            %not RESS data.
            if id<3
                
                [s,f] = mtspectrumc(datause', param_spctrm);
                
                s= log(s);
                mS = squeeze(nanmean(s,2));
                %take SNR,
%                 kernelw = [-1/6 -1/6 -1/6   0 0 1  0 0 -1/6  -1/6 -1/6];
                kernelw = [-1/4 -1/4   0 0 1  0 0   -1/4 -1/4];
                nSNR = conv(mS, kernelw, 'same');
                
                %store each freq
                [~, hz20] = min(abs(f-20));
                [~, hz40] = min(abs(f-40));
                
                
                acrossALLRawvsRESS(icounter,1,nPROC,nEXP)=nSNR(hz20);
                acrossALLRawvsRESS(icounter,2,nPROC,nEXP)=nSNR(hz40);
                
            else %take separate SNRs per RESS component TS.
                
                for ifreq=1:2
                    dtmp= squeeze(datause(ifreq,:,:));
                    
                    [s,f] = mtspectrumc(dtmp', param_spctrm);
                    
%                     s= 10*log10(s.^2);
                    s= log(s);
                    mS = squeeze(nanmean(s,2));
                    %take SNR,
%                     kernelw = [-1/6 -1/6 -1/6   0 0 1  0 0 -1/6  -1/6 -1/6];
                    kernelw = [-1/4 -1/4   0 0 1  0 0   -1/4 -1/4];
                    nSNR = conv(mS, kernelw, 'same');
                    
                    %store each freq
                    switch ifreq
                        case 1
                            [~, hzid] = min(abs(f-20));
                        case 2
                            [~, hzid] = min(abs(f-40));
                    end
                    
                    acrossALLRawvsRESS(icounter,ifreq,nPROC,nEXP) = nSNR(hzid);
                    
                    
                end
                
                
            end
            
        end
    icounter=icounter+1;
    end
    
    % export to desktop for JASP analysis.
    
    %%
    TableforExport=table();
    thisentry=1; %initialize counter.
    %
    peakfreqsare=[20,40];
    for iPROC = 1:2
        for ihz=1:2
            for iEXP=1:2
            %loc data =
            ppantdata = squeeze(acrossALLRawvsRESS(:, ihz, iPROC, iEXP));
            
            %append horizontally to table.
            TableforExport=[TableforExport, table(ppantdata)];
            TableforExport.Properties.VariableNames(thisentry)={['dataRawRess' num2str(iPROC) '_' num2str(peakfreqsare(ihz)) 'Hz_CATCHPFI' num2str(iEXP)]};
            thisentry=thisentry+1;
            end
        end
    end
    
    try cd('/Users/matthewdavidson/Desktop')
    catch
        cd('/Users/mattdavidson/Desktop')
    end
    
    writetable(TableforExport, 'RESSSanityCheck_BGs.csv')
    %%
    
    
    % take mean across both
    %% ppant level subtraction
    [nppants, nhz, ndataPROCESSEDtype, ndataEXPtype]=size(acrossALLRawvsRESS);
    clf

    %organise across all for plot
% for easier calcs:
    plDiff = [squeeze(acrossALLRawvsRESS(:,1,1,1)), squeeze(acrossALLRawvsRESS(:,1,1,2)), squeeze(acrossALLRawvsRESS(:,2,1,1)), squeeze(acrossALLRawvsRESS(:,2,1,2)),...
        squeeze(acrossALLRawvsRESS(:,1,2,1)),squeeze(acrossALLRawvsRESS(:,1,2,2)),  squeeze(acrossALLRawvsRESS(:,2,2,1)),squeeze(acrossALLRawvsRESS(:,2,2,2))];
    
    stE = std(plDiff)/sqrt(size(plDiff,1));
    
    
    mAll = squeeze(nanmean(plDiff,1));
    %all raw vs RESS
    
    plotmean= [mAll(1:4);mAll(5:8)];
    plotErr = [stE(1:4);stE(5:8)];
    
    subplot(211)
    b1=bar(plotmean);
    b1(1).FaceColor= 'r';
    b1(2).FaceColor= 'r';
    b1(2).EdgeColor= [.2 .2 .2];
    b1(2).LineWidth= 5;
    b1(3).FaceColor= [0 .5 0];
    b1(4).FaceColor= [0 .5 0];
    b1(4).EdgeColor= [.2 .2 .2];
    b1(4).LineWidth= 5;
    
    %
    hold on
     numgroups = 2;
    numbars = 4;%length(peakfreqsare);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    % determine locations for errorbars:
    
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, plotmean(:,i), plotErr(:,i), 'k', 'linestyle', 'none');
    end
%   title({['Comparing SNR changes introduced by RESS'];['timewindow = ' num2str(tcheck(1)) ':' num2str(tcheck(2)) 's']});
%   title({['Comparing SNR changes introduced by RESS'];['timewindow = ' num2str(tcheck(1)) ':' num2str(tcheck(2)) 's']});
    lg=legend({'20 Hz pre PFI','20 Hz pre Catch',  '40 Hz pre PFI','40 Hz  pre Catch'});
    set(lg, 'location', 'northeast')
    ylabel('log(SNR)')
    set(gca, 'xtick', 1:2, 'xticklabel', {'before RESS', 'after RESS'}, 'fontsize', 20)
    xlim([.5 3.5])
    ylim([0 4])
%
    subplot(212); %plot diff: PFI - catch 
    title({['Difference in SNR '];['[pre PFI - pre Catch]']})
    hold on
    ressdiffAll = squeeze(acrossALLRawvsRESS(:,:,2,1))-squeeze(acrossALLRawvsRESS(:,:,2,2));
    rawdiffAll= squeeze(acrossALLRawvsRESS(:,:,1,1))-squeeze(acrossALLRawvsRESS(:,:,1,2));
    %calculate stErr
    plDiff = [rawdiffAll, ressdiffAll];
    
    stE = std(plDiff)/sqrt(size(plDiff,1));
    
    %rearrange for grouped barchart:
    mDiff = [squeeze(nanmean(rawdiffAll(:,1),1)),squeeze(nanmean(rawdiffAll(:,2),1)); ...
        squeeze(nanmean(ressdiffAll(:,1),1)), squeeze(nanmean(ressdiffAll(:,2),1))];
    
    stEe= [stE(1), stE(2);stE(3), stE(4)];
    
   b2=bar(mDiff);
   b2(1).FaceColor = 'r';
   b2(2).FaceColor = [0 .5 0];
   %plot errorbars
   hold on;
  numgroups = 2;
    numbars = 2;%length(peakfreqsare);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    % determine locations for errorbars:
    
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, mDiff(:,i), stEe(:,i), 'k', 'linestyle', 'none');
    end
%     eb=errorbar(1:4, squeeze(mean(plDiff,1)), stE, 'LineStyle', 'none');
xlim([.5 3.5])
ylim([ -.15 .3])
legend({' 20 Hz ', ' 40 Hz'})   

set(gca, 'xtick', 1:2,'xticklabel', {'before RESS', 'after RESS'}, 'fontsize', 20)
   ylabel('\Delta log(SNR)')
   set(gcf, 'color', 'w')
   shg
    %%
     cd(basefol)
    cd('newplots-MD')
    cd('figures')
    print('-dpng', 'RESS-sancheck2')
end


