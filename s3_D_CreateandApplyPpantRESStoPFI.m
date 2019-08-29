
clear all
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/MatthewDavidson/Desktop/FromIrenes Mac/Data Experiment copy')
end

basefol=pwd;
clearvars -except basefol allppants
dbstop if error

%%%%%% NOTE that at this stage, we sort to include ONLY relevant
%%%%%% disappearances and reappearances, where the RESS filter (HzxLoc) matches
%%%%%% the disap type (HzxLoc).

%
job.applyRESSperppant_toPFI=1;
%
% job.constructSNRpreandpostPFIperppant=0;
%
%
% job.crunchtostaticSNRperppant=0;
% job.concatstaticSNRacrossppants=0;
%
% job.plotstaticSSVEPSNR=1;
%
%
% %% also dynamic SSVEP
%
% job.constructdynSSVEPperppant=0;
% job.concatacrossppants=1;
%
% job.plotPFIdynSSVEP=1;



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
baserem1 = [-2.99 -2]; % baseline to remove in secs for ONSET of catch
baserem2 = baserem1;%[2 2.99]; % baseline to remove in secs for offset of catch
epochdur = sum(abs(window))*250;
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/250:60;


% which frequencies to analyze?/apply?
        peakfreqsare=[8,13,15,18,20,40];
%         peakfreqsare=[16,26,30,36]; % can change to 2*f of TGs
        

if job.applyRESSperppant_toPFI==1
    for ifol =allppants
        
        
        %%
        cd(basefol)
        cd(num2str(ifol))
        
        
        %% load the relevant PFI data.
        load('ppant_PFI_Epoched');
        if length(peakfreqsare)==4
        load('RESSfilterdata_2fTG')
        else
            load('RESSfilterdata')
        end
        load('TrialIndicesbyLocationandHz.mat')

        dataIN=[];
        
         %%
                for id=1:6% Use all epochs so as not to bias condition comparisons.
                    
                    trialsbyHz=[];
                    switch id
                        case 1
                            dataIN=ppant_SNREEG_PFI_0_1;
                            searchtrials = [Freqwas.dir0_1.Trialind];
                            trialsbyHz = Freqwas.dir0_1;
%                             durscheck=durs0_1;
                            ress_PFI_0_1_HzxLoc=[]; %we don't know how many of each type!
                            BPstokeep= BPs0_1;
                        case 2
                            dataIN=ppant_SNREEG_PFI_1_0;
                            searchtrials = [Freqwas.dir1_0.Trialind];
                            trialsbyHz = Freqwas.dir1_0;
%                             durscheck=durs1_0;
                            ress_PFI_1_0_HzxLoc=[]; 
                            
                            BPstokeep= BPs1_0;
                            
                        case 3
                            dataIN=ppant_SNREEG_PFI_1_2;
                            searchtrials = Freqwas.dir1_2.Trialind;
                            % need to make it easy for the trials to be
                            % recognised: 
                            
                            trialsbyHz = Freqwas.dir1_2;
%                             durscheck=durs1_2;
                            ress_PFI_1_2_HzxLoc=[]; %we don't know how many of each type!
                            BPstokeep= BPs1_2;
                        case 4
                            dataIN=ppant_SNREEG_PFI_2_1;
                            searchtrials = Freqwas.dir2_1.Trialind;
                            trialsbyHz = Freqwas.dir2_1;
%                             durscheck=durs2_1;
                            ress_PFI_2_1_HzxLoc=[]; %we don't know how many of each type!
                            BPstokeep= BPs2_1;
                        case 5
                            dataIN=ppant_SNREEG_PFI_2_3;
                            searchtrials = [Freqwas.dir2_3.Trialind];
                            trialsbyHz = Freqwas.dir2_3;
%                             durscheck=durs2_3;
                            ress_PFI_2_3_HzxLoc=[]; %we don't know how many of each type!
                            BPstokeep= BPs2_3;
                            
                        case 6
                           
                            dataIN=ppant_SNREEG_PFI_3_2;                            
                            searchtrials = [Freqwas.dir3_2.Trialind];
                            trialsbyHz = Freqwas.dir3_2;
%                             durscheck=durs3_2;
                            ress_PFI_3_2_HzxLoc=[]; %we don't know how many of each type!
                            BPstokeep= BPs3_2;
                           
                       
                            
                    end
                    %%
            %compare eeg 
%             clf
%                     plot(squeeze(mean(ppant_SNREEG_PFI_0_1(:,32,:),1)));
%                     hold on;
%                     plot(squeeze(mean(ppant_SNREEG_PFI_1_2(:,32,:),1)));
%                     plot(squeeze(mean(ppant_SNREEG_PFI_2_3(:,32,:),1)));
% %                     
%                     plot(squeeze(mean(ppant_SNREEG_PFI_1_0(:,32,:),1)));
%                     plot(squeeze(mean(ppant_SNREEG_PFI_2_1(:,32,:),1)));
%                     plot(squeeze(mean(ppant_SNREEG_PFI_3_2(:,32,:),1)));
%                     shg
                    %%
                    
            for ifreq=1:length(peakfreqsare)
                usehz=peakfreqsare(ifreq);
                
                for iloc=1:4
                    % apply the spatial filter per trial
                    %%
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
                    
                    
%which trials? (1:24) are the appropriate hz x location for RESS.

                    trialtypesperTargPresent = find(ismember(searchtrials, relevanttrials));
%                    
% if trialtypesperTargPresent >0 % if we have trials available.
if ifreq<=4
    %critically for TG Hz, we only want to retain the
    %disappearances (ie. BP that involved a particular HZ)
    %                       keeptrials=[]
    
    %first we need to restrict the DISAP to only those in the correct trial
    %types outlined above.
    
    %all Disaps:
    
    
    keeptrials=[];
    for item = 1:size(trialsbyHz,2) %this is needed, as some disaps had simultaneous BP, messes with indexing.
        
        try hzc=trialsbyHz(item).d; % are there trials of this type?  (some ppants did not have enough disapperances of >2 targets for example).
        catch
            continue
        end
        if peakfreqsare(1)==8; %first harmonic 
            if any(hzc==usehz)
                keeptrials=[keeptrials, (item)];
                
            end
        else
            if any(hzc==usehz/2)
                keeptrials=[keeptrials, (item)];
                
            end
        end
            
        
        
    end
    
    %now restrict the 'saved' trials to only those listed.
    trialtypesperTargPresent = find(ismember(trialtypesperTargPresent, keeptrials));
        
end


                    if ~isempty(trialtypesperTargPresent)
                        %reduce size.
                        if length(trialtypesperTargPresent)<2
                            datast(1,:,:)= squeeze(dataIN(trialtypesperTargPresent,:,:));
                        else
                            datast= squeeze(dataIN(trialtypesperTargPresent,:,:));
                        end
                        
                        
                        
                        % remove bad trials.
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
%                     % also skip the trials which were only transient button
%                     % presses. (less than one second).                    
%                     shorttrials = find(durscheck<60);
%                     badtrials = [badtrials, shorttrials];
%                     
%                     % remove these from consideration.
%                     keeptrials(badtrials)=[];
%                     datast=datast(keeptrials,:,:);
%                      

 

                        
                    BPstosave = BPstokeep(trialtypesperTargPresent,:);
                        
                        %now we have the correct trials, get the appropriate
                        %filter
                        evecs = squeeze(ressEVEC_byHzxLoc(ifreq,iloc,:)); %
                        
                        
                        
                        ress_ts1=zeros(size(datast,1), size(datast,3));
                        for ti=1:size(datast,1)
                            ress_ts1(ti,:) = evecs'*squeeze(datast(ti,:,:));
                        end
                        
                    else
                        ress_ts1=nan(1,1501);
                        BPstosave=nan(1,1501);
                    end
                    
                    %Save the ress data per type.
                    %also keep BP data which we will need for plotting.
                    switch id
                        case 1
                            ress_PFI_0_1_HzxLoc(ifreq,iloc).ressTS= ress_ts1;
                            ress_PFI_0_1_HzxLoc(ifreq,iloc).BPs = BPstosave;
                            ress_PFI_0_1_HzxLoc(ifreq,iloc).retainedtrialIndex= trialtypesperTargPresent;
                            
                            
                        case 2
                            
                            ress_PFI_1_0_HzxLoc(ifreq,iloc).ressTS= ress_ts1;
                            ress_PFI_1_0_HzxLoc(ifreq,iloc).BPs = BPstosave;
                            ress_PFI_1_0_HzxLoc(ifreq,iloc).retainedtrialIndex= trialtypesperTargPresent;
                            
                            
                        case 3
                            
                            ress_PFI_1_2_HzxLoc(ifreq,iloc).ressTS= ress_ts1;
                            ress_PFI_1_2_HzxLoc(ifreq,iloc).BPs = BPstosave;
                            ress_PFI_1_2_HzxLoc(ifreq,iloc).retainedtrialIndex= trialtypesperTargPresent;
                            
                            
                            
                        case 4
                            
                            ress_PFI_2_1_HzxLoc(ifreq,iloc).ressTS= ress_ts1;
                            ress_PFI_2_1_HzxLoc(ifreq,iloc).BPs = BPstosave;
                            ress_PFI_2_1_HzxLoc(ifreq,iloc).retainedtrialIndex= trialtypesperTargPresent;
                            
                        case 5
                            
                            ress_PFI_2_3_HzxLoc(ifreq,iloc).ressTS= ress_ts1;
                            ress_PFI_2_3_HzxLoc(ifreq,iloc).BPs = BPstosave;
                            ress_PFI_2_3_HzxLoc(ifreq,iloc).retainedtrialIndex= trialtypesperTargPresent;
                            
                        case 6
                            
                            ress_PFI_3_2_HzxLoc(ifreq,iloc).ressTS= ress_ts1;
                            ress_PFI_3_2_HzxLoc(ifreq,iloc).BPs = BPstosave;
                            ress_PFI_3_2_HzxLoc(ifreq,iloc).retainedtrialIndex= trialtypesperTargPresent;
                            
                            
                    end
                end
            end
            clearvars dataIN
            
                end
        if length(peakfreqsare)==4
            savename='ppant_PFI_Epoched_RESS_2fTG';
        else
            savename='ppant_PFI_Epoched_RESS';
        end
            
        save(savename, 'ress_PFI_0_1_HzxLoc', 'ress_PFI_1_0_HzxLoc',...
            'ress_PFI_1_2_HzxLoc','ress_PFI_2_1_HzxLoc', 'ress_PFI_2_3_HzxLoc', 'ress_PFI_3_2_HzxLoc')
        display([' finished ppant ' num2str(ifol)])
    end
end

% 
% 
% if job.constructSNRpreandpostPFIperppant==1
%     
%     % Compute the average SNR for PFI pre PFI BP
%     % Store and compare to post.
%     
%     
%     srate=250;
%     % this bit hard coded, 6 s window length
%     epochdur = sum(abs([-3 3]));
%     timeid = 0:1/srate:epochdur;
%     timeid= timeid-3;
%     
%     windowsmall=[];
%     %     windowsmall(1,:) = [-3 -0.05] ;% short window targets present
%     %     windowsmall(2,:) = [.05 3] ;% short window targets present
%     windowsmall(1,:) = [-3 3] ;% short window targets present
%     windowsmall(2,:) = [-3 3] ;% short window targets present
%     
%     nfft = ceil( srate/.1 ); % .1 Hz resolution
%     hz    = linspace(0,srate,nfft);
%     
%     markert='o';
%     
%     for ifol=allppants
%         cd(basefol)
%         cd(num2str(ifol))
%         %%
%         
%         load('ppant_PFI_RESS')
%         
%         
%         %%
%         
%         
%         
%         for dtype=1:2%
%             switch dtype %obtain correct RESS timeseries.
%                 case 1
%                     allTS=ressPFI_0_1.ressTS;
%                     freqstocheck = Freqwas.dir0_1;
%                     ddir='0 -> 1';
%                     dtype1= 'visible';
%                     dtype2= 'PFI';
%                     
%                     cols(1).c = 'k';
%                     cols(2).c = 'b';
%                 case 2
%                     allTS=ressPFI_1_0.ressTS;
%                     freqstocheck = Freqwas.dir1_0;
%                     ddir='1 -> 0';
%                     dtype1= 'PFI';
%                     dtype2= 'visible';
%                     
%                     cols(1).c = 'b';
%                     cols(2).c = 'k';
%                     %                 case 3
%                     %                     dataIN=ppant_SNREEG_PFI_1_2;
%                     %                 case 4
%                     %                     dataIN=ppant_SNREEG_PFI_2_1;
%                     %                 case 5
%                     %                     dataIN=ppant_SNREEG_PFI_2_3;
%                     %                 case 6
%                     %                     dataIN=ppant_SNREEG_PFI_3_2;
%             end
%             
%             
%             peakfreqsare = [8,13,15,18, 20,40];
%             
%             
%             ppantRESSsanitycheck =[];%zeros(2,length(peakfreqsare), 2500); %hard coded. change if window changes.
%             mbar=[];barstore=[];
%             
%             
%             
%             
%             for ifreq=1:length(peakfreqsare)
%                 peakfreq1= peakfreqsare(ifreq);
%                 %%
%                 
%                 
%                 
%                 % for each frequency, we want ONLY the trials that involved PFI
%                 % at that frequency.
%                 
%                 % all hz involved:
%                 hzinvolved= [freqstocheck(:).d];
%                 
%                 crunchtrials=find(hzinvolved==peakfreq1);
%                 
%                 if ifreq>4
%                     %use all trials for 20Hz construction?
%                     
%                     crunchtrials = 1:size(freqstocheck,2);
%                     
%                 end
%                 
%                 
%                 
%                 for wind=1:2
%                     
%                     
%                     %correct color.
%                     col=cols(wind).c ;
%                     
%                     
%                     
%                     windcheck= windowsmall(wind,:);
%                     tidx=dsearchn(timeid', windcheck');
%                     
%                     
%                     
%                     
%                     % take fft and SNR of pre disappearance, pre performed RESS
%                     % over whole EPOCH.
%                     crunchdata = squeeze(allTS(crunchtrials,ifreq, tidx(1):tidx(2)))';
%                     
%                     
%                     
%                     %%
%                     ressxALL =(abs( fft(crunchdata,nfft,1)/diff(tidx) ).^2);
%                     ressx=nanmean(ressxALL,2);
%                     
%                     % compute SNR>
%                     
%                     
%                     % compute SNR for this trial:
%                     [snrR,snrE] = deal(zeros(size(hz)));
%                     
%                     
%                     
%                     %%% very important, compare again to the correct frequencies
%                     %we have used
%                     %low shoulder
%                     
%                     skipbinsLo =  dsearchn(hz', [peakfreq1-neighbour(ifreq).snrlow(2)]'); % hard-coded in Hz, how many away before start shoulder.
%                     numbinsLo  = dsearchn(hz', [peakfreq1-neighbour(ifreq).snrlow(1)]'); % hard-coded in Hz, how many away before start shoulder.
%                     
%                     
%                     skipbinsHi =  dsearchn(hz', [neighbour(ifreq).snrhigh(1)- peakfreq1]'); % hard-coded in Hz, how many away before start shoulder.
%                     numbinsHi  = dsearchn(hz', [neighbour(ifreq).snrhigh(2) - peakfreq1]'); % hard-coded in Hz, how many away before start shoulder.
%                     
%                     % loop over frequencies and compute SNR
%                     for hzi=numbinsLo+1:length(hz)-numbinsHi-1
%                         
%                         numer = ressx(hzi);
%                         
%                         denom = nanmean( ressx([hzi-numbinsLo:hzi-skipbinsLo hzi+skipbinsHi:hzi+numbinsHi]) );
%                         
%                         snrR(hzi) = numer./denom;
%                         
%                     end
%                     
%                     
%                     
%                     
%                     
%                     
%                     %%% plot new SNR and save.
%                     figure(1);
%                     
%                     if wind==1 && ifreq==1
%                         %                     clf
%                     end
%                     plot(hz, snrR, [col markert '-' ])
%                     hold on
%                     xlim([0 45])
%                     % also collect SNR across trials for SD/variance estimate.
%                     %             [~,fis]=min(abs(hz-peakfreq1));
%                     %
%                     %             mbar(wind,ifreq,:) = nanmean(snrRall(fis,:));
%                     %             barstore(wind,ifreq,:) = nanstd(snrRall(fis,:));
%                     
%                     % now store
%                     ppantRESSsanitycheck(wind, ifreq,:)= snrR;
%                     
%                     
%                     
%                     % also collect SNR across trials for SD/variance estimate.
%                     [~,fis]=min(abs(hz-peakfreq1));
%                     
%                     mbar(wind,ifreq,:) = nanmean(snrR(fis));
%                     %                 barstore(wind,ifreq,:) = nanstd(snrRall(fis,:));
%                     
%                 end
%                 %          shg
%             end
%             set(gcf,'color', 'w')
%             legend({['Target ' dtype1], ['Target ' dtype2]}, 'location', 'northwest')
%             xlim([ 0 45])
%             %     ylim([0 50])
%             set(gca, 'fontsize', 20)
%             xlabel('Hz')
%             ylabel('SNR')
%             title({['Instances of PFI '];[ddir ',ntrials= ' num2str(length(crunchtrials))]})
%             %     print snr
%             print('-dpng', ['figure_PFI_SNR_BPzerod_' ddir '.png'])
%             %also print log/linear var
%             %     %%
%             %     figure(2)
%             %     icounter=1;
%             %
%             %     for ifreq=1:length(peakfreqsare)
%             %        subplot(6,2,icounter)
%             %        % pre orpost?
%             %        if mod(icounter,2)
%             %            %even numbers
%             %            windp=2;
%             %        else
%             %            windp=1;
%             %        end
%             %        hist(mbar(windp, ifreq));
%             %        hold on
%             %        hist(log(mbar(windp,ifreq,:)));
%             %        icounter=icounter+1;
%             %     end
%             hzUsed=hz;
%             switch dtype
%                 case 1
%                     snrSSVEP_PFI_0_1 = ppantRESSsanitycheck;
%                 case 2
%                     snrSSVEP_PFI_1_0 = ppantRESSsanitycheck;
%             end
%             %%
%             
%         end
%         save('ppant_PFI_SNR', 'snrSSVEP_PFI_0_1',...
%             'snrSSVEP_PFI_1_0','hzUsed', 'windowsmall');
%     end
% end
% 
% 
% if job.crunchtostaticSNRperppant==1
%     
%     for ifol=allppants
%         cd(basefol)
%         cd(num2str(ifol))
%         %%
%         tmp=[];
%         load('ppant_PFI_SNR')
%         
%         
%         tmp(1,:,:)= squeeze(snrSSVEP_PFI_0_1(1,:,:)); %targets present
%         tmp(2,:,:)= squeeze(snrSSVEP_PFI_1_0(2,:,:)); %targets present!
%         
%         ppant_PFI_SNR_0= squeeze(nanmean(tmp,1));
%         
%         
%         tmp(1,:,:)= squeeze(snrSSVEP_PFI_0_1(2,:,:)); %targets absent
%         tmp(2,:,:)= squeeze(snrSSVEP_PFI_1_0(1,:,:)); %targets absent!
%         
%         
%         ppant_PFI_SNR_1= squeeze(nanmean(tmp,1));
%         
%         save('ppant_PFI_SNR.mat', 'ppant_PFI_SNR_1', 'ppant_PFI_SNR_0', '-append');
%     end
% end
% 
% 
% if job.concatstaticSNRacrossppants==1
%     
%     across_PFI_SNR_1=zeros(length(allppants),6,size(snrSSVEP_PFI_0_1,3));
%     across_PFI_SNR_0=across_PFI_SNR_1;
%     counter=1;
%     for ifol=allppants
%         cd(basefol)
%         cd(num2str(ifol))
%         %%
%         
%         load('ppant_PFI_SNR')
%         
%         across_PFI_SNR_0(counter,:,:)=ppant_PFI_SNR_0;
%         across_PFI_SNR_1(counter,:,:)=ppant_PFI_SNR_1;
%         counter=counter+1;
%         
%     end
%     
%     cd(basefol)
%     cd('newplots-MD')
%     save('EEG_GrFX_PFI_SNR.mat', 'across_PFI_SNR_1', 'across_PFI_SNR_0', 'hzUsed')
%     
%     
% end
% 
% 
% if job.plotstaticSSVEPSNR==1
%     peakfreqsare= [8,13,15,18,20,40];
%     figure(1);
%     clf
%     load('EEG_GrFX_PFI_SNR.mat')
%     hz=hzUsed;
%     mbar=zeros(2,21,6);
%     p=[];
%     for id=1:2
%         switch id
%             case 1
%                 dtype= across_PFI_SNR_1;
%                 col='b';
%                 markert='o';
%             case 2
%                 dtype=across_PFI_SNR_0;
%                 
%                 col='k';
%                 markert='x';
%         end
%         snrR = squeeze(nanmean(dtype,1));
%         
%         
%         
%         %store for barplots of variance.
%         for ifreq=1:length(peakfreqsare)
%             
%             subplot(211)
%             
%             ptmp=  plot(hz, snrR(ifreq,:), [col markert '-' ]);
%             p(id)=ptmp;
%             
%             hold on
%             xlim([0 45])
%             xlabel('Hz')
%             ylabel('SNR')
%             
%             peakfreq1=peakfreqsare(ifreq);
%             % also collect SNR across trials for SD/variance estimate.
%             [~,fis]=min(abs(hz-peakfreq1));
%             
%             mbar(id,:,ifreq) = squeeze(dtype(:,ifreq,fis));
%             
%             
%         end
%     end
%     %
%     legend([p(1) p(2)], {'PFI', 'noPFI'})
%     set(gca, 'fontsize', 15);
%     %%
%     subplot(212)
%     
%     mplot=squeeze(nanmean(mbar,2))';
%     bar(mplot)
%     %
%     stE= squeeze(nanstd(mbar,0,2))/sqrt(21);
%     stE=stE';
%     hold on
%     
%     numgroups = 6;
%     numbars = 2;
%     groupwidth = min(0.8, numbars/(numbars+1.5));
%     for i = 1:numbars
%         % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
%         x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
%         errorbar(x, mplot(:,i), stE(:,i), 'k', 'linestyle', 'none');
%     end
%     
%     set(gca,'xticklabel', {'8' '13' '15' '18' '20' '40'})
%     xlabel('Hz')
%     ylabel('SNR')
%     lga=legend('PFI', 'noPFI');
%     set(lga, 'location', 'northwest');
%     set(gca, 'fontsize', 15);
%     set(gcf, 'color', 'w')
% end
% 
% 
% if job.constructdynSSVEPperppant==1
%     
%     rmvbaseEPOCH=2; %0 for no removal,  1 for trial baseline, 2 for pooled freq baseline (compare across conds).
%     baserem = [-2.99 -.1]; % baseline to remove in secs for ONSET of catch
%     
%     epochdur = sum(abs(window))*250;
%     onsetc = ceil(epochdur)/2;
%     srate=250;
%     % this bit hard coded, 6 s window length
%     epochdur = sum(abs([-3 3]));
%     timeid = [0:1/srate:epochdur];
%     timeid= timeid-3;
%     
%     
%     
%     
%     for ifol=allppants
%         cd(basefol)
%         cd(num2str(ifol))
%         %%
%         
%         load('ppant_PFI_RESS')
%         load('TrialIndicesbyLocationandHz.mat')
%         %%
%         
%         
%         ppantPFI_dynSSVEP=zeros(6,length(timeid));
%         
%         for dtype=1:2%:4;%
%             clf
%             switch dtype %obtain correct RESS timeseries.
%                 case 1
%                     allTS=ressPFI_0_1.ressTS;
%                     freqstocheck = Freqwas.dir0_1;
%                     ddir='0 -> 1';
%                     dtype1= 'visible';
%                     dtype2= 'PFI';
%                     
%                     cols(1).c = 'k';
%                     cols(2).c = 'b';
%                 case 2
%                     allTS=ressPFI_1_0.ressTS;
%                     freqstocheck = Freqwas.dir1_0;
%                     ddir='1 -> 0';
%                     dtype1= 'PFI';
%                     dtype2= 'visible';
%                     
%                     cols(1).c = 'b';
%                     cols(2).c = 'k';
%                     
%                     
%             end
%             
%             
%             peakfreqsare = [8,13,15,18, 20,40];
%             colsare(1).c=[1 .1 .1];
%             colsare(2).c=[1 .9 0];
%             colsare(3).c=[0 .6 0];%'g';
%             colsare(4).c=[0 0 .6];%'b';
%             colsare(5).c=[0 0  0];%k';
%             colsare(6).c=[.5 .5 .5]; %grey
%             
%             
%             outgoingDATA=zeros(length(peakfreqsare), 1501); %hard coded. change if window changes.
%             
%             mbar=[];barstore=[];
%             lg=[];
%             
%             for ifreq=1:length(peakfreqsare)
%                 peakfreq1= peakfreqsare(ifreq);
%                 %%
%                 
%                 
%                 
%                 
%                 
%                 
%                 % for each frequency, we want ONLY the trials that involved PFI
%                 % at that frequency.
%                 
%                 % all hz involved:
%                 hzinvolved= [freqstocheck(:).d];
%                 
%                 crunchtrials=find(hzinvolved==peakfreq1);
%                 
%                 if ifreq>4
%                     %use all trials for 20Hz construction?
%                     
%                     crunchtrials = 1:size(freqstocheck,2);
%                     
%                 end
%                 
%                 %%
%                 
%                 
%                 if length(crunchtrials)>1
%                     peakwidt=3;
%                     % take dynSSEP
%                     % over whole EPOCH.
%                     crunchdata = squeeze(allTS(crunchtrials,ifreq, :))';
%                     
%                     %filter first.
%                     
%                     
%                     %filter, but inlcude side lobes for temporal dynamics
%                     if size(crunchdata,1)==1
%                         crunchdata=crunchdata';
%                     end
%                     fdatAt = filterFGx(crunchdata',srate,peakfreq1,peakwidt,0);
%                     
%                     
%                     % take absolute of hilbert transform, to show dynamic SSEP
%                     dynSSEP= abs(hilbert(fdatAt'))';
%                     
%                     
%                     
%                     %%%%% perform any epoch preprocessing?
%                     % remove baseline
%                     if rmvbaseEPOCH ==1
%                         
%                         %find relevant window iF downsampled
%                         rme= onsetc + baserem*srate;
%                         
%                         rme=ceil(rme);
%                         %removes trial by trial baseline:
%                         
%                         rmdynSSEP=zeros(size(dynSSEP));
%                         for itrial=1:size(dynSSEP,1)
%                             %base line norm
%                             dynSSEPt= squeeze(dynSSEP(itrial,:));
%                             rmv = nanmean(dynSSEPt(1,rme(1):rme(2)));
%                             %subtract
%                             rmdynSSEP(itrial,:) = dynSSEPt - repmat(rmv, [1 size(dynSSEP,2)]);
%                             
%                         end
%                         dynSSEP=rmdynSSEP;
%                     elseif rmvbaseEPOCH==2
%                         
%                         % or remove baseline for all this freq..
%                         
%                         %find relevant window iF downsampled
%                         rme= onsetc + baserem*srate;
%                         rme=ceil(rme);
%                         rmdynSSEP=zeros(size(dynSSEP));
%                         
%                         rmv = mean(nanmean(dynSSEP(:,rme(1):rme(2)),2));
%                         for itrial=1:size(dynSSEP,1)
%                             %base line norm
%                             dynSSEPt= squeeze(dynSSEP(itrial,:));
%                             
%                             %subtract
%                             rmdynSSEP(itrial,:) = dynSSEPt - repmat(rmv, [1 size(dynSSEP,2)]);
%                             
%                         end
%                         dynSSEP=rmdynSSEP;
%                         
%                     end
%                     
%                     
%                     
%                     %Store mean per freq.
%                     outgoingDATA(ifreq,  :) = squeeze(nanmean(dynSSEP,1));
%                 else
%                     outgoingDATA(ifreq,  :) = nan(1,size(dynSSEP,2));
%                 end
%                 
%                 %%% plot new SNR and save.
%                 figure(1);
%                 
%                 mplot = squeeze(nanmean(dynSSEP,1));
%                 %
%                 %
%                 %             adjust for withinsub errorbars.
%                 
%                 step1 = squeeze(nanmean(dynSSEP,2)); %within ppant mean (across timepoints);
%                 step2 = mean(step1); %total exp mean.
%                 
%                 %subtract ppant mean from ppant., and group mean from all
%                 newd = dynSSEP- repmat(step1, 1, size(dynSSEP,2)) - repmat(step2, size(dynSSEP));
%                 
%                 %             newd=tmpd;
%                 %calc new stErr
%                 stErr = nanstd(newd)/sqrt(size(newd,1));
%                 
%                 %plot
%                 lgt=shadedErrorBar(timeid, mplot, stErr, ['k'],[1]);
%                 lgt.mainLine.Color = colsare(ifreq).c;
%                 lgt.patch.FaceColor = colsare(ifreq).c;
%                 lgt.edge(1).Color = colsare(ifreq).c;
%                 lgt.edge(2).Color = colsare(ifreq).c;
%                 
%                 lg(ifreq)=lgt.mainLine;
%                 %                 lg(idir)=plot(windowplot,mplot) ;
%                 hold on
%                 
%                 %          shg
%             end
%             set(gcf,'color', 'w')
%             lgt=legend([lg(1), lg(2), lg(3), lg(4), lg(5), lg(6)], {'8', '13', '15', '18', '20', '40'});
%             set(lgt, 'location', 'northwest')
%             set(gca, 'fontsize', 20)
%             xlabel('Time from Catch')
%             ylabel('dynSSVEP')
%             %             ylim([-2 2])
%             title(['DynSSVEP during ' ddir])
%             %     print snr
%             print('-dpng', ['figure_PFI_dynSSEP_' ddir '.png'])
%             
%             switch dtype
%                 case 1
%                     ppantPFI_dynSSVEP_0_1 = outgoingDATA;
%                 case 2
%                     ppantPFI_dynSSVEP_1_0 = outgoingDATA;
%                     
%                     
%             end
%             
%             
%             
%             
%         end
%         save('ppant_PFI_dynSSVEP', 'ppantPFI_dynSSVEP_0_1',...
%             'ppantPFI_dynSSVEP_1_0','timeid')
%     end
% end
% %%
% if job.concatacrossppants==1
%     cd(basefol);
%     peakfreqsare=[8,13,15,18,20,40];
%     
%     
%     
%     
%     icounter=1;
%     %
%     for id=1:2
%         storeall=zeros(length(allppants), 6, size(ppantPFI_dynSSVEP_0_1,2));
%         icounter=1;
%         for ifol=allppants
%             
%             cd(basefol)
%             cd(num2str(ifol));
%             load('ppant_PFI_dynSSVEP.mat')
%             
%             switch id
%                 case 1
%                     storeall(icounter,:,:) = ppantPFI_dynSSVEP_0_1;
%                 case 2
%                     storeall(icounter,:,:) = ppantPFI_dynSSVEP_1_0;
%                     
%             end
%             icounter=icounter+1;
%         end
%         
%         
%         switch id
%             case 1
%                 acrossp_PFI_dynSSVEP_0_1= storeall;
%             case 2
%                 acrossp_PFI_dynSSVEP_1_0= storeall;
%                 
%         end
%     end
%     
%     
%     cd(basefol)
%     cd('newplots-MD')
%     save('EEG_GrFX_PFI_dynSSVEP', 'acrossp_PFI_dynSSVEP_0_1','acrossp_PFI_dynSSVEP_1_0', 'timeid')
% end
% 
% 
% %%
% if job.plotPFIdynSSVEP==1;
%     cd(basefol)
%     cd('newplots-MD')
%     load('EEG_GrFX_PFI_dynSSVEP.mat')
%     
%     
%     
%     figure(1)
%     clf
%     hold on
%     ipl=1;
%     rlegend=[];
%     lg=[];
%     
%     linew=5;
%     colsare(1).c=[1 .1 .1];
%     colsare(2).c=[1 .9 0];
%     colsare(3).c=[0 .6 0];%'g';
%     colsare(4).c=[0 0 .6];%'b';
%     colsare(5).c=[0 0  0];%k';
%     colsare(6).c=[.5 .5 .5]; %grey
%     
%     for ihz=1:4%6%[1, 5, 6] %20Hz 40Hz
%         
%         
%         hzw=[num2str(peakfreqsare(ihz)) ' Hz'];
%         col=colsare(ihz).c;
%         
%         for idir=1:2
%             switch idir
%                 case 1 %show catch onset (removal of target)
%                     
%                     tmpd= squeeze(acrossp_PFI_dynSSVEP_0_1(:,ihz,:));
%                     
%                     
%                     lgwas = 'PFI Onset';
%                     linew=5;
%                     
%                     linst='-';
%                 case 2
%                     
%                     tmpd= squeeze(acrossp_PFI_dynSSVEP_1_0(:,ihz,:));
%                     lgwas = 'PFI offset';
%                     
%                     linst='-.';
%                     %                                         col='b';
%             end
%             
%             
%             
%             %for single subj. plotting
%             if size(tmpd,1)<size(tmpd,2)
%                 mplot = squeeze(nanmean(tmpd,1));
%                 
%                 
%                 %             adjust for withinsub errorbars.
%                 
%                 step1 = squeeze(nanmean(tmpd,2)); %within ppant mean (across timepoints);
%                 step2 = nanmean(step1); %total exp mean.
%                 
%                 %subtract ppant mean from ppant., and group mean from all
%                 newd = tmpd- repmat(step1, 1, size(tmpd,2)) - repmat(step2, size(tmpd));
%                 
%                 %             newd=tmpd;
%                 %calc new stErr
%                 stErr = nanstd(newd)/sqrt(size(newd,1));
%                 
%                 %plot
%                 %plot
%                 lgt=shadedErrorBar(timeid, mplot, stErr, 'k',[1]);
%                 lgt.mainLine.Color = colsare(ihz).c;
%                 lgt.mainLine.LineStyle = linst;
%                 lgt.patch.FaceColor = colsare(ihz).c;
%                 lgt.edge(1).Color = colsare(ihz).c;
%                 lgt.edge(2).Color = colsare(ihz).c;
%                 lg =[lg lgt.mainLine];
%                 hold on
%                 
%                 
%             else
%                 %                 mplot = tmpd';
%                 %                 lg(1) = plot(windowplot, mplot);
%             end
%             
%             ipl=ipl+1;
%             rlegend = [rlegend {[lgwas ' ' hzw ]}];
%             counter=counter+1;
%             
%         end
%     end
%     xlim([-3 3])
%     lga=legend([lg], rlegend);
%     ylim([-.4 .7])
%     set(lga, 'location', 'northwest')
%     xlabel('Time s')
%     ylabel({['dynamic SSVEP power'];['baseline normalised']})
%     title([''])
%     set(gca, 'fontsize', 15)
%     set(gcf, 'color', 'w');
% end
