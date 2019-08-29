% Follows the format of s3_CE... only now applying to RESS timeseries

%having epoched tgs and catch, POST RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.



clear all
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/MatthewDavidson/Desktop/FromIrenes Mac/Data Experiment copy')
end

basefol=pwd;
clearvars -except basefol allppants
dbstop if error




job.calcppantSNRperfreq=0; %sorts by hz x location. %topos in preivous script (s3_Da)
job.calcppantSNRperfreq_2fTG=0; %appends for 2fTGs/
job.concatacrossppanst=0;
job.plotHzsacrossppants=1; % this is averaging across all periods (excluding catch).

job.compareHzxLocations_bar=0;

getelocs;

cd(basefol)
cd('newplots-MD')
%% load data to determine physical catch timing.

%remaining participants after behavioral data analysis/exclusion
allppants=1:29;
badppants = [8, 15,28, 4 , 7,5 6,10];
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
tt = 0:1/srate:60;



if job.calcppantSNRperfreq==1
    %Collect all DATA
    %%
    windowsmall=[];
    windowsmall(1,:) = [-3 -.1] ;%  window targets present
    windowsmall(2,:) = [.1 3] ;%  window targets present after BP
%     windowsmall(1,:) = [-3 -1] ;%  window targets present
%     windowsmall(2,:) = [1 3] ;%  window targets present after BP
    %%
    
    
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.fpass= [0 50];
    param_spctrm.trialave=0;
    for ifol =allppants
        
        
        %%
        cd(basefol)
        cd(num2str(ifol))
        
        
        %% load the relevant PFI data.
        load('ppant_PFI_Epoched_RESS');
        load('ppant_Catch_Epoched_RESS')
        load('ppant_PFI_Epoched') %for info of types missing.
        load('TrialIndicesbyLocationandHz.mat')
        peakfreqsare=[8,13,15,18,20,40];
        %%
        
        RESS_staticSNR_byHzxLoc=zeros(6,4,205);
%         RESS_staticSNR_byHzxLoc=zeros(6,4,103);

        
        for ifreq=1:6
            
            
            switch ifreq
                case 1
                    usehz=8;
                case 2
                    usehz=13;
                case 3
                    usehz=15;
                case 4
                    usehz=18;
                case 5
                    usehz=20;
                case 6
                    usehz=40;
                    
            end
            
            for iloc=1:4
                
                
                
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
                
                %%
                for id=1:14% Use all epochs so as not to bias condition comparisons.
                    
                    switch id
                        case 1
                            dataIN=ress_PFI_0_1_HzxLoc;
                            searchtrials = [Freqwas.dir0_1.Trialind];
                            durscheck=durs0_1;
                            usewindow=1:2;
                        case 2
                            dataIN=ress_PFI_1_0_HzxLoc;
                            searchtrials = [Freqwas.dir1_0.Trialind];
                            durscheck=durs1_0;
                            usewindow=1:2;
                        case 3
                            dataIN=ress_PFI_1_2_HzxLoc;
                            searchtrials = [Freqwas.dir1_2.Trialind];
                            durscheck=durs1_2;
                            usewindow=1:2;
                        case 4
                            dataIN=ress_PFI_2_1_HzxLoc;
                            searchtrials = [Freqwas.dir2_1.Trialind];
                            durscheck=durs2_1;
                            usewindow=1:2;
                        case 5
                            dataIN=ress_PFI_3_2_HzxLoc;
                            searchtrials = [Freqwas.dir3_2.Trialind];
                            durscheck=durs3_2;
                            usewindow=1:2;
                        case 6
                            dataIN=ress_PFI_2_3_HzxLoc;
                            searchtrials = [Freqwas.dir2_3.Trialind];
                            durscheck=durs2_3;
                            usewindow=1:2;
                        case 7
                            dataIN=ress_BPcatchonsetTGs;
                            
                            usewindow=1;
                        case 8
                            dataIN=ress_BPcatchonsetBGs;
                            
                            usewindow=1;
                        case 9
                            dataIN=ress_catchonsetTGs;
                            
                            usewindow=1;
                        case 10
                            dataIN=ress_catchonsetBGs;
                            
                            usewindow=1;
                        case 11
                            dataIN=ress_catchoffsetTGs;
                            
                            usewindow=2;
                            
                        case 12
                            dataIN=ress_catchoffsetBGs;                            
                            usewindow=2;
                            
                        case 13
                            dataIN=ress_BPcatchoffsetTGs;                            
                            usewindow=2;
                        case 14
                            dataIN=ress_BPcatchoffsetBGs;                            
                            usewindow=2;
                            
                    end
                    
              
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    for windch=1:length(usewindow)
                        wind=usewindow(windch);
                        windcheck= windowsmall(wind,:);
                        
                        tidx=dsearchn(timeid', [windcheck]');
                        
                        
                        %reduce size. 
                        if id<7
                            datanow=squeeze(dataIN(ifreq,iloc).ressTS);
                            
                        datast= squeeze(datanow(:,tidx(1):tidx(2)));
                        else %catch trials
                            
                            datast= squeeze(dataIN(ifreq,iloc,:,tidx(1):tidx(2)));
                        end
                        
                        %% check for bad trials (noisy)
                        %std per trial(average over timepoints)
%                         
%                         datastSD = nanstd(datast,0,2);
%                         
%                         %remove those with 2.5*std from mean.
%                         trialSD=nanstd(datastSD);
%                         mSD=nanmean(datastSD);
                        keeptrials=1:size(datast,1);
%                         
%                         %remove the trials with excessive noise.
badtrials=[];
%                         badtrials=find(datastSD>mSD+2.5*trialSD)';
%                         
%                         % also skip the trials which were only transient button
%                         % presses. (less than one second).
                        shorttrials=[]; 
%                         
%                         %also remove NANs:
                        nantrials = find(isnan(datast(:,1)));
                        zerotrials= find(datast(:,1)==0);
                        badtrials = [badtrials, shorttrials, nantrials', zerotrials'];
%                         
%                         
%                         % remove these from consideration.
if sum(badtrials)>0
    keeptrials(badtrials)=[];
    datast=datast(keeptrials,:);
end
                        %%
%                         %%
                        
                        if id==1 && windch==1%start here
                            if ndims(datast)<2
                                dataOUT(1,:,:)=datast;
                            else
                                dataOUT=datast;
                            end
                            
                        else
                            % append
                            if ndims(datast)<2
                                tmp(1,:,:)=datast;
                                dataOUT=cat(1,dataOUT,tmp);
                                
                            else
                                dataOUT=cat(1,dataOUT,datast);
                            end
                            
                        end
                        
                    end
                    
                end
                
                %now we have ALL the data this freqxLoc.
                %%
                
                %
                    
                    [s,f]=mtspectrumc(dataOUT', param_spctrm)    ;
                    
                    %TAKE mean
                    sMEAN = squeeze(nanmean(log(s),2));
                    
                    %% comp SNR;
%                     kernelw = [-1/6 -1/6 -1/6   0 0 1  0 0 -1/6  -1/6 -1/6];
                    kernelw = [-1/4 -1/4 0 0 1  0 0 -1/4  -1/4 ];
                    %                 
%                     %take mean
                    sNOW= conv(sMEAN, kernelw, 'same');
                %
%                     figure(1); clf
%                     hold on;
%                     plot(f,sNOW);
%             
            %%
            
            
            RESS_staticSNR_byHzxLoc(ifreq,iloc,:) = sNOW; %channs x hz
            
            end
        end
        
        
        save('RESS_statSNR', 'RESS_staticSNR_byHzxLoc', 'f', 'kernelw')
    end
    
end

if job.calcppantSNRperfreq_2fTG==1
    %Collect all DATA
    %%
    windowsmall=[];
    windowsmall(1,:) = [-3 -.1] ;%  window targets present
    windowsmall(2,:) = [.1 3] ;%  window targets present after BP
%     windowsmall(1,:) = [-3 -1] ;%  window targets present
%     windowsmall(2,:) = [1 3] ;%  window targets present after BP
    %%
    
    
    param_spctrm.tapers = [1 1];
    param_spctrm.Fs= [250];
    param_spctrm.fpass= [0 50];
    param_spctrm.trialave=0;
    for ifol =allppants
        
        
        %%
        cd(basefol)
        cd(num2str(ifol))
        
        
        %% load the relevant PFI data.
        load('ppant_PFI_Epoched_RESS_2fTG');
        load('ppant_Catch_Epoched_RESS_2fTG')
        load('ppant_PFI_Epoched') %for info of types missing.
        load('TrialIndicesbyLocationandHz.mat')
        peakfreqsare=[8,13,15,18,20,40];
        %%
        
        RESS_staticSNR_byHzxLoc_2fTG=zeros(4,4,205);
%         RESS_staticSNR_byHzxLoc=zeros(6,4,103);

        
        for ifreq=1:4
            
            
            switch ifreq
                case 1
                    usehz=16;
                case 2
                    usehz=26;
                case 3
                    usehz=30;
                case 4
                    usehz=36;
                
            end
            
            for iloc=1:4
                
                
                
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
                
                %%
                for id=1:14% Use all epochs so as not to bias condition comparisons.
                    
                    switch id
                        case 1
                            dataIN=ress_PFI_0_1_HzxLoc;
                            searchtrials = [Freqwas.dir0_1.Trialind];
                            durscheck=durs0_1;
                            usewindow=1:2;
                        case 2
                            dataIN=ress_PFI_1_0_HzxLoc;
                            searchtrials = [Freqwas.dir1_0.Trialind];
                            durscheck=durs1_0;
                            usewindow=1:2;
                        case 3
                            dataIN=ress_PFI_1_2_HzxLoc;
                            searchtrials = [Freqwas.dir1_2.Trialind];
                            durscheck=durs1_2;
                            usewindow=1:2;
                        case 4
                            dataIN=ress_PFI_2_1_HzxLoc;
                            searchtrials = [Freqwas.dir2_1.Trialind];
                            durscheck=durs2_1;
                            usewindow=1:2;
                        case 5
                            dataIN=ress_PFI_3_2_HzxLoc;
                            searchtrials = [Freqwas.dir3_2.Trialind];
                            durscheck=durs3_2;
                            usewindow=1:2;
                        case 6
                            dataIN=ress_PFI_2_3_HzxLoc;
                            searchtrials = [Freqwas.dir2_3.Trialind];
                            durscheck=durs2_3;
                            usewindow=1:2;
                        case 7
                            dataIN=ress_BPcatchonsetTGs;
                            
                            usewindow=1;
                        case 8
                            dataIN=ress_BPcatchonsetBGs;
                            
                            usewindow=1;
                        case 9
                            dataIN=ress_catchonsetTGs;
                            
                            usewindow=1;
                        case 10
                            dataIN=ress_catchonsetBGs;
                            
                            usewindow=1;
                        case 11
                            dataIN=ress_catchoffsetTGs;
                            
                            usewindow=2;
                            
                        case 12
                            dataIN=ress_catchoffsetBGs;                            
                            usewindow=2;
                            
                        case 13
                            dataIN=ress_BPcatchoffsetTGs;                            
                            usewindow=2;
                        case 14
                            dataIN=ress_BPcatchoffsetBGs;                            
                            usewindow=2;
                            
                    end
                    
              
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    for windch=1:length(usewindow)
                        wind=usewindow(windch);
                        windcheck= windowsmall(wind,:);
                        
                        tidx=dsearchn(timeid', [windcheck]');
                        
                        
                        %reduce size. 
                        if id<7
                            datanow=squeeze(dataIN(ifreq,iloc).ressTS);
                            
                        datast= squeeze(datanow(:,tidx(1):tidx(2)));
                        else %catch trials
                            
                            datast= squeeze(dataIN(ifreq,iloc,:,tidx(1):tidx(2)));
                        end
                        
                        %% check for bad trials (noisy)
                        %std per trial(average over timepoints)
%                         
%                         datastSD = nanstd(datast,0,2);
%                         
%                         %remove those with 2.5*std from mean.
%                         trialSD=nanstd(datastSD);
%                         mSD=nanmean(datastSD);
                        keeptrials=1:size(datast,1);
%                         
%                         %remove the trials with excessive noise.
badtrials=[];
%                         badtrials=find(datastSD>mSD+2.5*trialSD)';
%                         
%                         % also skip the trials which were only transient button
%                         % presses. (less than one second).
                        shorttrials=[]; 
%                         
%                         %also remove NANs:
                        nantrials = find(isnan(datast(:,1)));
                        zerotrials= find(datast(:,1)==0);
                        badtrials = [badtrials, shorttrials, nantrials', zerotrials'];
%                         
%                         
%                         % remove these from consideration.
if sum(badtrials)>0
    keeptrials(badtrials)=[];
    datast=datast(keeptrials,:);
end
                        %%
%                         %%
                        
                        if id==1 && windch==1%start here
                            if ndims(datast)<2
                                dataOUT(1,:,:)=datast;
                            else
                                dataOUT=datast;
                            end
                            
                        else
                            % append
                            if ndims(datast)<2
                                tmp(1,:,:)=datast;
                                dataOUT=cat(1,dataOUT,tmp);
                                
                            else
                                dataOUT=cat(1,dataOUT,datast);
                            end
                            
                        end
                        
                    end
                    
                end
                
                %now we have ALL the data this freqxLoc.
                %%
                
                %
                    
                    [s,f]=mtspectrumc(dataOUT', param_spctrm)    ;
                    
                    %TAKE mean
                    sMEAN = squeeze(nanmean(log(s),2));
                    
                    %% comp SNR;
                    kernelw = [-1/6 -1/6 -1/6   0 0 1  0 0 -1/6  -1/6 -1/6];
                    kernelw = [-1/4 -1/4    0 0 1  0 0 -1/4 -1/4];
                    %                 
%                     %take mean
                    sNOW= conv(sMEAN, kernelw, 'same');
                %
%                     figure(1); clf
%                     hold on;
%                     plot(f,sNOW);
%             
            %%
            
            
            RESS_staticSNR_byHzxLoc_2fTG(ifreq,iloc,:) = sNOW; %channs x hz
            
            end
        end
        
        
        save('RESS_statSNR', 'RESS_staticSNR_byHzxLoc_2fTG', 'f', 'kernelw', '-append')
    end
    
end


if job.concatacrossppanst==1
    
    %%
    
    acrossRESS_SNR_freqs=zeros(length(allppants), 10,4,205);
    acrossRESS_SNR=zeros(length(allppants), 10,4);
    
    icounter=1;
    peakfreqsare=[8,13,15,18,20,40,16,26,30,36];
    for ippant=allppants
        cd(basefol)
        cd(num2str(ippant))
        
        load('RESS_statSNR')
        
        acrossRESS_SNR_freqs(icounter,1:6,:,:)=RESS_staticSNR_byHzxLoc;
        acrossRESS_SNR_freqs(icounter,7:10,:,:)=RESS_staticSNR_byHzxLoc_2fTG;
        
        for ifreq= 1:length(peakfreqsare)
            [~, id]= min(abs(f-peakfreqsare(ifreq)));
            if ifreq <7
            acrossRESS_SNR(icounter,ifreq,:)=squeeze(RESS_staticSNR_byHzxLoc(ifreq,:,id));
            else
                acrossRESS_SNR(icounter,ifreq,:)=squeeze(RESS_staticSNR_byHzxLoc_2fTG(ifreq-6,:,id));
            end
        end
        icounter=icounter+1;
    end
    %
    cd(basefol)
    cd('newplots-MD')
    save('GFX_ressSNR_static', 'acrossRESS_SNR_freqs','f', 'acrossRESS_SNR')
    
end



if job.plotHzsacrossppants==1
    %%
    cd(basefol)
    cd('newplots-MD')
    %load raw for comparison
    load('GFX_rawSNR_static.mat');
      load('GFX_ressSNR_static.mat'); 
      
    %%
    figure(1)
    clf
    
    cols(5).c='r';
    cols(6).c='r';
    
    peakfreqsare=[8,13,15,18,20,40,16,26,30,36];
    ressd=[];

    plotpeaks = [1,2,3,4,7,8,9,10];
    %%
    clf
    
    for ipl=1:2
        switch ipl
            case 1
                usehz=1:4;
                col= 'r';
                pname = '1st harmonic';
            case 2
                usehz=7:10;
                
                col= [0 .5 0];
                pname = '2nd harmonic';
                clf
            case 3
                usehz=5:6;
                
                col= [0 .5 0];
                
                pname = 'BGs';
%                 clf
        end
        pcounter=1;
    for ihz=usehz
        
        subplot(2,2,pcounter)
        if ihz==6            
            subplot(2,2,3)
            col='r';
        end
        % first plot the RESS result        
    meanHz=mean(squeeze(nanmean(acrossRESS_SNR_freqs(:,ihz,:,:),1)),1);
    
    
     ressd(ihz)=plot(f, meanHz, 'color', col, 'linewidth', 2);
    hold on
    
    xlim([ 5 20])   
    ylim([ -.25 .5])
    
    if  ihz>7
        xlim([20 40])
    elseif ipl==3
        xlim([peakfreqsare(ihz)-5 peakfreqsare(ihz)+5])
        ylim([ -1.5 5]) 
    end
    ylabel('log(SNR)')
    xlabel('Hz')
%     then plot the traditional
 usechan=30;
    meanHz=squeeze(mean(acrossRAWSNR_freqs(:,:,:,usechan,:),1));
    
    plotme=squeeze(mean(meanHz(1,:,:),2))';
    rawd=plot(f, plotme, 'k', 'linewidth', 3);
set(gca, 'fontsize', 15)
    legend(['RESS-SNR ' num2str(peakfreqsare(ihz)) ' Hz'])
    pcounter=pcounter+1;
%     axis square
    end
    set(gcf, 'color', 'w')
    %
    cd(basefol)
    cd('newplots-MD')
    cd('Figures')
    %
    print('-dpng', ['RESS target SNR, ' pname])
    end
     %%
  
   
   %% construct hz x loc for Raw
% now BAR equivalent to test for sig:
clf
acrossRAWSNR =[];
for ihz = 1:6
    for iloc = 1:4

 hzid = peakfreqsare(ihz);
        [~,hzid]=min(abs(f-hzid));
        
        acrossRAWSNR(:,ihz,iloc)=squeeze(acrossRAWSNR_freqs(:,ihz,iloc,usechan,hzid));
    end
end

%%
% construct data:
mRAW = squeeze(nanmean(acrossRAWSNR,1));
mRESS= squeeze(nanmean(acrossRESS_SNR,1));
tmpB=[];
tmpB(1,:,:)= mRAW;
tmpB(2,:,:)= mRESS;
clf
subplot(211)
bar(mRAW);% hold on; bar(mRESS)
subplot(212)
bar(mRESS)
shg
end

if job.compareHzxLocations_bar==1
 cd(basefol)
    cd('newplots-MD')
    %load raw for comparison
%     load('GFX_rawSNR_static.mat');
      load('GFX_ressSNR_static.mat'); 
      
    %%
    usefreqs=[1:4, 7:10]; %all TGs
% usefreqs=[1:4];

   storeall= squeeze(acrossRESS_SNR(:,[usefreqs],:));% TGs, ignoring BG.
    
   figure(1)   ; clf

    mbar= squeeze(nanmean(storeall,1));
    bh=bar(mbar') %use transpose to display with LOC on x axis
%     for ib=1:4
%     bh(ib).LineStyle=':';
%     bh(ib).LineWidth=4;
%     end
    ylabel('SNR(dB/Hz)')
%     xlabel('Hz')
%     legend TL TR BL BR
xlabel('Target Location')
% legend({'8 Hz', '13 Hz', '15 Hz', '18 Hz'})
legend({'8 Hz', '13 Hz', '15 Hz', '18 Hz', '16 Hz', '26 Hz', '30 Hz', '36 Hz'})
% legend({'16 Hz', '26 Hz', '30 Hz', '36 Hz'})
    shg
    %%
% end
%     storeall=
    X=storeall;
    pX=squeeze(mean(nanmean(storeall,2),3));
    meanOverall = mean(pX);
    
    tmpPX= ones(size(storeall,1),size(storeall,2), size(storeall,3));
    for i=1:size(pX)
        tmpPX(i,:,:) = tmpPX(i,:,:).*pX(i);
    end
    
    tmpGX= ones(size(storeall,1),size(storeall,2), size(storeall,3));
    adjustedST = X - tmpPX + (tmpGX.*meanOverall);
    %
    barstore=squeeze(nanstd(adjustedST))./sqrt(size(adjustedST,1));
    %
    barstore=barstore;
    hold on
    %
    numgroups = 4;
    numbars = size(storeall,2);%length(peakfreqsare);
    groupwidth = min(0.8, numbars/(numbars+1.5));
    % determine locations for errorbars:
    
    for i = 1:numbars
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
%         % use this row for Hz on x axis,%         
%                 errorbar(x, mbar(:,i), barstore(:,i), 'k', 'linestyle', 'none');
% use this row for LOC on x axis
        errorbar(x, mbar(i,:), barstore(i,:), 'k', 'linestyle', 'none');
    end
    %%
    %
%     title({['Pre and post ' dtypep ' onset']})
%     set(gca,'xticklabel',{'8 Hz', '13 Hz', '15 Hz' ,'18 Hz'}, 'fontsize', 20)
    set(gca,'xticklabel',{'Top Left', 'Top Right', 'Bottom Left' ,'Bottom Right'}, 'fontsize', 25)
    ylabel('SNR(dB/Hz)')

%      set(lga, 'location', 'northeast', 'fontsize', 15);
    set(gcf, 'color', 'w')
    ylim([ 0 6])
     shg
     cd(basefol)
     cd('newplots-MD')
     cd('figures')
     %%
     print('-dpng', 'ress stat-TG SNR by Hz x Loc')
     
     %% export for JASP
       %% export for stats analysis
    %
    TableforExport=table();
    thisentry=1; %initialize counter.
    %
    peakfreqsare=[8,13,15,18,16,26,30,36];
%     st
    for iloc = 1:4
        for ihz=1:length(peakfreqsare)
            
            %loc data =
            ppantdata = squeeze(storeall(:, ihz, iloc));
            
            %append horizontally to table.
            TableforExport=[TableforExport, table(ppantdata)];
            TableforExport.Properties.VariableNames(thisentry)={['LOC' num2str(iloc) '_' num2str(peakfreqsare(ihz)) 'Hz']};
            thisentry=thisentry+1;
        end
    end
    
    try cd('/Users/matthewdavidson/Desktop')
    catch
        cd('/Users/mattdavidson/Desktop')
    end
    
    writetable(TableforExport, 'RESS-TG hz x loc.csv')
end
