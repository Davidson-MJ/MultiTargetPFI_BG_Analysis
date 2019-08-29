% s3_Da_computeRESSSNR_TargandCatch

% s3_CE_computeRAWSNR_TargandCatch

%having epoched tgs and catch, pre RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.



clear all
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/MatthewDavidson/Desktop/FromIrenes Mac/Data Experiment copy')
end

basefol=pwd;
clearvars -except basefol allppants
dbstop if error





job.concatTOPOacrossppanst=1;
job.plotTOPOsacrossppants=0;
job.plotTOPOsacrossppants_crunched=1;




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
tt = 0:1/srate:60;



if job.concatTOPOacrossppanst==1
    
    %%
    
    acrossRESSSNRTOPO=zeros(length(allppants), 10,4,64); %now 10 freqs f1 TGs x4, BG, 2BG, 2fTGsx4


    icounter=1;
    for ippant=allppants
        cd(basefol)
        cd(num2str(ippant))
        
        load('RESSfilterdata')
        acrossRESSSNRTOPO(icounter,1:6,:,:)= ressEVEC_byHzxLoc_MAPS;
        ressEVEC_byHzxLoc_MAPS=[];
        load('RESSfilterdata_2fTG')
        
        acrossRESSSNRTOPO(icounter,7:10,:,:)= ressEVEC_byHzxLoc_MAPS;
        
        icounter=icounter+1;
    end
    %
    cd(basefol)
    cd('newplots-MD')
    save('GFX_RESSmaps', 'acrossRESSSNRTOPO')
    
    
    
end
%%

if job.plotTOPOsacrossppants==1
    %%
    cd(basefol)
    cd('newplots-MD')
    load('GFX_RESSmaps.mat');
    %%
    figure(1)
    clf
    % normalize each topo for plotting.
    for ippant = 1:size(acrossRESSSNRTOPO,1)
        for iHz = 1:6%size(acrossRESSSNRTOPO,2)
            for iloc=1:size(acrossRESSSNRTOPO,3)
                acrossRESSSNRTOPO(ippant, iHz, iloc,:) = acrossRESSSNRTOPO(ippant, iHz, iloc,:)./ max(acrossRESSSNRTOPO(ippant, iHz, iloc,:));
            end
            
        end
        subplot(8,3,ippant)
        topoplot(squeeze(mean(acrossRESSSNRTOPO(ippant,4,:,:),3)), elocs(1:64));
        title(num2str(ippant))
    end
    %%
    figure(2);
    topoplot(squeeze(mean(acrossRESSSNRTOPO(7,:,iloc,:),2)), elocs(1:64));
    c=colorbar;
    ylabel(c, 'RESS filter weights')
    caxis([-.4 .4])
            set(gca, 'fontsize', 25)
    set(gcf, 'color', 'w'); 
    %%
    meanTOPOs=squeeze(mean(acrossRESSSNRTOPO,1));
    %%
    for it=1%:2
        figure(2);
        switch it
            case 1
                titlesare={'TL' , 'TR', 'BL', 'BR'};
            case 2                
                titlesare={'8 Hz' , '13 Hz', '15 Hz', '18 Hz'};
%                 titlesare={'16 Hz' , '26 Hz', '30 Hz', '36 Hz'};
        end
        for ipl=1:4
            subplot(2,2,ipl)
            if it==1
                topoplot(squeeze(mean(meanTOPOs(1:4,ipl,:),1)), elocs(1:64));
            else
                topoplot(squeeze(mean(meanTOPOs(ipl,:,:),2)), elocs(1:64));
            end
            c=colorbar;
            title(titlesare{ipl})
            ylabel(c, 'RESS filter weights')
            caxis([-.2 .2])
            set(gca, 'fontsize', 25)
        end
        
        set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('Figures')
%         %%
%         if it==1
%             print('-dpng', ['RESS TOPOs by target Locs.png'])
%         else
%             print('-dpng', ['RESS TOPOs by target frequencies.png'])
%         end
    end
    %%
    figure(2)
    clf
    titlesare={'20 Hz' , '40 Hz'};
    for ipl=1:2
        subplot(2,1,ipl)
        
        topoplot(squeeze(mean(meanTOPOs(ipl+4,:,:),2)), elocs(1:64), 'emarker2', {30, '*', 'w', 20,5});
        hold on
        topoplot(squeeze(mean(meanTOPOs(ipl+4,:,:),2)), elocs(1:64), 'emarker2', {30, 'o', 'k', 25, 5});
         c=colorbar;
        title(titlesare{ipl})
        ylabel(c, 'RESS filter weights')
        caxis([-.2 .2])
        set(gca, 'fontsize', 25)
    end
    set(gcf, 'color', 'w')
    %
    print('-dpng', ['RESS TOPOs by Background frequencies.png'])
    %
end

if job.plotTOPOsacrossppants_crunched==1
    %%
    cd(basefol)
    cd('newplots-MD')
    load('GFX_RESSmaps.mat');
    %%
    figure(1)
    clf
    % normalize each topo for plotting.
    for ippant = 1:size(acrossRESSSNRTOPO,1)
        for iHz = 1:6%size(acrossRESSSNRTOPO,2)
            for iloc=1:size(acrossRESSSNRTOPO,3)
                acrossRESSSNRTOPO(ippant, iHz, iloc,:) = acrossRESSSNRTOPO(ippant, iHz, iloc,:)./ max(acrossRESSSNRTOPO(ippant, iHz, iloc,:));
            end
            
        end
%         subplot(7,3,ippant)
%         topoplot(squeeze(mean(acrossRESSSNRTOPO(ippant,4,:,:),3)), elocs(1:64));
%         title(num2str(ippant))
    end
    %%
%     figure(2);
%     topoplot(squeeze(mean(acrossRESSSNRTOPO(7,:,iloc,:),2)), elocs(1:64));
%     c=colorbar;
%     ylabel(c, 'RESS filter weights')
%     caxis([-.4 .4])
%             set(gca, 'fontsize', 25)
%     set(gcf, 'color', 'w'); 
    %%
    meanTOPOs=squeeze(mean(acrossRESSSNRTOPO,1));
    %%
    %mean over LOC
    figure(2);
    clf
    for it=1:2
        
        switch it
            case 1
                % Tg freqs, all locs,
                dis = squeeze(mean(squeeze(meanTOPOs(1:4,:,:)),1));
                titleis={['Target flicker'];['1st harmonic']};
            case 2
                dis = squeeze(mean(squeeze(meanTOPOs(7:10,:,:)),1));
                titleis={['Target flicker'];['2nd harmonic']};
            case 3
                dis = squeeze(mean(squeeze(meanTOPOs(5,:,:)),1));                
                titleis={['Background flicker'];['1st harmonic']};
             case 4
                dis = squeeze(mean(squeeze(meanTOPOs(6,:,:)),1));                
                titleis={['Background flicker'];['2nd harmonic']};
            
        end
        %then mean over LOCS.
        dplot = squeeze(mean(dis,1));
        %plot the mean
        subplot(2,2,it);
        topoplot(dplot, elocs(1:64), 'conv', 'on')
        c=colorbar;
            title(titleis)
            ylabel(c, 'RESS filter weights')
            if it<3
                caxis([-.1 .1])
            else
            caxis([-.2 .2])
            end
            set(gca, 'fontsize', 20)
        end
        %%
        set(gcf, 'color', 'w')
        cd(basefol)
        cd('newplots-MD')
        cd('Figures')
%         %%
%         if it==1
%             print('-dpng', ['RESS TOPOs by target Locs.png'])
%         else
            print('-dpng', ['RESS TOPOs by Harmonics.png'])
%         end
    
    %%
    figure(2)
    clf
    titlesare={'20 Hz' , '40 Hz'};
    for ipl=1:2
        subplot(2,1,ipl)
        
        topoplot(squeeze(mean(meanTOPOs(ipl+4,:,:),2)), elocs(1:64), 'emarker2', {30, '*', 'w', 20,5});
        hold on
        topoplot(squeeze(mean(meanTOPOs(ipl+4,:,:),2)), elocs(1:64), 'emarker2', {30, 'o', 'k', 25, 5});
         c=colorbar;
        title(titlesare{ipl})
        ylabel(c, 'RESS filter weights')
        caxis([-.2 .2])
        set(gca, 'fontsize', 25)
    end
    set(gcf, 'color', 'w')
    %
    print('-dpng', ['RESS TOPOs by Background frequencies.png'])
    %
end
