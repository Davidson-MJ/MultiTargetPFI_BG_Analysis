
clear all
try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
catch
    cd('/Users/MatthewDavidson/Desktop/FromIrenes Mac/Data Experiment copy')
end

basefol=pwd;

clearvars -except basefol allppants
dbstop if error



job.crunchandplotppantstatSNR=0; % per ppant, using BG freqs to check for changes.


job.concatacrossppants=0;

job.plotacrossppants=1; %paper figure?
job.exportforJASP=0; % for stats.

job.checkPFIsubtractCATCH=0;


allppants=1:29;
badppants = [8, 15,28, 4 , 7,5 6,10];
allppants(badppants)=[];

% window=[-2 4];
window=[-3 3];

srate=250;

if job.crunchandplotppantstatSNR==1
%sanity check
%% plot mean topo by location
% compare precatch and prePFI (same data vs RESS)
% clear all




printing=0;



for ifol=allppants;
    
    cd(basefol)
cd(num2str(ifol))
load('ppant_PFI_Epoched_RESS.mat')
load('ppant_Catch_Epoched_RESS.mat')
%%
% load('RESSfilterdata.mat');
load('TrialIndicesbyLocationandHz.mat');
%
%%
% which sections to compute static SNR over?
windowsmall=[];
windowsmall(1,:) = [-3 -0.1] ;% longer window targets present
windowsmall(2,:) = [.1 3] ;% short window targets absent

%how big is the whole epoch?
%(sue for finding the above sections).
window=[-3 3];

srate=250;

epochdur = sum(abs(window))*srate;

timeid = [0:1/srate:epochdur];
timeid= timeid-3;



% FFT parameters
nfft = ceil( srate/.1 ); % .1 Hz resolution
hz    = linspace(0,srate,nfft);
%%
mbarout=[];
mbarerr=[];
counter=1;
%
% clf
%


logscale=1; %0 for linear, 1 for log.

%%%% save and store both types
for ihz=1%:2
    
%for plotting.
    mbarout=[];
    mbarerr=[];

    switch ihz
        case 1
            ifreqt=5; %20Hz
            findhz=20;
        case 2
            ifreqt=6; %40Hz
            findhz=40;
    end
    
    for id=1:2 % data type = pfi or catch
        
        switch id
            case 1
                dis = ppant_SNREEG_PFI_0_1(:,usechans,:);
                trialswere=[Freqwas.dir0_1.Trialind];
                %
                durscheck= durs0_1;
                nPFI0_1=size(dis,1);
            case 2
                dis = ppant_SNREEG_catchramponset(:,usechans,:);
                durscheck=[];
        end
        
        
        
        % remove bad trials and very short.
        
        
        
        datastSD = nanstd(squeeze(nanmean(dis,2)),0,2);
        
        %remove those with 2.5*std from mean.
        trialSD=nanstd(datastSD);
        mSD=nanmean(datastSD);
        keeptrials=1:size(dis,1);
        
        %remove the trials with excessive noise.
        badtrials=find(datastSD>mSD+2.5*trialSD)';
        
        % also skip the trials which were only transient button
        % presses. (less than one second).
        shorttrials = find(durscheck<60);
        badtrials = [badtrials, shorttrials];
        
        % remove these from consideration.
        keeptrials(badtrials)=[];
        datast=dis(keeptrials,:,:);
        
        if id==1
        nPFI0_1=size(datast,1);
        end
        
        
        
        
        
        
        data=permute(datast,[2 3 1]);
        
        for prep=1:2 %with RESS, without.
%             if prep==1 %apply RESS
%                 %apply RESS to this trial data
%                 col='r';
%                 markert='o-';
%                 
%                 %% no specific location needed (since bg hz only).
%                 ressM=squeeze(ressEVEC_byHzxLoc(ifreqt,1,:));
%                 
%                 
%                 
%                 % reconstruct RESS component time series
%                 ress_ts1 = zeros(size(data,2),size(data,3));
%                 
%                 
%                 %%
%                 for ti=1:size(data,3)
%                     ress_ts1(:,ti) =ressM'*squeeze(data(:,:,ti));
%                 end
%             else
%                 col='k';
%                 %just use best chan (62)
%                 pchan = find(usechans==62);
%                 ress_ts1=squeeze(data(pchan, :,:));
%             end
            
            %
            % now take FFT of pre and post BP window.
            for iwin=1:2
                windcheck= windowsmall(iwin,:);
                
                tidx=dsearchn(timeid', [windcheck]');
                
                
                % OK, now with pre/post RESS data, take FFT of window:
                
                ressxT = (abs( fft(ress_ts1(tidx(1):tidx(2),:),nfft,1)/diff(tidx) ).^2);
                %
                
                if logscale==1
                    ressxT=10*log10(ressxT);
                end
                %takeSNR
                skipbinsHi=  5; % .5 Hz, hard-coded!
                skipbinsLo=5;
                numbinsLo  = 20+skipbinsLo; %  2 Hz, also hard-coded!
                numbinsHi=numbinsLo;
                 % compute SNR for this trial:
                
                [snrR,~] = deal(zeros(size(hz)));
                %%
                snrR=zeros(size(ressxT,2),length(hz));
                
                %careful when constructing SNR,
                
%                 skipbinsLo =  dsearchn(hz', [findhz-neighbour(ifreqt).snrlow(2)]'); % hard-coded in Hz, how many away before start shoulder.        
%         numbinsLo  = dsearchn(hz', [findhz-neighbour(ifreqt).snrlow(1)]'); % hard-coded in Hz, how many away before start shoulder.        
% %         
% %         
%         skipbinsHi =  dsearchn(hz', [neighbour(ifreqt).snrhigh(1)- findhz]'); % hard-coded in Hz, how many away before start shoulder.        
%         numbinsHi  = dsearchn(hz', [neighbour(ifreqt).snrhigh(2) - findhz]'); % hard-coded in Hz, how many away before start shoulder.        
% % % 
                
                if printing==1 %collect SNR per trial to plot SD.
                for itrial= 1:size(ressxT,2)
                    ressx=ressxT(:,itrial);
                    
                    for hzi=numbinsLo+1:length(hz)-numbinsHi-1
                        
                        numer = ressx(hzi);
                        
                        denom = nanmean(ressx([hzi-numbinsLo:hzi-skipbinsLo hzi+skipbinsHi:hzi+numbinsHi]) );
                        
                        if logscale==1
                        snrR(itrial,hzi) = numer- denom;
                        else
                        snrR(itrial,hzi) = numer./denom;
                        end
                    end
                end
                else %just take mean SNR.
                    ressx=squeeze(nanmean(ressxT,2));
                    for hzi=numbinsLo+1:length(hz)-numbinsHi-1
                        
                        numer = ressx(hzi);
                        
                        denom = nanmean(ressx([hzi-numbinsLo:hzi-skipbinsLo hzi+skipbinsHi:hzi+numbinsHi]) );
                        
                        if logscale==1
                        snrR(1,hzi) = numer- denom;
                        else
                        snrR(1,hzi) = numer./denom;
                        end
                    end
                    
                    
                    
                end
                %%
%                 figure(1);
%                 plot(hz, squeeze(nanmean(snrR,1)),[ col markert])
%                 hold on
%                 plot(hz, mean(ressxT,2), 'k')
%                 xlim([0 45])
%                 hold on
                
                [~, idH]= min(abs(hz-findhz));
                
                mbarout(id,prep,iwin)=nanmean(snrR(:,idH));
                mbarouterr(id,prep,iwin)=nanstd(snrR(:,idH));
                %             [~, idH]= min(abs(hz-40));
                %             mbar40= [mbar40,snrR(idH)];
                
                %             couner=counter+1;
            end
            
            
        end
    end
    switch ihz
        case 1
            mbar20=mbarout;
            mbar20err=mbarouterr;
        case 2
            mbar40=mbarout;
            mbar40err=mbarouterr;
    end
end
%%
if printing==1
for iplot=2:3
figure(iplot)
clf
switch iplot
    case 2
        mbar= mbar20;
        mbarerr=mbar20err;
        hzwas='20';
    case 3
        
        mbar= mbar40;
        mbarerr=mbar40err;
        hzwas='40';
end
%sort
plotbar= [mbar(1,1,1), mbar(1,2,1), mbar(2,1,1), mbar(2,2,1);
    mbar(1,1,2), mbar(1,2,2), mbar(2,1,2), mbar(2,2,2)];

errbar=[mbarerr(1,1,1), mbarerr(1,2,1), mbarerr(2,1,1), mbarerr(2,2,1);
    mbarerr(1,1,2), mbarerr(1,2,2), mbarerr(2,1,2), mbarerr(2,2,2)];

bar(plotbar)
hold on

numgroups = 2;
numbars =4;
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, plotbar(:,i), errbar(:,i), 'k', 'linestyle', 'none');
end


if logscale==1
ylabel('10log10(SNR)');
else    
ylabel('SNR');
end

%
lg=legend('PFI RESSsnr','PFI rawsnr', 'Catchtrial RESSsnr', 'Catchtrial rawsnr');
set(lg, 'location', 'northwest', 'fontsize', 15)
set(gca,'xticklabel', {'Targ present', 'Targ absent'})
title([hzwas 'Hz ppant ' num2str(ifol) ', M and SD, nPFI 0->1=' num2str(nPFI0_1)])
set(gcf, 'color', 'w')
set(gca, 'fontsize', 15)
ylim([-10 50])
if logscale==1
print('-dpng', ['figure_logsnr ' hzwas])
else
print('-dpng', ['figure_linearsnr ' hzwas])
end


end
end

dimsare = {'Datatype PFI/Catch', 'prep RESS/raw', 'window pre/post'};


save('RESS_logsnrPFIvsCATCH', 'mbar20', 'mbar20err', 'mbar40', 'mbar40err', 'dimsare')

end
end


if job.concatacrossppants==1;
   
    acrossp20=[];
    acrossp40=[];
    counter=1;
    
    logscale=1;
    
    for  ifol=allppants
        cd(basefol)
        cd(num2str(ifol))
        if logscale==1
        load('RESS_logsnrPFIvsCATCH.mat');
        else
        load('RESS_linearsnr.mat');
        end
        
        acrossp20(counter,:,:,:)= mbar20;
        acrossp40(counter,:,:,:)= mbar40;
        counter=counter+1;
    end
    %%
    cd(basefol)
    cd('newplots-MD')
    if logscale==1
        savename='GFX_logSNR';
    else
       savename='GFX_linearSNR';
    end
    
    save(savename, 'acrossp20', 'acrossp40', 'dimsare')
end
%%
if job.plotacrossppants==1; %paper figure?
    %%
    cd(basefol)
    cd('newplots-MD')
    logscale=1;
    if logscale==1
    load('GFX_logSNR');
    else
    load('GFX_linearSNR');
    end
    %%
    for ifig=1%1%:2
        figure(ifig)
        clf
        switch ifig
            case 1
                databar = acrossp20;
                hzwas='20';
            case 2
                databar = acrossp40;
                hzwas='40';
        end
        %%
        
        mbar = squeeze(nanmean(databar,1));
        mbarerr = squeeze(nanstd(databar,1))./sqrt(size(databar,1));
        %%
%         %dimsare = pfi/catch, ress/raw, pre/post
        plotbar= [mbar(1,1,1), mbar(1,1,2), mbar(1,2,1), mbar(1,2,2);
        mbar(2,1,1), mbar(2,1,2), mbar(2,2,1), mbar(2,2,2)];
       
errbar=[mbarerr(1,1,1), mbarerr(1,1,2), mbarerr(1,2,1), mbarerr(1,2,2);
    mbarerr(2,1,1), mbarerr(2,1,2), mbarerr(2,2,1), mbarerr(2,2,2)];
        
% 
% %dimsare = pfi/catch, ress/raw, pre/post
%         plotbar= [mbar(1,1,1), mbar(2,1,1); mbar(1,2,1), mbar(2,2,1)];        
%        
% errbar=[mbarerr(1,1,1), mbarerr(2,1,1); mbarerr(1,2,1), mbarerr(2,2,1)];
%     



bar(plotbar)
hold on

numgroups = 2;
numbars =4;
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
    errorbar(x, plotbar(:,i), errbar(:,i), 'k', 'linestyle', 'none');
end

if logscale==1
ylabel('10log10(SNR)')
else
ylabel('SNR')
end


%
lg=legend('RESS-TGpres','RESS-TGabent', 'Raw-TGpres', 'Raw-TGabsent');
% lg=legend('PFI-TGpres', 'Catch-TGpres');
set(gca,'xticklabel', {'RESS', 'raw'})
set(lg, 'location', 'northeast', 'fontsize', 15)
% set(gca,'xticklabel', {'PFI', 'Catch'})
set(gca,'xticklabel', {'RESS', 'raw'})
title([hzwas 'Hz acrossppants, M and Sterr'])
set(gcf, 'color', 'w')
set(gca, 'fontsize', 15)
ylim([5 20])
axis tight
shg
%%

if logscale==1
print('-dpng', ['figure_GFX logsnr ' hzwas ])    
else
% print('-dpng', ['figure_GFX linearsnr ' hzwas ])
end
    end
            
    
    
end

if job.exportforJASP==1; % for stats.


logscale=1;

  cd(basefol)
    cd('newplots-MD')
    %%
    if logscale==0
  load('GFX_linearSNR.mat')
    else
        load('GFX_logSNR.mat')
    end
        
  %%
    hzare = [20, 40];
       
    % we will export targ hz x Loc
    % ppants x Targ Loc
    
    %%
    TableforExport=table();
    thisentry=1; %initialize counter.
    %
    for id=1:2
        switch id
            case 1
                dtype='PFI';
            case 2
                dtype='Catch';
        end
    for iprep = 1:2
        switch iprep
            case 1
                ptype='RESS';
            case 2
                ptype='raw';
        end
        for iwin=1:2
            
        for ihz=1:2
            switch ihz
                case 1
            ppantdata = acrossp20(:,id,iprep,iwin);
            hzwas='20';
                case 2
            ppantdata = acrossp40(:,id,iprep,iwin);
            hzwas='40';
            end
            
            %append horizontally to table.
            TableforExport=[TableforExport, table(ppantdata)];
            TableforExport.Properties.VariableNames(thisentry)={[dtype '_' ptype '_' hzwas 'Hz_win' num2str(iwin)]};
            thisentry=thisentry+1;
        end
        end
    end
    end

    try cd('/Users/matthewdavidson/Desktop')
    catch
        cd('/Users/mattdavidson/Desktop')
    end
        
    writetable(TableforExport, 'PFIvsCATCHtable.csv')
end




if job.checkPFIsubtractCATCH==1
    %%
    cd(basefol)
    cd('newplots-MD')
    load('GFX_logSNR.mat')
    %         %dimsare = ppants, pfi/catch, ress/raw, pre/post
    % 
    diffPFICATCH= squeeze(acrossp20(:,1,:,1))-squeeze(acrossp20(:,2,:,1));
    
    
       plotbar= [nanmean(diffPFICATCH(:,1)) , nanmean(diffPFICATCH(:,2))];        
       
       errbar= [nanstd(diffPFICATCH(:,1))./sqrt(size(diffPFICATCH,1)) , nanstd(diffPFICATCH(:,2))./sqrt(size(diffPFICATCH,1))];        
    
       clf
bar(plotbar)
hold on
errorbar(1:2, plotbar, errbar)
    %%
set(gca, 'xticklabel', {'RESS', 'raw'}, 'fontsize', 25) 
title('PFI - catch')
ylabel('10log10SNR')
end



