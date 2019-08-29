% Follows the format of s3_CE... only now applying to RESS timeseries

%having epoched tgs and catch, POST RESS. apply FFT and SNR to windows when
%all targets are present/ away from Catch stimulus onset.




try cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
    addpath('/Users/MattDavidson/Desktop/SSVEP-PFI_BackgroundTag-master')
catch
    cd('/Users/mDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy')
end

basefol=pwd;
clearvars -except basefol allppants tgrm
dbstop if error

%%

%%%%% using nPFI based analysis from here:
job.calcppantDYNSNRperfreq=0; %sorts by hz x location. %topos in preivous script (s3_Da)

appendtradSNRtoDYN=0; % runs through the above, but performs traditional SNR on RESStimecourse for comparison.
job.concatacrossppanst=0;
job.plot_dynSSEP_acrossppants=1;

job.calcPpantDYNSSVEP_crosspoint=0;
job.plotPpant_crosspoints=0;
job.FirstSIGintimecourse=1; %no longer in use.

% separate analysis at image level.



%%%%% using image based analysis from here:
job.erpimagePpantlevel=0; %using dynRESS SSVEP or SNR 

job.concaterpdataacrossppants=0;

job.erpimageacrossppants=0;
job.gradedchangesinERPimage=0;

job.concaterpdataacrossppants_keepnumberSeparate=0; %trying to replicate plots, now separates by either side of median PFI


job.BPandSSVEPtimecourseacrossppants_numSeparate=0; 
job.BPandSSVEPtimecourseacrossppants_group=1; %same format as previous script (epoch PFI etc.)



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

epochdur = sum(abs(window));

timeidDYN = [0:1/srate:epochdur];
timeidDYN= timeidDYN-3;
%%
onsetc = ceil(epochdur)/2;
% peakfreqsare=[20,40]; %hz
%timing
tt = 0:1/srate:60;



if job.calcppantDYNSNRperfreq==1
    %Collect all DATA
    %%
    %%
    
    
    for ifol =allppants
        
        
        %%
        cd(basefol)
        cd(num2str(ifol))
        for iver=1%:2
            switch iver
                case 1
        peakfreqsare=[8,13,15,18,20,40];
                case 2
        peakfreqsare=[16,26,30,36];
            end
        
        %% load the relevant PFI data.
        if peakfreqsare(1)==16
        load('ppant_PFI_Epoched_RESS_2fTG');
        load('ppant_Catch_Epoched_RESS_2fTG');
        else
            load('ppant_PFI_Epoched_RESS');
        load('ppant_Catch_Epoched_RESS');
        end
        load('ppant_PFI_Epoched') %for info of types missing.
        load('TrialIndicesbyLocationandHz.mat')

        %%
        
        RESS_dynamicSNR_byTYPExHzxLoc=nan(10,length(peakfreqsare),4,1501);
        RESS_traditionalSNR_byTYPExHzxLoc=nan(10,length(peakfreqsare),4,33);
        
        
        for ifreq=1:length(peakfreqsare)
            
            
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
            
            if peakfreqsare(1)==16 %account for second harmonic.
                usehz=usehz*2;
%                 used = 6
            end
            for iloc=1:4
                
                
                
                % collect relevant trials for each type of spatial
                % configuration/filter construction.
                
                %%
                
                for id=1:10% Use all epochs so as not to bias condition comparisons.
                    
                    switch id
                        case 1
                            dataIN=ress_PFI_0_1_HzxLoc;
                            
                            dirsi= 1; %increasing PFI
                            
                        case 2
                            dataIN=ress_PFI_1_0_HzxLoc;
                           
                            dirsi= 0; %decre PFI
                            
                        case 3
                            dataIN=ress_PFI_1_2_HzxLoc;
                           
                            
                            dirsi=1;
                        case 4
                            dataIN=ress_PFI_2_1_HzxLoc;

                            
                            dirsi=0;
                        case 5
                            dataIN=ress_PFI_2_3_HzxLoc;
                            
                      
                            
                            dirsi=1;
                            
                        case 6
                            dataIN=ress_PFI_3_2_HzxLoc;
                          
                      
                            dirsi=0;
                            
                        case 7
                            if ifreq<5
                            dataIN=ress_BPcatchonsetTGs;
                            else
                                dataIN=ress_BPcatchonsetBGs;
                            end
                            usewindow=1;
                            dirsi=1;
                        
                        case 8
                            if ifreq<5
                            dataIN=ress_catchonsetTGs;
                            else
                                dataIN=ress_catchonsetBGs;                            
                            end
                            dirsi=1;
                            usewindow=1;
                        
                        case 9
                            if ifreq<5
                                dataIN=ress_BPcatchoffsetTGs;
                            else
                                dataIN=ress_BPcatchoffsetBGs;
                            end
                            dirsi=0;
                            usewindow=2;
                        
                        case 10
                            if ifreq<5
                                dataIN=ress_catchoffsetTGs;
                            else
                            dataIN=ress_catchoffsetBGs;    
                            end
                            dirsi=0;
                            usewindow=2;
                    end
                    
                    
                    % now that we have a datatype, get the right epochs that share stimulus configuration
                    
                    % also correct window (pre or post indication of target
                    % presence)
                    %reduce size.
                    if id<7
                        datast=squeeze(dataIN(ifreq,iloc).ressTS);
                        
                        
                    else %catch trials
                        
                        datast= squeeze(dataIN(ifreq,iloc,:,:));
                    end
                    
                    if size(datast,1)>1
                        %% check for bad trials (noisy)
                        %std per trial(average over timepoints)
                        %
                        datastSD = nanstd(datast,0,2);
                        %
                        %                     %remove those with 2.5*std from mean.
                        trialSD=nanstd(datastSD);
                        mSD=nanmean(datastSD);
                        keeptrials=1:size(datast,1);
                        %%
                        %remove the trials with excessive noise.
                        %                     badtrials=[];
                        badtrials=find(datastSD>mSD+2.5*trialSD)';
                        
                        % also skip the trials which were only transient button
                        % presses. (less than one second).
                        shorttrials=[];
                        
                        %also remove NANs:
                        nantrials = find(isnan(datast(:,1)));
                        zerotrials= [];%find(datast(:,1)==0);
                        badtrials = [badtrials, shorttrials, nantrials', zerotrials'];
                        
                        
                        % remove these from consideration.
                        keeptrials(badtrials)=[];
                        datast=datast(keeptrials,:);
                        %%
                                                                     
                        %now we have ALL the data this freqxLoc.
                        %%
                        if appendtradSNRtoDYN~=1 % then continue with
                            % dyn SSVEP (as per Cohen &
                            % Gulbinaite, 2017). or straight SNR.
                            
                            
                            %dynamic SSVEP given by abs(hilbert transform!
                            %first filter:
                            peakwidt=2;
                            %filter, but inlcude side lobes for temporal dynamics
                            fdatAt = filterFGx(datast,srate,usehz,peakwidt,0);
                            
                            
                            % take absolute of hilbert transform, to show dynamic SSEP
                            try dynSSEP= abs(hilbert(fdatAt'))';
                                
                                
                                %                 %TAKE mean per trial, careful with baseline subtraction
                                %                 timeid = [0:1/srate:epochdur];
                                %                 timeid= timeid-3;
                                %                 if dirsi==1 %increasing PFI / Target absence
                                %                     if ifreq<5 %targets, so low point post event.
                                %                         windcheck= [2 3];
                                %                     else %BG hz, so
                                %                         windcheck= [-3 -1];
                                %                     end
                                %                 elseif dirsi==0 %oppoite case, decreasing PFI/ targets are appearing
                                %                     if ifreq<5 %targets, so low point pre event.
                                %                         windcheck= [-3 -2];
                                %                     else %BG hz, so
                                %                         windcheck= [2 3];
                                %                     end
                                %                 end
                                %remove one second baseline
                                windcheck=[-2.9 -1.5];
                                %                 windcheck=[-3 3]; %whole epoch.
                                %
                                tidx=dsearchn(timeidDYN', [windcheck]');
                                %
                                %                %remove base
                                for itrial=1:size(dynSSEP,1)
                                    tmp= dynSSEP(itrial,:);
                                    mtmp = nanmean(tmp(tidx(1):tidx(2)));
                                    %                    mtmp = mean(tmp);
                                    
                                    dynSSEP(itrial,:)= tmp - repmat(mtmp, 1, length(tmp));
                                    
                                end
                                
                                
                                
                            catch
                                dynSSEP=nan(1,1501);
                            end
                            
                            
                            % store this type
                            RESS_dynamicSNR_byTYPExHzxLoc(id,ifreq,iloc,:) = squeeze(nanmean(dynSSEP,1));
                            
                        else % spectrogram of this RESS ts for comparison
                            %copied params from previous (EPOCH scripts).
                            
                            
                            
                            param_spcgrm.tapers = [1 1];
                            param_spcgrm.Fs= [250];
                            param_spcgrm.Fpass= [0 50];
                            param_spcgrm.trialave=0;

                                param_spcgrm.pad = 2;
                            
%                             movingwin=[2.5,.15];

                             movingwin=[1,.15];

                            [sgrm ,tgrm, fgrm] = mtspecgramc(datast', movingwin, param_spcgrm);
                            %%
                            %conv SNR
                            snr_sgrm =zeros(size(sgrm));
                            
                            %adjust for HBW
                            k = ((param_spcgrm.tapers(1,1)));%
                            hbw = (k+1)./ (2.*[movingwin(1,1)]);
                            neighb = 2; %hz
                            
                            %space out kernel appropriately
                            distz = dsearchn(fgrm', [hbw, neighb, neighb*2+hbw]');
                   
tmps= squeeze(nanmean(sgrm,3));
                                
                                
%                                 % SNR by convolution.
%    kernelw= [-1/4 -1/4  0 0 1 0 0 -1/4 -1/4 ];
%                             
%                                 snr_sgrm=[];
%                                 for itime= 1:size(tmps,1)
%                                     checkput = conv(log(tmps(itime,:)), kernelw,'same');
%                                     if ~isreal(checkput)
%                                         snr_sgrm(itime,:)= nan(1, size(tmps,2));
%                                     else
%                                         snr_sgrm(itime,:)= conv(log(tmps(itime,:)), kernelw,'same');
%                                     end
%                                 end


% 
% loop over frequencies and compute SNR
numbins = distz(2); 
skipbins = distz(1);
snr_sgrm= zeros(size(tmps));
for itime=1:size(tmps,1)
    for hzi=numbins+1:length(fgrm)-numbins-1
        numer = tmps(itime,hzi);
        denom = nanmean( tmps(itime,[hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snr_sgrm(itime,hzi) = numer./denom;
        
    end
end
                            
                          
                                %                             %store just stim freq.
                            [~, hzid]= min(abs(fgrm-usehz));
                            %
                            %reduce size.
                            mdynSSEP=squeeze(snr_sgrm(:,hzid,:));

                            RESS_traditionalSNR_byTYPExHzxLoc(id,ifreq,iloc,:) =mdynSSEP;
%                            

timeidGRM=tgrm-3;
%                             
                        end
                        
                        %% store mean per type per ppant.
                    else %else store nan if no PFI this type.
                        dynSSEP=datast;
                    end
                        
                

                end
            end
        end
        
        
        DIMSare = {'0_1', '1_0', '1_2', '2_1', '2_3', '3_2',...
            'BPonsets', 'catchonset', 'BPoffsets', 'catchoffset'};
         
        
        % check that the dims look right:
        if iver==1
        %%
%         clf
%             tmp1=squeeze(mean(RESS_traditionalSNR_byTYPExHzxLoc([1,3,5], 5,:,:),3)); %increase (m over locs)
%         tmp2=squeeze(mean(RESS_traditionalSNR_byTYPExHzxLoc([2,4,6], 5,:,:),3)); %decrease(m over locs)
%         
%          plot(tmp1', 'linestyle', '-'); hold on
%         plot(tmp2', 'linestyle', ':');
%         legend(DIMSare([1,3,5,2,4,6]))
        end
        
        
        RESS_traditionalSNR_byTYPExHzxLoc_snrMXC=RESS_traditionalSNR_byTYPExHzxLoc;
        
        
        %save whether traditional or dynamic amplitude, first or second
        %harmonic.
        
        if peakfreqsare(1)==16
            if appendtradSNRtoDYN~=1
                RESS_dynamicSNR_byTYPExHzxLoc_2fTG = RESS_dynamicSNR_byTYPExHzxLoc;
                save('RESS_dynamicSNR_byTYPExHzxLoc', 'RESS_dynamicSNR_byTYPExHzxLoc_2fTG', '-append')
            else
%                 RESS_traditionalSNR_byTYPExHzxLoc_2fTG = RESS_traditionalSNR_byTYPExHzxLoc;                
%                 save('RESS_dynamicSNR_byTYPExHzxLoc', 'RESS_traditionalSNR_byTYPExHzxLoc_2fTG', '-append')
RESS_traditionalSNR_byTYPExHzxLoc_snrMXC_2fTG = RESS_traditionalSNR_byTYPExHzxLoc;

                save('RESS_dynamicSNR_byTYPExHzxLoc', 'RESS_traditionalSNR_byTYPExHzxLoc_snrMXC_2fTG', '-append')
            end
        else
            if appendtradSNRtoDYN~=1
                save('RESS_dynamicSNR_byTYPExHzxLoc', 'RESS_dynamicSNR_byTYPExHzxLoc', 'timeidDYN', 'DIMSare')
            else
%                 save('RESS_dynamicSNR_byTYPExHzxLoc', 'RESS_traditionalSNR_byTYPExHzxLoc', 'timeidGRM', 'DIMSare','-append')
                    save('RESS_dynamicSNR_byTYPExHzxLoc', 'RESS_traditionalSNR_byTYPExHzxLoc_snrMXC', 'timeidGRM', 'DIMSare','-append')
            end
            
            
            
        end
        disp(['fin ppant ' num2str(ifol)]);
        end
    end
end


if job.concatacrossppanst==1
    
    %%
    
    acrossRESS_DYNSNR_freqs=zeros(length(allppants),10, 10,4,1501);
    acrossRESS_TRADSNR_freqs=zeros(length(allppants),10, 10,4,33);
%     [nppants, ndims, nfreqs,nlocs,nsamps]=size(acrossRESS_DYNSNR_freqs);
    icounter=1;
    for ippant=allppants
        cd(basefol)
        cd(num2str(ippant))
        
        load('RESS_dynamicSNR_byTYPExHzxLoc')
        
%         acrossRESS_DYNSNR_freqs(icounter,:,1:6,:,:)=RESS_dynamicSNR_byTYPExHzxLoc;
%         acrossRESS_DYNSNR_freqs(icounter,:,7:10,:,:)=RESS_dynamicSNR_byTYPExHzxLoc_2fTG(:,1:4,:,:);
        
        
        acrossRESS_TRADSNR_freqs(icounter,:,1:6,:,:)=RESS_traditionalSNR_byTYPExHzxLoc_snrMXC;
%         acrossRESS_TRADSNR_freqs(icounter,:,7:10,:,:)=RESS_traditionalSNR_byTYPExHzxLoc_2fTG(:,1:4,:,:);
        icounter=icounter+1;
    end
    %
    
    acrossRESS_TRADSNR_freqs_snrMXC=acrossRESS_TRADSNR_freqs;
    cd(basefol)
    cd('newplots-MD')
    
    save('GFX_ressSNR_Dynamic', 'acrossRESS_DYNSNR_freqs','acrossRESS_TRADSNR_freqs', 'timeidDYN', 'timeidGRM', 'DIMSare')
% save('GFX_ressSNR_Dynamic', 'acrossRESS_TRADSNR_freqs_snrMXC','-append')
    
end



if job.plot_dynSSEP_acrossppants==1
    %% using above format.
    cd(basefol)
    cd('newplots-MD')
    
    % ThiS SNR is average across subjects
    
    load('GFX_ressSNR_Dynamic.mat');
    %else we can load the average after normalization- but shouldn't
    %matter?
      
%     dbstop if error
    %
    
    useTRADorDYNSNR = 1; % plot traditional or dynamic SNR after RESS.
    
    interpCatch=1; % to plot with interpolated region for catch onsets
    
    switch useTRADorDYNSNR
        case 1
%             usedata=acrossRESS_TRADSNR_freqs_snrMXC;
            usedata=acrossRESS_TRADSNR_freqs;
            timeid=timeidGRM;
        case 2
            usedata=acrossRESS_DYNSNR_freqs;
            timeid=timeidDYN;
    end
    %
    figure(1)
    clf
    
    
    fontsize=15; 
    
    hold on
%     clf
    peakfreqsare= [8,13,15,18,20,40, 16,26,30,36];
    
    %different dimensions for different event types:
%     disp(DIMSare)

    PFIincrease = [1,3,5]; %
    PFIdecrease = [2,4,6];
%     PFIincrease = [5]; %
%     PFIdecrease =6;
    
    CatchONSET= 8;
    CatchOFFSET= 10;
    BPCatchONSET= 7;
    BPCatchOFFSET= 9;
    
    checkcluster=0; % temporal shuffle for sig.
    counter=1;
  
pl=[];
legendPRINT={};
       ttestdata=[];
       
       
       printOVERLAY=0; % set if want to print ontop (same subplot columns).
%        figure(1); clf
%          figure(2); clf
       %move ttest data if comparing types.
       ttestdata=[];
       
       
       storeBARSalso = []; % for extra interp.
       
       barcounter=1; %increase bar counter with each plot!
       
       for ihz=[2,3]%[2,3]%[1,4]%,4]%  2 3]% each HZ separately.
        
        mrks = '-';
          counter=1;
          
          itypech=[5,6];%,6]; %used for indexing ttests
          
          for itype=itypech
                
    
            switch itype
                case 1
                    EVENTdata = squeeze(nanmean(usedata(:,PFIincrease,:,:,:),2));
%                     EVENTdata = squeeze(usedata(:,PFIincrease,:,:,:));
%                     col = [.5 .7 .5];% dark green
                    
%                     col='r'
locis='NorthEast';
                    chtype = {['PFI Increase'];['Target Disappearance']};
                    xlabelis = 'Time from target invisible';
                    
                    mrks='-';
                    
                    %load BP data.
                    load('GFX_PFIperformance_withSNR_20_min0.mat', 'storeacrossPpant_onsetBP') 
                    useBP=storeacrossPpant_onsetBP;
                case 2
                    EVENTdata = squeeze(nanmean(usedata(:,PFIdecrease,:,:,:),2));
                    col='r';
                    chtype = {['PFI Decrease'];['Target Reappearance']};
                    locis='SouthEast';
                    xlabelis = 'Time from target visible';
                    load('GFX_PFIperformance_withSNR_20_min0.mat', 'storeacrossPpant_offsetBP') 
                    useBP=storeacrossPpant_offsetBP;
                    mrks=':';
                
                case 3
                    EVENTdata = squeeze(nanmean(usedata(:,CatchONSET,:,:,:),2));
                    col='g';
                    chtype = 'catch onset';
                    xlabelis = 'Time from PMD onset';
                    load('GFX_Catchperformance_withSNR_20.mat', 'storeacrossPpant_onsetBP') 
                    useBP=storeacrossPpant_onsetBP;
                case 4
                    EVENTdata = squeeze(nanmean(usedata(:,CatchOFFSET,:,:,:),2));
                    col='r';
                    chtype = 'catch offset';
                    xlabelis = 'Time from PMD offset';
                    load('GFX_Catchperformance_withSNR_20.mat', 'storeacrossPpant_offsetBP') 
                    useBP=storeacrossPpant_offsetBP;
                case 5
                    EVENTdata = squeeze(nanmean(usedata(:,BPCatchONSET,:,:,:),2));
                    col = [.5 .7 .5];
                    chtype = {['Catch'];['Target Disappearance']};
                    xlabelis = 'Time from reporting PMD onset';
                    mrks='-';
                    load('GFX_Catchperformance_withSNR_20,BPaligned.mat', 'storeacrossPpant_onsetBP') 
                    useBP=storeacrossPpant_onsetBP;
                case 6
                    EVENTdata = squeeze(nanmean(usedata(:,BPCatchOFFSET,:,:,:),2));
                    col='r';
%                     mrks=':';
                    chtype = {['Catch'];['Target Reappearance']};
                    xlabelis = 'Time from reporting PMD offset';
                    load('GFX_Catchperformance_withSNR_20,BPaligned.mat', 'storeacrossPpant_offsetBP') 
                    useBP=storeacrossPpant_offsetBP;
            end
            
            
               

            switch ihz
                case 1                    
            plotd= squeeze(nanmean(EVENTdata(:,1:4,:,:),3)); %mean over locations.
              plotd = squeeze(nanmean(plotd,2)); %mean over targets
%               titleis = {['Targets 1st harm.']};
                    col=['r']; 
                    
                    sigheight=.275;
                    ylimsare = [-.35 .55];
                case 2
                    plotd= squeeze(nanmean(EVENTdata(:,5,:,:),3)); %mean over locations.
% titleis = {['Background 1st harm.']};
              
              sigheight=1.15;
              ylimsare = [-.6 1];
              col='b';
%               mrks = '-';
                    case 3
                    plotd= squeeze(nanmean(EVENTdata(:,6,:,:),3)); %mean over locations.
% titleis = {['Background 2nd harm.']};
% sigheight=.75;
%                     mrks = '-';
col=[0 .5 0];
col='m';
sigheight=1.25;

                case 4
                                        
                    plotd= squeeze(nanmean(EVENTdata(:,7:10,:,:),3)); %mean over locations.
                    plotd = squeeze(nanmean(plotd,2));%mean over targets
col=[ 0 .5 0];
                                        mrks = '-';
%                     titleis = {['Targets 2nd harm.']};
            end
%               plotd = squeeze(nanmean(plotd,2));
            



%We can interpolate if using catch onSet data:
if (itype==3 || itype==5) && interpCatch==1 && (ihz ==3)
   %may need to interpolate at trial level. see what this looks like.
   %time vector 
   tvector=timeid;
   points= 1:length(tvector);
   %remove bad section from trace.
   %we used a one second sliding window, so:
   if itype==3
                             
                        badsec = dsearchn(tvector', [-.4 .4]');
   else
       badsec = dsearchn(tvector', [-1.45 -.36]');
   end
       
   tmp=points;
   tmp(badsec(1):badsec(2))=[];
   %input x vector with points missing.
   xIN=tmp;
   for iptrial=1:size(plotd,1)
      reald= plotd(iptrial,:);
      
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

plotd(iptrial,:)=reald2;

   end
    
end


%to place on the same figure!
if printOVERLAY==1 && mod(counter,2)==0
    counter=counter-1;
end

figure(1);
             subplot(2,2,counter)
            %plot BP.
            hold on
            tvector=[-3:1/60:3];
            
             ppantMeanBP = squeeze(mean(useBP,2));
             stE = std(ppantMeanBP)/sqrt(size(ppantMeanBP,1));
             sh=shadedErrorBar(tvector, mean(ppantMeanBP,1),stE,['-'],[1]);
              
             if itype==1 && printOVERLAY==1
             sh.mainLine.LineWidth=6;
             else
                 sh.mainLine.LineWidth=3;
             end
                    
              if printOVERLAY==11 && (itype==5 || itype==6)
                  colBP=[.3 .3 .3];
              else
                  colBP='r';
              end
                    sh.mainLine.Color = colBP;
                    sh.patch.FaceColor = colBP;
                    sh.edge(1).Color = colBP;
                    sh.edge(2).Color = colBP;
                    
             ylabel(['Buttons pressed'])
             set(gca, 'fontsize', fontsize)    
             
             xlim([-3 3])
             ylim([0 2.5])
             if printOVERLAY~=1
             xlabel(xlabelis)
             else
             xlabel('Time from subjective report')
             end
             
             
             
             
%%%%%%%%%         %%%%%%%%%       %%%%%%%%% 
%%%%%%%%%        %%%%%%%%%         %%%%%%%%%        NOW plot SNR
             subplot(2,2,counter+2)
             
            



            %SUBRTACT overall MEAN
%                       plotd=plotd-nanmean(plotd(:));
   
            %now plot %mean over ppants.
            mplot = squeeze(nanmean(plotd,1)); 
            
            
            %adjust errorbars.
            X=plotd;
            pX = nanmean(plotd,2);
            mX= nanmean(pX);
            
            newX = X - repmat(pX, [1, size(X,2)]) + repmat(mX, [size(X,1), size(X,2)]);
            
            stE = nanstd(newX,0,1)./sqrt(size(newX,1));
            
            sh=shadedErrorBar(timeid, mplot, stE, [ mrks], 1);
            
           
            %change plot properties if needed
            if printOVERLAY==11 &&  (itype==5  ||itype==6 )
                
                sh.mainLine.Color='k';
                sh.patch.FaceColor=col;
                sh.edge(1).Color='k';
                sh.edge(2).Color='k';
                
            else
                sh.mainLine.Color=col;
                sh.patch.FaceColor=col;
                sh.edge(1).Color=col;
                sh.edge(2).Color=col;
                sh.mainLine.LineWidth=3;
            end
            
            pl(itype)=sh.mainLine;
            hold on
            %
            set(gca, 'fontsize', fontsize)
            
            xlim([-3 3])
            
            if printOVERLAY~=1
                xlabel(xlabelis)
            else
                %              xlabel('Time from subjective report')
                xlabel('Time from PMD report')
            end
            
            
            ylabel({['RESS log(SNR)']})
            ttestdata(itype,:,:)=plotd;
            
             if ihz~=2 && ihz~=3
                 ylimsare=[0 .45];
             else
                 ylimsare=[1 2.2];
             end

%             ylim(ylimsare)
            
        
        if itype==3 && interpCatch~=1 % plot patch to cover frame slip
            yt=ylim;
            xV = [-.5 -.5 .5 .5]; %x coordinates of vertices.
            yV = [yt(1) yt(2) yt(2) yt(1)]; %y coordinates of vertices.
            pch=patch(xV,yV,[.5 .5 .5]);
            %                     pch.FaceAlpha =(.75);
            shg
            hold on
        end
%         ylim([-.05 .5])
%         xlim([-2.5 2.5])
axis tight

if ihz==2 || ihz ==3        
ylim([1.37 2.25])
else
    ylim([0 .5])
end
hold on
plot([0 0], ylim, ['k:']);

             
%%%%%%%%%         %%%%%%%%%       %%%%%%%%% 
%%%%%%%%%        %%%%%%%%%         %%%%%%%%%        NOW plot BAR CHARTS
         
            %                      storeBARSalso(counter,:,:) = plotd;
%             barlims=[-2.5 -.5 .5 2.5];
            barlims=[-2 -.1 .1 2];
            
            inspectimes = dsearchn(timeid', [barlims]');
            %shrink
            
            figure(2);
            if ihz==3 || ihz==4 
            subplot(2,2,counter+2)
            else
                subplot(2,2,counter)
            end
            %shrink data to timw window of interest
            bardata=[mean(plotd(:, inspectimes(1):inspectimes(2)),2),  mean(plotd(:, inspectimes(3):inspectimes(4)),2)];
            allbar = mean(bardata,1); %first plotted relationship
            stbar = std(bardata,1)./sqrt(length(allppants)); %first plotted relationship
            
            bh=bar(allbar); hold on;
            bh.FaceColor=col;
            
            errorbar(1:2, allbar, stbar,'linestyle', 'none', 'color', 'k');
            xlabel(xlabelis)
%             set(gca, 'xticklabels', {[num2str(barlims(1)) ':' num2str(barlims(2)) ' s'],[num2str(barlims(3)) ':' num2str(barlims(4)) ' s']})
            set(gca, 'xticklabels', {['before'], ['after']})
            ylabel('RESS log(SNR)')
                
            
            %check sig!
            [h,p,~,stats]= ttest(bardata(:,1), bardata(:,2));
            if p<.05;
%                 title(['significant, t= ' num2str(stats.tstat) ', p=', num2str(p)])
            end
            set(gca,'fontsize', fontsize)
            if ihz==2 %20 hz
                ylim([1.35 2])
            elseif ihz ==3 % 40 hz
                ylim([1.6 2.4])
            else
                ylim([0 .4])
            end
            
                
            
            
        counter=counter+1;
%           axis tight
          end            
       
     
          
%        lg=legend([pl(1) pl(2) ], legendPRINT);    
%          set(lg, 'Location','NorthEastOutside', 'fontsize',18)
        
%      check for significance from zero
%check for sig
if checkcluster==1;
    tvals=zeros(1,size(ttestdata,3));
    for itime = 1:size(ttestdata,3)
        
        try [h,pvals(itime),~,stat]=ttest(ttestdata(itypech(1),:,itime), ttestdata(itypech(2),:,itime));
            shuffType=1;
        catch
            %                 [h,pvals(itime),~,stat]=ttest(ttestdata(ihz,:,itime)); %compares to zero.
            shuffType=2; %whether or not to skip the non-parametric test for sig.
        end
        
        tvals(itime)= stat.tstat;
    end
    sigs=find(pvals<.05);
else
    sigs=[];
end
%         sigs=[];
%             %perform cluster based correction.
            if length(sigs)>2 && checkcluster==1
                % find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                %grab largest
%                 ignore bad points.
                
                
                  % find biggest cluster:
                %finds adjacent time points
                sigs = find(pvals<.05);
                
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
               
                %%
%                 sigheight=.45
                for icl=1:size(clusterSTandEND,1)

                    %start and end are now:
                    % change icl to maxClust if we only want the largest
                    % cluster.
                    STC=sigs(clusterSTandEND(icl,1));
                    ENDC=sigs(clusterSTandEND(icl,2)+1);
                    checktimes =STC:ENDC;
                    observedCV = sum(abs(tvals(checktimes)));
                    % now shuffle condition labels to see if this cluster is
                    % sig (compared to chance).
                    nshuffs=500;
                    sumTestStatsShuff = zeros(1,nshuffs);
                    
                    for irand = 1:nshuffs
                        %testing the null that it isn't mismatched - matched at time 2
                        % which creates a diff. so select from either!
                        shD= zeros(2,length(checktimes),size(ttestdata,2));
                        
                        %change shuffle parameters based on test of
                        %interest. (ie between conditions, or temporal
                        %null).
                        if shuffType==1 %null is that no condition differences.
                        for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    if mod(randi(100),2)==0 %if random even number
                                        pdata = ttestdata(itypech(1),randi(size(ttestdata,2)), checktimes(itime)); %select both chans
                                    else
                                        pdata = ttestdata(itypech(2),randi(size(ttestdata,2)), checktimes(itime));
                                    end
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                        else %null is that there are no temporal coincident sig values.
                            for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
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
                    
                        clf
                    
%                     subplot(2,1, plcount-1)
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
                          [~, c2] = min(abs(H.Data-observedCV)); %closest value. 
                        pvalis= 1-cdf(c2);
                    title(['sum tvals = ' num2str(observedCV), 'p=' num2str(pvalis)]);
%                     title(['sum tvals = ' num2str(observedCV)]);
                    %              title('Spatial Cluster  Significant!')
                    figure(1); hold on;
                    %space out the clusterfor plotting
                    if length(checktimes)>30
                    tryt= downsample(checktimes,30);
                    else
                        tryt=checktimes;
                    end
                    timeid(checktimes);
                    figure(1); hold on;
                        for itime=tryt
                            figure(useTRADorDYNSNR);
                            if ihz==3 && itime < 16 % skip the catch pre onset.
                            
                            else
                            hold on
%                             plot(timeid(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
                            plot(timeid(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', col)
                            end
                        end
                    end
            
%            
            end
            elseif length(sigs)>2
                
%                     tryt= downsample(sigs,30);
%                     
%                     
%                         for itime=tryt
%                             figure(1);
%                             hold on
%                             plot(timeidDYN(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'k')
% %                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
%                         end
            end
%    axis tight

    end
        
    %%
    
%     hold on; plot([0 0 ], ylim, ['k-'])
set(gcf, 'color', 'w')

%     xlim([timeid(1) timeid(end)])
cd(basefol)
cd('newplots-MD')
cd('Figures')

% cd('"X" over time')
%%

set(gcf, 'color', 'w')
%%
figure(1); set(gcf, 'color', 'w'); xlim([-2.5 2.3])
%%
print('-dpng', ['Dynamic SSVEP during catch '])
print('-dpng', ['Dynamic SSVEP during catch _bar'])
%     print('-dpng', ['Dynamic SSVEP during CATCH_BG.png'])

%     print('-dpng', ['Dynamic SSVEP during PFI_BG.png'])


end



if job.calcPpantDYNSSVEP_crosspoint==1
    
    
    
    
    cd(basefol)
    cd('newplots-MD')
    load('GFX_ressSNR_Dynamic.mat');
    
    usedata=acrossRESS_TRADSNR_freqs;
    timeid=timeidGRM;
    
    %
    figure(1)
    clf
    peakfreqsare= [8,13,15,18,20,40, 16,26,30,36];
    
    %different dimensions for different event types:
    %     disp(DIMSare)
    
    PFIincrease = [1,3,5]; %
    PFIdecrease = [2,4,6];
    CatchONSET= 8;
    CatchOFFSET= 10;
    BPCatchONSET= 7;
    BPCatchOFFSET= 9;
    
    checkcluster=0;
    counter=1;
    
    pl=[];
    legendPRINT={};
    ttestdata=[];
    
    
    printOVERLAY=1; % set if want to print ontop (same columns).
    clf
    storeBOTH_HZandPFIdirs=zeros(2,2,22,33);
    Xmarksthespot = zeros(2,22); %hz by ppants.
    
    typesch=[1,2];
    for itype=[typesch]
        
        switch itype
            case 1 % PFI Increase (target disappears)
                EVENTdata = squeeze(nanmean(usedata(:,PFIincrease,:,:,:),2));
               EVENTdata1= squeeze(usedata(:,PFIincrease,:,:,:));
                locis='NorthEast';
                chtype = {['PFI Increase'];['Target Disappearance']};
                xlabelis = 'Time from target invisible';
                
                
                
                %load BP data.
                load('GFX_PFIperformance_withSNR_20_min0.mat', 'storeacrossPpant_onsetBP')
                useBP=storeacrossPpant_onsetBP;
            case 2
                EVENTdata = squeeze(nanmean(usedata(:,PFIdecrease,:,:,:),2));
                col='r';
                chtype = {['PFI Decrease'];['Target Reappearance']};
                locis='SouthEast';
                xlabelis = 'Time from target visible';
                load('GFX_PFIperformance_withSNR_20_min0.mat', 'storeacrossPpant_offsetBP')
                useBP=storeacrossPpant_offsetBP;
                
            case 3
                EVENTdata = squeeze(nanmean(usedata(:,CatchONSET,:,:,:),2));
                col='g';
                chtype = 'catch onset';
                xlabelis = 'Time from catch onset';
                load('GFX_Catchperformance_withSNR_20.mat', 'storeacrossPpant_onsetBP')
                useBP=storeacrossPpant_onsetBP;
            case 4
                EVENTdata = squeeze(nanmean(usedata(:,CatchOFFSET,:,:,:),2));
                col='r';
                chtype = 'catch offset';
                xlabelis = 'Time from catch offset';
                load('GFX_Catchperformance_withSNR_20.mat', 'storeacrossPpant_offsetBP')
                useBP=storeacrossPpant_offsetBP;
            case 5
                EVENTdata = squeeze(nanmean(usedata(:,BPCatchONSET,:,:,:),2));
                col = [.5 .7 .5];
                chtype = {['Catch'];['Target Disappearance']};
                xlabelis = 'Time from reporting catch onset';
                
                load('GFX_Catchperformance_withSNR_20,BPaligned.mat', 'storeacrossPpant_onsetBP')
                useBP=storeacrossPpant_onsetBP;
            case 6
                EVENTdata = squeeze(nanmean(usedata(:,BPCatchOFFSET,:,:,:),2));
                col='r';
                chtype = {['Catch'];['Target Reappearance']};
                xlabelis = 'Time from reporting catch offset';
                load('GFX_Catchperformance_withSNR_20,BPaligned.mat', 'storeacrossPpant_offsetBP')
                useBP=storeacrossPpant_offsetBP;
        end
        
        
        counter=1;
        for ihz=[2,3]%  2 3]% each HZ separately.
            
            ttestdata=[];
            
            
            %               plotd= squeeze(nanmean(EVENTdata(:,ihz,:,:),3)); %mean over locations.
            switch ihz
                case 1 %targets
                    plotd= squeeze(nanmean(EVENTdata(:,1:4,:,:),3)); %mean over locations.
                    plotd = squeeze(nanmean(plotd,2)); %mean over targets
                    %               titleis = {['Targets 1st harm.']};
                    col=['r'];
                    mrks='-';
                    sigheight=.4;
                    ylimsare = [-.35 .55];
                case 2 %20 hz
                    plotd= squeeze(nanmean(EVENTdata(:,5,:,:),3)); %mean over locations.
                    % titleis = {['Background 1st harm.']};
                    mrks=':';
                    sigheight=.85;
                    ylimsare = [-.6 1];
                    col='r';
                case 3 %40 hz
                    plotd= squeeze(nanmean(EVENTdata(:,6,:,:),3)); %mean over locations.
                    % titleis = {['Background 2nd harm.']};
                    % sigheight=.75;
                    mrks = ':';
                    col=[0 .5 0];
                case 4 %targets 2nd harm.
                    
                    plotd= squeeze(nanmean(EVENTdata(:,7:10,:,:),3)); %mean over locations.
                    plotd = squeeze(nanmean(plotd,2));%mean over targets
                    col=[ 0 .5 0];
                    mrks = '-';
                    %                     titleis = {['Targets 2nd harm.']};
            end
            %               plotd = squeeze(nanmean(plotd,2));
            
            
            
            storeBOTH_HZandPFIdirs(counter,itype,:,:)= plotd;
            counter=counter+1;
            
        end
    end 
        
        % now we have all the data, calc "X" time per ppant, per hz.
        
        for ippant=1:size(storeBOTH_HZandPFIdirs,3)
           for ihz=1:2
               
               
            tmp1= squeeze(storeBOTH_HZandPFIdirs(ihz,typesch(1), ippant,:));
           tmp2=squeeze(storeBOTH_HZandPFIdirs(ihz,typesch(2),ippant,:));
           
           tmp1=tmp1-mean(tmp1);
           tmp2=tmp2-mean(tmp2);
           
           %find nearest SNR closest time point:
           %interp to find closest point (ms)
           samps = timeid(1):1/1000:timeid(end);
           
           yd1= pchip(timeid, tmp1,samps);
           yd2= pchip(timeid, tmp2,samps);
           
           % now find closest point
           diffV= abs(yd1-yd2);
           
           [~,Xmark] = find(diffV<.005);
           
           %include array endpoint so we don't skip last "X"
           Xmark= [Xmark, length(samps)];
           %restrict to independent "X" events.
           Xmarkclust=find(diff(Xmark)>1);
           
           Xmark=Xmark(Xmarkclust);
           
           falseX=1;
           Xcounter=1;
           skiptrial=0;
           % we want the first in a stretch:
           while falseX==1
               try
                   Xmarknow = Xmark(Xcounter);
                   
                   %test to see if this "X" continues.
                   %                % the difference for the next 2 seconds.
                   %                tmpD = yd1-yd2;
                   %                tmpDiff= tmpD(Xmarknow:Xmarknow+500);
                   %                % if we don't cross over again:
                   %                if min(tmpDiff)>=-.005 && mean(tmpDiff)>0 % protects against flips in wrong direction.
                   %                    falseX=0;
                   %                else
                   %                Xcounter=Xcounter+1;
                   %                end
                   %
                   %
                   
                   % or select
                   
                   %            sanity check
                   clf;
                   plot(samps, yd1, 'r'); hold on;
                   plot(samps, yd2,'k'); hold on;
                   legend('onset', 'offset')
                   title(num2str(ihz))
                   plot([samps(Xmark(Xcounter)),samps(Xmark(Xcounter))], ylim);
                   INP=input('Select this time?y/n','s');
                   
                   if strcmp(INP,'n')
                       Xcounter=Xcounter+1;
                   else
                       falseX=0;
                   end
                   
               catch % in case no "X" for this trial.
                   skiptrial=1;
                   falseX=0;
               end
               
               
               
               
               
           end
              
           
               if skiptrial~=1
           timeX = samps(Xmark(Xcounter));
           
           %            sanity check
%            clf;
%            plot(samps, yd1); hold on;
%            plot(samps, yd2); hold on;
%            plot([samps(Xmark(Xcounter)),samps(Xmark(Xcounter))], ylim);
%            
          

               else
                   timeX=nan;
               end

%store
Xmarksthespot(ihz, ippant) = timeX;
           
           end
            
            
            
        end
        %%
    MeanandSD=[nanmean(Xmarksthespot,2),nanstd(Xmarksthespot,0,2)];
    
    [h,p,~,tstat]=ttest(Xmarksthespot(1,:),Xmarksthespot(2,:));
        %%
    
%     Xmarksthespot_catch= Xmarksthespot;
    
    save('GFX_PFIdirCROSSPOINT', 'Xmarksthespot_catch', 'timeid','-append')
    
    %%%% some plotting:
    %%
    cd(basefol)
    cd('newplots-MD')
    load('GFX_PFIdirCROSSPOINT','timeid')
    %%
     figure(1); clf ;hold on
                
     Xmarksthespot=Xmarksthespot_catch        
%      Xmarksthespot=ppantFIRSTSIG;
    plot(timeid, ones(1, length(timeid)), 'color', 'w'); hold on; ylim([0 1])
    %
    %what are the means to plot?
    m1=(nanmean(Xmarksthespot(1,:),2)); 
    m2=(nanmean(Xmarksthespot(2,:),2)); 
    % SEM
    %adjust error bars.
    

    std1=nanstd(Xmarksthespot(1,:))/ sqrt(size(Xmarksthespot,1));
    std2=nanstd(Xmarksthespot(2,:))/ sqrt(size(Xmarksthespot,1));
    
%     std1=nanstd(Xmarksthespot(1,:))%/ sqrt(size(Xmarksthespot,1));
%     std2=nanstd(Xmarksthespot(2,:))%/ sqrt(size(Xmarksthespot,1));
    
    %
%      mXppant =squeeze( nanmean(Xmarksthespot,1)); %mean across conditions we are comparing (within ppant ie. hztypes).
%             mXgroup = nanmean(mXppant); %mean overall (remove b/w sub differences
%             
%             %for each observation, subjtract the subj average, add
%             %the group average.
%             NEWdata1 = m1x - mXppant + repmat(mXgroup, length(m1x),1);
%             NEWdata2 = m2x - mXppant + repmat(mXgroup, length(m2x),1);
%             
%             %             %compute new stErr %which version?
%             stE1 = nanstd(NEWdata1)/sqrt(size(m1x,1));
%             stE2 = nanstd(NEWdata2)/sqrt(size(m2x,1));
%     
%     std1=nanstd(Xmarksthespot(1,:))/ sqrt(size(Xmarksthespot,1));
%     std2=nanstd(Xmarksthespot(2,:))/ sqrt(size(Xmarksthespot,1));
    
%
% try plot all
% subplot(3,2,3)
clf
plot(timeid, ones(1, length(timeid)), 'color', 'w'); hold on; ylim([0 1])
 
for ihz=1:2
 h=histogram(Xmarksthespot(ihz,:),15); 
    
 switch ihz
        case 1
            colb=['b'];
        case 2
            colb=['m'];
 end
    h.FaceColor=colb;
    hold on; 
%     pl(plcounter)=plot([nanmean(ppantFIRSTSIG(ihz,:),2) nanmean(ppantFIRSTSIG(ihz,:),2)], [ylim], 'linew', 4, 'color', colb);
end
%
% subplot(3,2,5)
    %
    placeat=7.5;
    plot([m1 m1], [placeat placeat], 'marker', '*', 'linew', 400, 'color', 'b'); hold on
    plot([m2 m2], [placeat placeat], 'marker', '*', 'linew', 400, 'color', ['m']);
    hold on; %place errbars
    plot([m1-std1/2 m1+std1/2], [placeat placeat], 'marker', '.', 'linew', 1, 'color', ['b']);
    plot([m2-std2/2 m2+std2/2], [placeat placeat], 'marker', '.', 'linew', 1, 'color', ['m']);
        xlim([-2.5 2.5]);
        xlabel('Time from subjective report')
        xlabel('Time from catch report')
%         set(gca, 'fontsize', 25, 'ytick', [])
        set(gca, 'fontsize', 25)
        set(gcf, 'color', 'w')
    axis tight;
    ylabel('# of participants')
    xlim([timeid(1) timeid(end)])
    ylim([0 placeat*1.5])
%     subplot(3,2,6)
    
end


if job.FirstSIGintimecourse==1 %this uses
%     
    
    submean =0; % subtract mean from epochs for comparison.
    
    
    clf
    plcount=1;
    ppantFIRSTSIG=nan(2,22); %hz by ppants.
    ttestdata=[];
    
     
    icounter=1; %ppant counter.
    
    for ifol=allppants
        cd(basefol)
      
        cd(num2str(ifol))
        
    for hzis=1:2
        
        %%
        

        switch hzis
            case 1
                load('PFIperformance_withSNR_20_RESS.mat');
%                 load('Catchperformance_withSNR_20
                
            case 2
                load('PFIperformance_withSNR_40_RESS.mat');
                
    
        end
        
        % we likely have different trial numbers. so Downsample:
        offsettrials=size(Ppant_offsetSNR,1);
        onsettrials=size(Ppant_onsetSNR,1);
        ptmp1=[];
        ptmp2=[];
        if onsettrials >=offsettrials
            % sub-select with replacement to equate trial numbers
            ptmp=zeros(800,offsettrials,size(Ppant_onsetSNR,2));
            
            for i=1:800 % perform shuff for robustness
            
                %select offset trials (length), but up to max onset trials.
                dstrials1= randi([1,onsettrials], [1, offsettrials]);
                %select offset trials (length), but up to max offset trials.
                dstrials2 = randi([1, offsettrials], [1 offsettrials]);
            
            
            
%             ptmp(i,:,:)= Ppant_onsetSNR(dstrials,:);

            ptmp1(i,:,:)= Ppant_onsetSNR(dstrials1,:);
            ptmp2(i,:,:)= Ppant_offsetSNR(dstrials2,:);
            end
            
%             ppantONSET=squeeze(mean(ptmp,1));
%             ppantOFFSET=Ppant_offsetSNR;
            
            
            
            
        elseif onsettrials <offsettrials
             ptmp=zeros(800,onsettrials,size(Ppant_onsetSNR,2));
             for i=1:800 % perform shuff for robustness
                 dstrials = randi([1, onsettrials], [1 onsettrials]);
                 
                 
                 
                 
                %select offset trials (length), but up to max onset trials.
                dstrials1= randi([1,onsettrials], [1, onsettrials]);
                %select offset trials (length), but up to max offset trials.
                dstrials2 = randi([1, offsettrials], [1 onsettrials]);
            
            
            
                 
%                  ptmp(i,:,:)= Ppant_offsetSNR(dstrials,:);

                 ptmp1(i,:,:)= Ppant_onsetSNR(dstrials1,:);
                 ptmp2(i,:,:)= Ppant_offsetSNR(dstrials2,:);
             end
            
%             ppantOFFSET = squeeze(mean(ptmp,1));
%             ppantONSET=Ppant_onsetSNR;
            
            
        end
                    ppantONSET=squeeze(mean(ptmp1,1));
                        ppantOFFSET=squeeze(mean(ptmp2,1));
            
        
        %
if submean==1
    ppantONSET=ppantONSET-squeeze(mean(ppantONSET(:)));
    ppantOFFSET=ppantOFFSET-squeeze(mean(ppantOFFSET(:)));
end


  % we also need to make sure the SNR have crossed (!)
        %find nearest SNR closest time point:
           %interp to find closest point (ms)
           samps = timeidDYN(1):1/1000:timeidDYN(end);
           
           yd1= pchip(timeidDYN, squeeze(mean(ppantONSET,1)),samps);
           yd2= pchip(timeidDYN, squeeze(mean(ppantOFFSET,1)),samps);
           
           % now find closest point
           diffV= abs(yd1-yd2);
           
           [~,Xmark] = find(diffV<.005);
           
           Xmarkclust=find(diff(Xmark)>1);
           
           if ~isempty(Xmarkclust)
           
           Xmark=Xmark(Xmarkclust);
           else
               Xmark=Xmark(1);
           end
        
           crosspoints = samps(Xmark);
           crosspoint_first = dsearchn(timeidDYN', [crosspoints(1)']);
           if crosspoint_first==1
               crosspoint_first = dsearchn(timeidDYN', [crosspoints(2)']);
           end
%% plot for sanity check.
%                     clf
% %                     plot(timeidDYN,mean(ppantONSET,1), 'r'); hold on;
% %                     plot(timeidDYN,mean(ppantOFFSET,1)); hold on;
%                     
%                     sh=shadedErrorBar(timeidDYN, mean(ppantONSET,1), std(ppantONSET,1),[],1); hold on
%                     sh.mainLine.Color = 'b';
%                     sh.mainLine.LineWidth=3;
%                     sh.patch.FaceColor = 'b';
%                      
%                     sh=shadedErrorBar(timeidDYN, mean(ppantOFFSET,1), std(ppantONSET,1),[],1);
%                     sh.mainLine.LineStyle='--';
%                     sh.mainLine.Color = 'b';
%                     sh.mainLine.LineWidth=3;
%                     sh.patch.FaceColor = 'b';
%                     xlabel( 'Time from subjective report')
%                     set(gcf, 'color', 'w')
%                     set(gca, 'fontsize', 30);
%                     axis tight
% %                     hold on; %plot crosspoints
% %                     for ic=1:length(crosspoints)
% %                        ds= dsearchn(timeidDYN', crosspoints(ic)');
% %                        plot([timeidDYN(ds),timeidDYN(ds)], [1 1.5], 'b')
% %                         
% %                     end
%                     shg
                    %
        
        
        %now we perform ttest
            Tstats=[];
            pvals=[];
        for itime=1:size(ppantONSET,2)
            
            [~, pvals(itime),~,tstat]=ttest(ppantONSET(:,itime), ppantOFFSET(:,itime)) ;
            
            Tstats(itime)=tstat.tstat;
            
            % only look at correct direction
            diffNow = squeeze(mean(ppantONSET(:,itime))) -squeeze(mean(ppantOFFSET(:,itime)));
            
            % if the SNR is not correct direction, or pre cross-point,
            % ignore.
             if diffNow<0 || itime < crosspoint_first
                            pvals(itime)=nan;
             end
                 
        end
                
                
                sigs=find(pvals<.05);
                
% % plot these sigs.?
% ylimsare= get(gca, 'ylim');
% for itime = 1:length(sigs)
%     
%     plot(timeidDYN(sigs(itime)), [ylimsare(1)], '*k')
% end
   %%             
                
                
%             %perform cluster based correction.
            
%                 find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                
                
                sigtimes = timeidDYN(sigs);
                v1 = (vect1(:)==1);
                if sum(v1)>0
                    d = diff(v1);
                    clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                    
                    % use max cluster:
                    
%                     [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));


                        % or first!
                    [mr,~] = find((clusterSTandEND(:,2)-clusterSTandEND(:,1))>2);                        
                    maxClust=min(mr);
                    if isempty(maxClust)
                    [mr,~] = find((clusterSTandEND(:,2)-clusterSTandEND(:,1))>=1);                        
                    maxClust=min(mr);
                    end
                    
                    %start and end are now:                    
                    STC=sigs(clusterSTandEND(maxClust,1));
                    
                    
                    
                    if STC<1
                        errror('check')
                    end
                       
                    
                 
                    %%
%                     %%
                    hold on
                    plot([timeidDYN(STC) timeidDYN(STC)], ylim, ['k-'])
                    title(num2str(hzis))
                    if (STC)<20 % +.4 s
                    ppantFIRSTSIG(hzis,icounter)= timeidDYN(STC);
                    end
                    
                else %no sig cluster
                    if (sigs)>0
%                     ppantFIRSTSIG(hzis,icounter)= timeidDYN(min(sigs));
                    end
                end
%  ppantFIRSTSIG(hzis,itimezero,ippant)= min(sigs);
                
        end
%     disp(['chans used hz ' num2str(hzis) ',= ' num2str(usechans)])    
    
    icounter=icounter+1;
    end
    
   %% 
    figure(1);   
    clf
    plcounter=1;
    
        figure(1); hold on
        for ihz=1:2
            subplot(2,1,ihz)
            switch ihz
                case 1
                    
                    colb='r';
                    
                case 2
                    colb=[0 .5 0];
            end
    h=histogram(ppantFIRSTSIG(ihz,:),22); 
    h.FaceColor= colb;
    
    hold on; 
    pl(plcounter)=plot([nanmean(ppantFIRSTSIG(ihz,:),2) nanmean(ppantFIRSTSIG(ihz,:),2)], [ylim], 'linew', 4, 'color', colb);
    plcounter=plcounter+1;
                xlim([-3 3])
        
        
        [h,p,~,stat]= ttest(squeeze(ppantFIRSTSIG(1,:)), squeeze(ppantFIRSTSIG(2,:)));
        title(['t=' num2str(stat.tstat), 'p=' num2str(p)]);
        
    end
    %%
        figure(2); clf ;hold on
        
        tbase=tgrm-3;
        
        itimezero=2;
        
        
        
    plot(tbase, ones(1, length(tbase)), 'color', 'w'); hold on; ylim([0 1])
    m1=tbase(floor(nanmean(ppantFIRSTSIG(1,itimezero,:),3))); 
    m2=tbase(floor(nanmean(ppantFIRSTSIG(2,itimezero,:),3))); 
    
    %adjust error bars.
    m1x=squeeze(ppantFIRSTSIG(1,itimezero,:));
    m2x=squeeze(ppantFIRSTSIG(2,itimezero,:));
    
    
     mXppant =squeeze( nanmean(ppantFIRSTSIG(:,itimezero,:),1)); %mean across conditions we are comparing (within ppant ie. hztypes).
            mXgroup = nanmean(mXppant); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata1 = m1x - mXppant + repmat(mXgroup, length(m1x),1);
            NEWdata2 = m2x - mXppant + repmat(mXgroup, length(m2x),1);
            
            %             %compute new stErr %which version?
            stE1 = nanstd(NEWdata1)/sqrt(size(m1x,1));
            stE2 = nanstd(NEWdata2)/sqrt(size(m2x,1));
    
    std1=nanstd(ppantFIRSTSIG(1,itimezero,:))/ sqrt(21);
    std2=nanstd(ppantFIRSTSIG(2,itimezero,:))/ sqrt(21);
    
    plot([m1 m1], [.4 .4], 'marker', '*', 'linew', 400, 'color', 'r');
    plot([m2 m2], [.4 .4], 'marker', '*', 'linew', 400, 'color', 'g');
    hold on; %place errbars
    plot([m1-stE1/2 m1+stE1/2], [.4 .4], 'marker', '.', 'linew', 1, 'color', 'r');
    plot([m2-stE2/2 m2+stE2/2], [.4 .4], 'marker', '.', 'linew', 1, 'color', ['k']);
        xlim([-2.5 2.5]);
        xlabel('Time from PFI decrease')
        set(gca, 'fontsize', 25, 'ytick', [])
        set(gcf, 'color', 'w')
    %%
    
%     legend(pl, {'PFIinc 20', 'PFIinc 40', 'PFIdec 20', 'PFIdec 40'})
%     histogram(ppantFIRSTSIG(2,1,:));
shg
    %%
end


if job.erpimagePpantlevel==1
    %%
    rmvbase=0;
    
    useSNRorhilbert=1;
    for ifol = allppants
        for ihz=5:6; %TGs1 TGs 2, 20hz 40 hz.
            switch ihz                                
                case 5
                    usehz=20;
                    iloc=1;
                case 6
                    usehz=40;
                    iloc=1; %doesnt matter, always the same RESS filter.
            end
            
            
            
            icounter=1;
            
            cd(basefol)
            cd(num2str(ifol))
            
            load(['ppant_PFI_Epoched_RESS'])
            
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).

                    datatouse = cat(1, ress_PFI_0_1_HzxLoc(ihz,iloc).ressTS,ress_PFI_1_2_HzxLoc(ihz,iloc).ressTS,ress_PFI_2_3_HzxLoc(ihz,iloc).ressTS);
%                     RTstouse = [durs0_1'; durs1_2'; durs2_3'];
                    BPstouse = cat(1, ress_PFI_0_1_HzxLoc(ihz,iloc).BPs,ress_PFI_1_2_HzxLoc(ihz,iloc).BPs,ress_PFI_2_3_HzxLoc(ihz,iloc).BPs);
                    ctype = 'PFI increase';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = cat(1, ress_PFI_1_0_HzxLoc(ihz,iloc).ressTS,ress_PFI_2_1_HzxLoc(ihz,iloc).ressTS,ress_PFI_3_2_HzxLoc(ihz,iloc).ressTS);
%                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
                    BPstouse= cat(1, ress_PFI_1_0_HzxLoc(ihz,iloc).BPs,ress_PFI_2_1_HzxLoc(ihz,iloc).BPs,ress_PFI_3_2_HzxLoc(ihz,iloc).BPs);
                    ctype = 'PFI decrease';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                
                
                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
                tidx=dsearchn(timeidDYN', [windcheck]');
                
                
                
                %restrict to certain PFI duration?
                allt= 1:size(datatouse,1);
                %                 baddurs = find(RTstouse<30);
                %remove those trials.
                %                 allt(baddurs)=[];
                
                
                
               
                
                %%
                %%
                clf
                colormap('jet')
                %plot BP for comparison
                subplot(3,2, [1 3])
                %sort trial index by longest RT first.
                %                 [sortedRTs, cid] = sort(RTstouse(allt), 'descend');
                
                %sort by accum PFI in window after onset.
                if itimezero==1
                    accumPeriod = sum(BPstouse(:,180:end),2);
                else
                    accumPeriod = sum(BPstouse(:,1:180),2);
                    
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                % MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
                datais=squeeze(datatouse(cid,:));
%                 RTstouse=RTstouse(cid);
                %%

                tBP=-3:1/60:3;
                imagesc(-3:1/60:3,1:length(cid),BPstouse)
                c=colorbar;
                ylabel(c, 'buttons pressed')
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
                ylabel('PFI events')
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted BP data for ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15)
                %%
                % show mean over time.
                subplot(3,2,5)
                plot(tBP, nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('nanmean BP ')
                %
                xlabel('Time secs')
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                
                
                set(gca, 'fontsize', 15)
                %%
                

                if useSNRorhilbert~=1
                %% now plot dynamic RESS:
                
                peakwidt=2;
                 %filter, but inlcude side lobes for temporal dynamics
                fdatAt = filterFGx(datais,srate,usehz,peakwidt,0);                
                
                % take absolute of hilbert transform, to show dynamic SSEP
                dynSSEP= abs(hilbert(fdatAt'))';
                
                
%                 %TAKE mean per trial, careful with baseline subtraction
%                 timeid = [0:1/srate:epochdur];
%                 timeid= timeid-3;
%                 if dirsi==1 %increasing PFI / Target absence
%                     if ifreq<5 %targets, so low point post event.
%                         windcheck= [2 3];
%                     else %BG hz, so
%                         windcheck= [-3 -1];
%                     end
%                 elseif dirsi==0 %oppoite case, decreasing PFI/ targets are appearing
%                     if ifreq<5 %targets, so low point pre event.
%                         windcheck= [-3 -2];
%                     else %BG hz, so
%                         windcheck= [2 3];
%                     end
%                 end

%                 
                else % use SNR
                     
                            
                            param_spcgrm.tapers = [1 1];
                            param_spcgrm.Fs= [250];
                            param_spcgrm.Fpass= [0 50];
                            param_spcgrm.trialave=0;
                            movingwin=[1,.15];
                            
                            
                            [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                            %%
                            %conv SNR
                            snr_sgrm =zeros(size(sgrm));
                            
                            
                            
                            kernelw= [-1/6 -1/6 -1/6 0 0 1 0 0 -1/6 -1/6 -1/6 ];
                            
                                
                                
%                                 % SNR by convolution.
   kernelw= [-1/4 -1/4  0 0 1 0 0 -1/4 -1/4 ];
%                             
%                             
                            %snr on trials or
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
%                             
                            %store just stim freq.
                            [~, hzid]= min(abs(fgrm-usehz));
                            %
% %                             %reduce size.
                            dynSSEP=squeeze(snr_sgrm(:,hzid,:))';
% %                             
                          timeidDYN=tgrm-3;

                            
                    
                    
                end
%                         
%                %remove base
if rmvbase==1
    tidx=dsearchn(timeidDYN', [bsrem]');
    for itrial=1:size(dynSSEP,1)
        tmp= dynSSEP(itrial,:);
        mtmp = mean(tmp(tidx(1):tidx(2)));
        %                    mtmp = mean(tmp);
        
        dynSSEP(itrial,:)= tmp - repmat(mtmp, 1, length(tmp));
        
    end
    
end
                
               %
               
                %%\
                snrgrm20=dynSSEP;
                % smooth across trials
                sm_snrgrm20=zeros(size(snrgrm20));
                for itime=1:size(snrgrm20,2)
                    sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),5);
                end
                %%
                hold on
                subplot(3,2,[2 4]);
                hold on
                imagesc(timeidDYN,  1:size(dynSSEP,1), sm_snrgrm20);
                c=colorbar;
                if useSNRorhilbert==1
                ylabel(c, 'log(SNR)')
                else
                ylabel(c, 'abs(hilbert)')
                end
                
                xlabel('Time secs')
%                 caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])                %                 set(gca, 'ytick', 1:24, 'yticklabel', round(sortedRTs./60,2))
caxis([ -1 3])
                axis tight
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted ' num2str(usehz) 'Hz RESS data (smoothed) for  ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15, 'yticklabel', [])
                
                %         ylim([0 20])
                %% show mean over time.
                subplot(3,2,6)
                plot(timeidDYN, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
%                 ylim([1 3])
                hold on;
                ylabel('RESS logSNR')
                %%
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                %
                xlabel('Time secs')
                set(gcf,'color', 'w');
                
%                         ylim([-.1 10])
                shg
                %%
                cd([basefol filesep 'newplots-MD' filesep 'Figures' filesep 'Participant Summaries'])
                %%
                print('-dpng', ['figure_PFI_' ctype '_summary_BPandSSVEP_' num2str(usehz) '_RESS_subj' num2str(ifol) '.png']);
                
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        Ppant_onsetBP=BPstouse;
                        Ppant_onsetSNR=snrgrm20; %sorted.
%                         Ppant_onsetRTs=RTstouse;
%                         Ppant_onsetTOPO=snr20;
                    case 2
                        
                        Ppant_offsetBP=BPstouse;
                        Ppant_offsetSNR=snrgrm20; %always sorted in descending order of PFI.
%                         Ppant_offsetRTs=RTstouse;
%                         Ppant_offsetTOPO=snr20;
                end
            end
            
            
            %%
            switch ihz
                case 5
                    savename='PFIperformance_withSNR_20_RESS';
                case 6
                    savename='PFIperformance_withSNR_40_RESS';
            end
            
            
            save(savename,...
                'Ppant_onsetBP','Ppant_offsetBP',...
                'Ppant_onsetSNR', 'Ppant_offsetSNR', 'timeidDYN')
            
            
        end
    end
end



if job.concaterpdataacrossppants==1
    
    
    dursMINIMUM= 0; %1 second. % or change to  %this works.
%     dursMINIMUM=123; % selects top third!.
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=1:2
        
        if ihz==1
            loadname='PFIperformance_withSNR_20_RESS';
        else
            loadname='PFIperformance_withSNR_40_RESS';
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
%                 Ppant_onsetSNR=Ppant_onsetSNR';
%                 Ppant_offsetSNR=Ppant_offsetSNR';
                
                %restrict to only 'long' PFIs if needed.
%                 try
                    if dursMINIMUM>0 && dursMINIMUM~=123
                        
                        %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        shrton= find(Ppant_onsetRTs<dursMINIMUM); %note switch
                        all_onset(shrton)=[];
                        keepON=all_onset;
                        
                        %same for offsets.
                        shrtoff=find(Ppant_offsetRTs<dursMINIMUM);
                        all_offset(shrtoff)=[];
                        keepOFF=all_offset;
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
                    elseif dursMINIMUM==123
                        
                        %select top third of disap durations.
                         %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_onsetRTs, 'descend');
                        %top third
                        keepON = nid(1:ceil(length(ndurs)/4));
                        
                        %same for offsets.
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_offsetRTs, 'descend');
                        %top third
                        keepOFF = nid(1:ceil(length(ndurs)/3));
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
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
%                             pdataN(1:16,:) = repmat(pData(1,:), [16,1]);
                            pdataN(1:16,:) = flipud(pData(1:16,:));
                            %and bottom;
%                             pdataN(end-15:end,:) = repmat(pData(end,:), [16,1]);
                            pdataN(end-15:end,:) = flipud(pData(end-15:end,:));
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
                    %                 storeacrossPpant_onsetRTs(icounter,:)=Ppant_onsetRTs;
%                     storeacrossPpant_onsetTOPO(icounter,:)=Ppant_onsetTOPO';
                    
                    
                    storeacrossPpant_offsetBP(icounter,:,:)=Ppant_offsetBP;
                    storeacrossPpant_offsetSNR(icounter,:,:)=Ppant_offsetSNR;
                    %                 storeacrossPpant_offsetRTs(icounter,:)=Ppant_offsetRTs;
%                     storeacrossPpant_offsetTOPO(icounter,:)=Ppant_offsetTOPO';
                    
                    
                    icounter=icounter+1;
%                 catch SETUP
%                     if ippant~=27
%                     rethrow(SETUP)
%                     end
%                 end
        end
        
        %save appropriately
        cd(basefol)
        cd('newplots-MD')
        
        
        switch ihz
            case 1
                savename=['GFX_PFIperformance_withSNR_20_min' num2str(dursMINIMUM) '_RESS'];
            case 2
                savename=['GFX_PFIperformance_withSNR_40_min' num2str(dursMINIMUM) '_RESS'];
        end
        
%         timeidDYN=tgrm-3;
        save(savename,...
            'storeacrossPpant_onsetBP','storeacrossPpant_offsetBP',...
            'storeacrossPpant_onsetSNR', 'storeacrossPpant_offsetSNR', '-append')
        
        
        
    end
    
    
    
end
    





if job.erpimageacrossppants==1
    %% now plot across ppants.
    %reshape manually.
    
    %%
    getelocs
    rmvbase=0;
    
    
    counter=1;
    for hzis=1:2
        cd(basefol)
        cd('newplots-MD')
        switch hzis
            case 1
                
                load('GFX_PFIperformance_withSNR_20_min0_RESS')
                
                usehz=20;
%                 ylimsare=[1.55 2.1];
ylimsare=[1.6, 2.1];
            case 2
                
                load('GFX_PFIperformance_withSNR_40_min0_RESS')
                
                usehz=40;
%                 ylimsare=[1.8 2.7];
                ylimsare=[1.9, 2.5];
        end
%         tgrm = timeidDYN;
        %%
        for itimezero=1:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
%                     useRTs=storeacrossPpant_onsetRTs;
%                     useTOPO=storeacrossPpant_onsetTOPO;
                    
                    ctype= 'target invisible';
                    chtype='button press';
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
%                     useRTs=storeacrossPpant_offsetRTs;
%                     useTOPO=storeacrossPpant_offsetTOPO;
                    
                    ctype= 'target visible';
                    chtype='button release';
                    
                    
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
%             acrBP= squeeze(nanmean(useBP(21,:,:),1));
            
            acrSNR=squeeze(nanmean(useSNR,1));
%             acrSNR=squeeze(nanmean(useSNR(21,:,:),1));
            
            
            
            
            %sort across all
            %sort trial index by longest RT first.
            %%
            figure(1)
%             clf
%% not topoplotting anymore.             
%             usechan=62;
%             subplot(2,2,1:2)
%             topoplot(acrTOPO, elocs(1:64), 'emarker2', {usechan, 'o', 'w'});
%             c=colorbar;
%             caxis([0 15])
%             ylabel(c, {['10log10(SNR)'];['-3000:-100ms']})
%             hold on
%             title({[num2str(usehz) ' Hz SSVEP']})
%             set(gca, 'fontsize', 20)
            %%
            subplot(3,2,itimezero)
            colormap('jet')
            %         [sortedRTs, cid] = sort(acrRTs, 'descend');
            
            %
            imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
            %
            title({['Buttons Pressed'] }, 'fontsize', 20)
            c=colorbar;
            ylabel(c, 'Total')
            %                 ylabel('resampled catch trials')
            xlim([-2.5 2.5])
            set(gca, 'fontsize', 20, 'ytick', [])
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
            axis square
            %%
            xlabel(['Time from ' ctype])
            
            %SNR
            placeis= itimezero+2 + 2*(hzis-1);
            subplot(3,2, placeis);
            
            %reorder SNR
            acrSNRsort=acrSNR;
            
            %%
            %                 %smooth across trials
            sm_snrgrm20=zeros(size(acrSNRsort));
            for itime=1:size(acrSNRsort,2)
                sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
            end
            
            
            %     imagesc(acrSNR);
            imagesc(timeidDYN, 1:size(acrSNRsort,1), acrSNRsort);
            c=colorbar;

            ylabel(c, {['RESS log(SNR)']})

            title([num2str(usehz) ' Hz (' num2str(hzis) 'F)'])
%             title(['BG ' num2str(hzis) 'F'])
            
            hold on
            plot([0 0] , ylim, ['k:'], 'linewidth', 3)
            plot([0 0] , ylim, ['w:'], 'linewidth', 2)
            xlim([-2.5 2.5])
            set(gca, 'fontsize', 15, 'ytick', [])
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
            xlabel(['Time from ' ctype])
            ylabel({'Normalized'; 'trial count'})
                caxis([ylimsare])
                axis square
                
            set(gca, 'fontsize', 20)
            cd(basefol)
            cd('newplots-MD')
            cd('Figures')
            set(gcf, 'color', 'w')
            shg
            %%
           colormap('viridis')
            shg
            counter=counter+1;
        end
        
    end
    %%
    shg
    %%
     print('-dpng', ['PFI SNR summary BG Hz, ' ctype '.png'])
end



if job.gradedchangesinERPimage==1
    %% now plot across ppants.
    % as above, with some tweaks.
    
    %%
    
    clf
    %can plot BARS next to trial-by-trial image, or snr-timecourse.
for     plotBarsorTimeseries=1:2 % plots bars beneath, then timeseries.

    
    getelocs
    rmvbase=0;
    
    counter2=1; % for plotting threeway PFI
    
    
    
    counter=1;
    for hzis=1:2
        cd(basefol)
        cd('newplots-MD')
        switch hzis
            case 1
                
                load('GFX_PFIperformance_withSNR_20_min0_RESS')
                
                usehz=20;
%                 ylimsare=[1.55 2.1];
ylimsare=[1.7, 2.1];
col='b';
            case 2
                
                load('GFX_PFIperformance_withSNR_40_min0_RESS')
                
                usehz=40;
%                 ylimsare=[1.8 2.7];
                ylimsare=[2.1, 2.5];
                col='m';
        end
%         tgrm = timeidDYN;
        %%
        for itimezero=1:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
%                     useRTs=storeacrossPpant_onsetRTs;
%                     useTOPO=storeacrossPpant_onsetTOPO;
                    
                    ctype= 'target invisible';
                    chtype='button press';
                    
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
%                     useRTs=storeacrossPpant_offsetRTs;
%                     useTOPO=storeacrossPpant_offsetTOPO;
                    
                    ctype= 'target visible';
                    chtype='button release';
                    
                    
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
%             acrBP= squeeze(nanmean(useBP(21,:,:),1));
            
            acrSNR=squeeze(mean(useSNR,1));
%             acrSNR=squeeze(nanmean(useSNR(21,:,:),1));
            
            
            
            
            %sort across all
            %sort trial index by longest RT first.
            %%
            figure(1)
            hold on
            %%
            subplot(4,2,itimezero)
            colormap('jet')
            %         [sortedRTs, cid] = sort(acrRTs, 'descend');
            
            %
            imagesc([-3:1/60:3], 1:size(acrBP,1), acrBP);%(cid,:));
            %
%             title({['PFI ' num2str(ctype)]}, 'fontsize', 15)
            c=colorbar;
            ylabel(c, 'buttons pressed')
            %                 ylabel('resampled catch trials')
            xlim([-2.5 2.5])
            set(gca, 'fontsize', 15, 'ytick', [])
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
            axis square
            %%
            xlabel(['Time from ' ctype])
            
            %SNR
            if itimezero==1 
                placeis = hzis + 4;
            elseif itimezero==2
                placeis =  hzis + 6;
            end
            subplot(4,4, placeis);
            
            %reorder SNR
            acrSNRsort=acrSNR;
            
            %%
            %                 %smooth across trials
            sm_snrgrm20=zeros(size(acrSNRsort));
            for itime=1:size(acrSNRsort,2)
                sm_snrgrm20(:,itime)= smooth(acrSNRsort(:,itime),15);
            end
            
            
            %     imagesc(acrSNR);
            imagesc(timeidDYN, 1:size(acrSNRsort,1), acrSNRsort);
            c=colorbar;

%             ylabel(c, {['RESS log(SNR)']})

            title([num2str(usehz) ' Hz (' num2str(hzis) 'f)'])
%             title(['BG ' num2str(hzis) 'F'])
            
            hold on
            plot([0 0] , ylim, ['k:'], 'linewidth', 3)
            plot([0 0] , ylim, ['w:'], 'linewidth', 2)
            xlim([-2.5 2.5])
            set(gca, 'fontsize', 15, 'ytick', [])
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
            xlabel(['Time from ' ctype])
            ylabel({'Normalized'; 'trial count'})
                caxis([ylimsare])
                axis square
                
            set(gca, 'fontsize', 15)
            cd(basefol)
            cd('newplots-MD')
            cd('Figures')
            set(gcf, 'color', 'w')
            shg
            
            % beside this, plot the bar!
            %grouped nPFI
            %%
            groupednPFISNR=nan(length(allppants), 3);  % for barchart
            groupednPFISNR_ts=nan(length(allppants), 3, size(acrSNR,2));  % for SNR time-seriestrace
            if itimezero==1 % interested in after BP                
            BPidx=181:361;
            else
                BPidx=1:180;
            end
                groupedBEHflpd = fliplr(squeeze(nanmean(acrBP(:,BPidx),2))'); %takes mean over all SNR timepoints.
            %%
            groupednPFISNR_thirds=nan(length(allppants),3);
            
            for ippant = 1:size(useBP,1); 
%                % per ppant, collect trial indices:
                ppant_BEH= squeeze(mean(useBP(ippant,:,BPidx),3)); %takes mean over all BP time points.
                ppant_SNR=squeeze(mean(useSNR(ippant,:,:),3)); %takes mean over all SNR timepoints.
                ppant_SNR_ts=squeeze(useSNR(ippant,:,:)); %takes mean over all SNR timepoints.
               maxnPFI=ceil(max(ppant_BEH));
%                
               %flip BEH /SNR to correct orientation
               ppantBEHflpd=fliplr(ppant_BEH);
               ppantSNRflpd=fliplr(ppant_SNR);
               ppant_SNR_ts=flipud(ppant_SNR_ts);
               % only collect bin width.
%              
%                 nPFIidx = dsearchn(ppantBEHflpd', [0:3]');
               nPFIidx = dsearchn(groupedBEHflpd', [0:3]');
               
               nPFIidx=unique(nPFIidx);
               
               if nPFIidx(1) ~=1
                   nPFIidx(1)=1;
               end
               
               %take mean for these trials (actually rows in trialxtrial)
               for iPFI=1:length(nPFIidx)-1;
                   
                   try indices = nPFIidx(iPFI):nPFIidx(iPFI+1)-1;
                   catch %if max
                       indices = nPFIidx(iPFI):100;
                   end
                   
                   if iPFI==4;
                       shg
                   end
                   groupednPFISNR(ippant, iPFI)= squeeze(nanmean(ppantSNRflpd(indices)));
                   
                   
                   groupednPFISNR_ts(ippant, iPFI,:)= squeeze(nanmean(ppant_SNR_ts(indices,:)));
%                    %
%                    figure(3)
%                    clf
%                    plot(ppantSNRflpd);
%                    hold on
%                    plot([indices(1) indices(1)], ylim, ['k:'])
%                    plot([indices(end) indices(end)], ylim, ['k:'])
%                    title([num2str(squeeze(nanmean(ppantSNRflpd(indices))))])
               end
               
% works well.  

groupednPFISNR_thirds(ippant,1)= squeeze(nanmean(ppant_SNR(1:33)));
groupednPFISNR_thirds(ippant,2)= squeeze(nanmean(ppant_SNR(34:66)));
groupednPFISNR_thirds(ippant,3)= squeeze(nanmean(ppant_SNR(67:100)));

%                
            end
            %% using raw thirds?
%             groupednPFISNR_thirds=fliplr(groupednPFISNR_thirds);
%                         groupednPFISNR=fliplr(groupednPFISNR_thirds);
%             groupednPFISNR
            
            %%
            figure(1); 
            if plotBarsorTimeseries==1
   %place Bar charts beneath SNR
            if itimezero==1 
                placeis = hzis + 8;
            elseif itimezero==2
                placeis =  hzis + 10;
            end
         mbar=(nanmean(groupednPFISNR,1));
         
         
         %adjust errorbars:
          %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = groupednPFISNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant)
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
        stB=nanstd(NEWdata,1)./sqrt(length(allppants));
        
%             stB=nanstd(groupednPFISNR,1)./sqrt(length(allppants));
        
subplot(4,4, placeis)

        b1=bar(mbar);        
         b1.FaceColor='b';
         hold on         
         errorbar(mbar, stB, 'color', 'k','LineStyle', 'none')
         
set(gca, 'fontsize', 15, 'xticklabels', {'0\leq1', '1\leq2', '2\leq3 or 4'})
xlabel('amount of PFI');
         ylabel('RESS log(SNR)')
         
         colormap('viridis')
ylim(ylimsare)            
shg
            title([num2str(usehz) ' Hz (' num2str(hzis) 'f)'])

            counter=counter+1;
            
            nppants=length(allppants);
            nconds = size(groupednPFISNR,2);
            %RMANOVA?
            dataan= reshape(groupednPFISNR, [nppants*nconds,1]);
            conds = [ones(nppants,1); ones(nppants,1)*2;ones(nppants,1)*3];
                subs = repmat([1:length(allppants)]', [nconds,1]);
            
            [anret]=rmanova(dataan, conds, subs);
%%
            %LME to be sure:            
          tblA=table(dataan, conds, subs);
            
             UNterm1='dataan ~  conds + (1|subs)';
%              UNterm2='dataan ~  conds';
             UNterm2='dataan ~  1+ (1|subs)';
             %creating models:
             lmeUN =fitlme(tblA, UNterm1);
             lmeUN2 =fitlme(tblA, UNterm2);
             
    LMEout=compare(lmeUN2, lmeUN);
            %%
% %%
            else % plot time series instead.
%place Bar charts beneath SNR
            if itimezero==1 
                placeis = hzis + 12;
            elseif itimezero==2
                placeis =  hzis + 14;
            end
            subplot(4,4,placeis);
                tbase = timeidDYN;
    lgsav=[];
% can also plot time-series:
for iPFI=1:3
    x= squeeze(groupednPFISNR_ts(:,iPFI,:));
    %adjust errorbars.
    
       
            mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = nanmean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            stE = nanstd(NEWdata)/sqrt(size(x,1));
            %
            hold on
            sh=shadedErrorBar(tbase, nanmean(x,1),stE,['-'],[1]);
                
                sh.mainLine.Color = col;
                sh.mainLine.LineWidth = iPFI*2;
                sh.patch.FaceColor = col;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                
   lgsav(iPFI)= sh.mainLine; 
end
%%
% legend([lgsav], [{'nPFI 0 \leq 1'}, {'nPFI 1 \leq 2'}, {'nPFI 2 \leq 3 or 4'}])
%%
set(gca,'fontsize', 15); set(gcf, 'color', 'w')
xlabel(['time from ' chtype])
ylabel('RESS logSNR')
axis tight
ylim([ylimsare(1)-.2 ylimsare(2)+.4])
            end
%%
counter2=counter2+1;
    
            %%
           
       
        end 
    end
end
    %%
    shg
    %%
    colormap('viridis')
    shg
%      print('-dpng', ['PFI SNR summary BG Hz, ' ctype '.png'])
end




if job.erpimagePpantlevel==1
    %%
    rmvbase=0;
    
    useSNRorhilbert=1;
    for ifol = allppants
        for ihz=5:6
            switch ihz
                case 5
                    usehz=20;
                    iloc=1;
                case 6
                    usehz=40;
                    iloc=1; %doesnt matter, always the same RESS filter.
            end
            
            
            
            icounter=1;
            
            cd(basefol)
            cd(num2str(ifol))
            
            load(['ppant_PFI_Epoched_RESS'])
            
            for itimezero = 1:2
                if itimezero==1
                    
                    %append them all. TARG-> Disappearing (more buttons
                    %pressed).

                    datatouse = cat(1, ress_PFI_0_1_HzxLoc(ihz,iloc).ressTS,ress_PFI_1_2_HzxLoc(ihz,iloc).ressTS,ress_PFI_2_3_HzxLoc(ihz,iloc).ressTS);
%                     RTstouse = [durs0_1'; durs1_2'; durs2_3'];
                    BPstouse = cat(1, ress_PFI_0_1_HzxLoc(ihz,iloc).BPs,ress_PFI_1_2_HzxLoc(ihz,iloc).BPs,ress_PFI_2_3_HzxLoc(ihz,iloc).BPs);
                    ctype = 'PFI increase';
                    bsrem = [-3 -1]; %seconds
                    
                else
                    
                    datatouse = cat(1, ress_PFI_1_0_HzxLoc(ihz,iloc).ressTS,ress_PFI_2_1_HzxLoc(ihz,iloc).ressTS,ress_PFI_3_2_HzxLoc(ihz,iloc).ressTS);
%                     RTstouse = [durs1_0'; durs2_1'; durs3_2'];
                    BPstouse= cat(1, ress_PFI_1_0_HzxLoc(ihz,iloc).BPs,ress_PFI_2_1_HzxLoc(ihz,iloc).BPs,ress_PFI_3_2_HzxLoc(ihz,iloc).BPs);
                    ctype = 'PFI decrease';
                    
                    
                    bsrem = [1 3]; %seconds
                    
                    
                end
                
                
                
                %plot the topo for sanity check: (pre disap)
                windcheck= [-3 -0.1];
                tidx=dsearchn(timeidDYN', [windcheck]');
                
                
                
                %restrict to certain PFI duration?
                allt= 1:size(datatouse,1);
                %                 baddurs = find(RTstouse<30);
                %remove those trials.
                %                 allt(baddurs)=[];
                
                
                
               
                
                %%
                %%
                clf
                colormap('jet')
                %plot BP for comparison
                subplot(3,2, [1 3])
                %sort trial index by longest RT first.
                %                 [sortedRTs, cid] = sort(RTstouse(allt), 'descend');
                
                %sort by accum PFI in window after onset.
                if itimezero==1
                    accumPeriod = sum(BPstouse(:,180:end),2);
                else
                    accumPeriod = sum(BPstouse(:,1:180),2);
                    
                end
                
                [checkBP, cid] = sort(accumPeriod, 'descend');
                
                
                % MOVE nan to bottom of order (not top)
                nanind= find(isnan(checkBP));
                %new start point
                cid=[cid(length(nanind)+1:end) ; cid(nanind)];
                %                     checkBP=[checkBP(length(nanind)+1:end); checkBP(nanind)]
                
                %rearrange.
                BPstouse=BPstouse(cid,:);
                datais=squeeze(datatouse(cid,:));
%                 RTstouse=RTstouse(cid);
                %
                
                tBP=-3:1/60:3;
                imagesc(-3:1/60:3,1:length(cid),BPstouse)
                c=colorbar;
                ylabel(c, 'buttons pressed')
                %                 set(gca, 'ytick', 1:length(cid), 'yticklabel', round(sortedRTs./60,2))
                ylabel('PFI events')
                xlabel('Time [sec]')
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted BP data for ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15)
                
                % show mean over time.
                subplot(3,2,5)
                plot(tBP, nanmean(BPstouse,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('nanmean BP ')
                %
                xlabel('Time secs')
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                
                
                set(gca, 'fontsize', 15)
                
                
                if useSNRorhilbert~=1
                %% now plot dynamic RESS:
                
                peakwidt=2;
                 %filter, but inlcude side lobes for temporal dynamics
                fdatAt = filterFGx(datais,srate,usehz,peakwidt,0);                
                
                % take absolute of hilbert transform, to show dynamic SSEP
                dynSSEP= abs(hilbert(fdatAt'))';
                
                
%                 %TAKE mean per trial, careful with baseline subtraction
%                 timeid = [0:1/srate:epochdur];
%                 timeid= timeid-3;
%                 if dirsi==1 %increasing PFI / Target absence
%                     if ifreq<5 %targets, so low point post event.
%                         windcheck= [2 3];
%                     else %BG hz, so
%                         windcheck= [-3 -1];
%                     end
%                 elseif dirsi==0 %oppoite case, decreasing PFI/ targets are appearing
%                     if ifreq<5 %targets, so low point pre event.
%                         windcheck= [-3 -2];
%                     else %BG hz, so
%                         windcheck= [2 3];
%                     end
%                 end

%                 
                else % use SNR
                     
                            
                            param_spcgrm.tapers = [1 1];
                            param_spcgrm.Fs= [250];
                            param_spcgrm.Fpass= [0 50];
                            param_spcgrm.trialave=0;
                            movingwin=[1,.15];
                            
                            
                            [sgrm ,tgrm, fgrm] = mtspecgramc(datais', movingwin, param_spcgrm);
                            %%
                            %conv SNR
                            snr_sgrm =zeros(size(sgrm));
                            
                            
                            
                            kernelw= [-1/6 -1/6 -1/6 0 0 1 0 0 -1/6 -1/6 -1/6 ];
                            
                            %snr on trials or
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
%                             
                            %store just stim freq.
                            [~, hzid]= min(abs(fgrm-usehz));
                            %
% %                             %reduce size.
                            dynSSEP=squeeze(snr_sgrm(:,hzid,:))';
% %                             
                          timeidDYN=tgrm-3;

                            
                    
                    
                end
%                         
%                %remove base
if rmvbase==1
    tidx=dsearchn(timeidDYN', [bsrem]');
    for itrial=1:size(dynSSEP,1)
        tmp= dynSSEP(itrial,:);
        mtmp = mean(tmp(tidx(1):tidx(2)));
        %                    mtmp = mean(tmp);
        
        dynSSEP(itrial,:)= tmp - repmat(mtmp, 1, length(tmp));
        
    end
    
end
                
               %
               
                %%\
                snrgrm20=dynSSEP;
                % smooth across trials
                sm_snrgrm20=zeros(size(snrgrm20));
                for itime=1:size(snrgrm20,2)
                    sm_snrgrm20(:,itime)= smooth(snrgrm20(:,itime),5);
                end
                %%
                hold on
                subplot(3,2,[2 4]);
                hold on
                imagesc(timeidDYN,  1:size(dynSSEP,1), sm_snrgrm20);
                c=colorbar;
                if useSNRorhilbert==1
                ylabel(c, 'log(SNR)')
                else
                ylabel(c, 'abs(hilbert)')
                end
                
                xlabel('Time secs')
%                 caxis([-1*max(max(sm_snrgrm20)) max(max(sm_snrgrm20))])                %                 set(gca, 'ytick', 1:24, 'yticklabel', round(sortedRTs./60,2))
caxis([ -1 3])
                axis tight
                hold on
                plot([0 0] , ylim, ['w:'], 'linewidth', 3)
                title({['sorted ' num2str(usehz) 'Hz RESS data (smoothed) for  ' ctype ','];['ppant' num2str(ifol)]})
                set(gca, 'fontsize', 15, 'yticklabel', [])
                
                %         ylim([0 20])
                %% show mean over time.
                subplot(3,2,6)
                plot(timeidDYN, nanmean(snrgrm20,1),['k-'], 'linewidth', 3)
                %         shadedErrorBar([-3:1/60:3], nanmean(acrBP,1),nanstd(acrBP,1))
                set(gca, 'fontsize', 15)
                hold on;
                ylabel('RESS logSNR')
                %%
                axis tight
                plot([0 0] , ylim, ['k:'], 'linewidth', 1)
                %
                xlabel('Time secs')
                set(gcf,'color', 'w');
                
                %         ylim([0 20])
                
                %%
                
                print('-dpng', ['figure_PFI_' ctype '_summary_BPandSSVEP_' num2str(usehz) '_RESS.png']);
                
                
                switch itimezero
                    case 1 %store for across ppant plots:
                        Ppant_onsetBP=BPstouse;
                        Ppant_onsetSNR=snrgrm20; %sorted.
%                         Ppant_onsetRTs=RTstouse;
%                         Ppant_onsetTOPO=snr20;
                    case 2
                        
                        Ppant_offsetBP=BPstouse;
                        Ppant_offsetSNR=snrgrm20; %always sorted in descending order of PFI.
%                         Ppant_offsetRTs=RTstouse;
%                         Ppant_offsetTOPO=snr20;
                end
            end
            
            
            %%
            switch ihz
                case 5
                    savename='PFIperformance_withSNR_20_RESS';
                case 6
                    savename='PFIperformance_withSNR_40_RESS';
            end
            
            
            save(savename,...
                'Ppant_onsetBP','Ppant_offsetBP',...
                'Ppant_onsetSNR', 'Ppant_offsetSNR', 'timeidDYN')
            
            
        end
    end
end




if job.concaterpdataacrossppants_keepnumberSeparate==1
    % as the previous, except now we will store the outgoing MEAN 
%SSVEP SNR for accumperiods in thirds 0-1 1-2 2and greater  
    
    
    dursMINIMUM= 0; %1 second. % or change to  %this works.
%     dursMINIMUM=123; % selects top third!.
    
    ppantsmoothing=1; % average across participants, after smoothing, or no.
    
    for ihz=1:2
        
        if ihz==1
            loadname='PFIperformance_withSNR_20_RESS';
        else
            loadname='PFIperformance_withSNR_40_RESS';
        end
        storeacrossPpant_onsetSNR_PFI1=[];%zeros(length(allppants),1501);
        storeacrossPpant_onsetSNR_PFI2=[];%zeros(length(allppants),1501);
        storeacrossPpant_onsetSNR_PFI3=[];%zeros(length(allppants),1501);
        
        storeacrossPpant_offsetSNR_PFI1=[];%zeros(length(allppants),1501);
        storeacrossPpant_offsetSNR_PFI2=[];%zeros(length(allppants),1501);
        storeacrossPpant_offsetSNR_PFI3=[];%zeros(length(allppants),1501);
        
        storeacrossPpant_offsetSNR_PFI_trialbytrial= zeros(length(allppants), 100,33);
        storeacrossPpant_onsetSNR_PFI_trialbytrial= zeros(length(allppants), 100,33);
        icounter=1;
        
        for ippant = allppants
             cd(basefol)
                cd(num2str(ippant));
                %onset types
                
                load(loadname)
%                 Ppant_onsetSNR=Ppant_onsetSNR';
%                 Ppant_offsetSNR=Ppant_offsetSNR';
                
                %restrict to only 'long' PFIs if needed.
                
                    if dursMINIMUM>0 && dursMINIMUM~=123
                        
                        %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        shrton= find(Ppant_onsetRTs<dursMINIMUM); %note switch
                        all_onset(shrton)=[];
                        keepON=all_onset;
                        
                        %same for offsets.
                        shrtoff=find(Ppant_offsetRTs<dursMINIMUM);
                        all_offset(shrtoff)=[];
                        keepOFF=all_offset;
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
                    elseif dursMINIMUM==123
                        
                        %select top third of disap durations.
                         %restrict to those trials of interest!
                        all_onset = 1:size(Ppant_onsetBP,1);
                        all_offset = 1:size(Ppant_offsetBP,1);
                        
                        
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_onsetRTs, 'descend');
                        %top third
                        keepON = nid(1:ceil(length(ndurs)/4));
                        
                        %same for offsets.
                        %find top third of distribution.
                        [ndurs,nid] = sort(Ppant_offsetRTs, 'descend');
                        %top third
                        keepOFF = nid(1:ceil(length(ndurs)/3));
                        
                        Ppant_onsetBP=Ppant_onsetBP(keepON,:);
                        Ppant_offsetBP=Ppant_offsetBP(keepOFF,:);
                        
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
                    
                    % now sort into chunks - median split? 
                    
                    for id=1:2 %onsets and offsets
                        switch id
                            case 1
                                dtype=Ppant_onsetBP;
%                                 accumPeriod = (sum(dtype(:, 180:end),2)/60)/3;

                            case 2
                                dtype=Ppant_offsetBP;
%                                 accumPeriod = (sum(dtype(:, 1:180),2)/60)/3;
                        end
                  
                        
                        switch id
                            case 1
                                p1D_on= mean(Ppant_onsetSNR(51:100,:),1);
                                p2D_on= mean(Ppant_onsetSNR(1:50,:),1);
%                              
                            case 2
                                p1D_off= mean(Ppant_offsetSNR(51:100,:),1);
                                p2D_off= mean(Ppant_offsetSNR(1:50,:),1);
%                                
                                
                        end
                                        
                    end
                    
                    %% sanity check:
%                     figure(1); clf
%                     plot(1:size(p1D_off,2), p1D_off), hold on;
%                     plot(1:size(p1D_off,2), p2D_off, 'k');
%                     title(num2str(ippant))
%                     shg
%                     figure(2);imagesc(Ppant_offsetSNR); colorbar, caxis([-5 5]); hold on, plot([xlim], [50 50], ['-k'])
%                   
%%
%                     storeacrossPpant_onsetSNR_PFI1(icounter,:)= p1D_on;
%                     storeacrossPpant_onsetSNR_PFI2(icounter,:)= p2D_on;
%                     
%                     storeacrossPpant_offsetSNR_PFI1(icounter,:)= p1D_off;
%                     storeacrossPpant_offsetSNR_PFI2(icounter,:)= p2D_off;

%or save trialbytrial image
storeacrossPpant_offsetSNR_PFI_trialbytrial(icounter,:,:)=Ppant_offsetSNR;
storeacrossPpant_onsetSNR_PFI_trialbytrial(icounter,:,:)=Ppant_onsetSNR;
                    icounter=icounter+1;
        end      
        %save appropriately
        cd(basefol)
        cd('newplots-MD')
        
        
        switch ihz
            case 1
                savename=['GFX_PFIperformance_withSNR_20_min' num2str(dursMINIMUM) '_RESS'];
            case 2
                savename=['GFX_PFIperformance_withSNR_40_min' num2str(dursMINIMUM) '_RESS'];
        end
        
        
        save(savename,...
'storeacrossPpant_onsetSNR_PFI_trialbytrial','storeacrossPpant_offsetSNR_PFI_trialbytrial','-append')
    %             'storeacrossPpant_onsetSNR_PFI1',...
%                     'storeacrossPpant_onsetSNR_PFI2',...
%                     'storeacrossPpant_onsetSNR_PFI3',...                    
%                     'storeacrossPpant_offsetSNR_PFI1',...
%                     'storeacrossPpant_offsetSNR_PFI2',...
%                     'storeacrossPpant_offsetSNR_PFI3','timeid','-append')
        
        
        
    end
    
    
    
end
if job.BPandSSVEPtimecourseacrossppants_numSeparate==1
    %%
  %%
    getelocs
    rmvbase=0;
    clf
    checkcluster=1;
    clf
    plcount=1;
    %
    superimposeBP=0; % for adding BP results (using median split also).
    
    for hzis=1:2
        plcount=1;
        cd(basefol)
        cd('newplots-MD')
        switch hzis
            case 1
                loadname='GFX_PFIperformance_withSNR_20';
%                        col=['r'];
col='b';
                       
                lint='-';
                usehz=20;
                sigheight=1.45;
            case 2
                
                loadname='GFX_PFIperformance_withSNR_40';
%                        col=[0 .5 0];
                       col='m';
                usehz=40;
                lint='-';
                sigheight=1.5;
        end
        
        %%
        anovatestdata=[];
        for itimezero=1:2
        
        
        for idur=1%:2
             cd(basefol)
             cd('newplots-MD')
            switch  idur
                case 1
                    load([loadname '_min0_RESS'])
                    
                case 2
                    load([loadname '_min30_RESS'])
                case 3
                    load([loadname '_min60_RESS'])
            end
            
            tbase = timeidDYN;
            useSNRNUM=[];
            switch itimezero
                case 1
                    
%                     useSNRNUM(1,:,:)=squeeze(storeacrossPpant_onsetSNR_PFI1);
%                     useSNRNUM(2,:,:)=squeeze(storeacrossPpant_onsetSNR_PFI2);
%                   
%                     useSNRNUM(1,:,:)=squeeze(nanmean(storeacrossPpant_onsetSNR_PFI_trialbytrial(:,1:25,:),2));
%                     useSNRNUM(2,:,:)=squeeze(nanmean(storeacrossPpant_onsetSNR_PFI_trialbytrial(:,26:50,:),2));
%                     useSNRNUM(3,:,:)=squeeze(nanmean(storeacrossPpant_onsetSNR_PFI_trialbytrial(:,51:75,:),2));
%                     useSNRNUM(4,:,:)=squeeze(nanmean(storeacrossPpant_onsetSNR_PFI_trialbytrial(:,76:100,:),2));
                    
                      
                    useSNRNUM(1,:,:)=squeeze(nanmean(storeacrossPpant_onsetSNR_PFI_trialbytrial(:,1:50,:),2));
                    useSNRNUM(2,:,:)=squeeze(nanmean(storeacrossPpant_onsetSNR_PFI_trialbytrial(:,51:100,:),2));
%                     useSNRNUM(3,:,:)=squeeze(nanmean(storeacrossPpant_onsetSNR_PFI_trialbytrial(:,67:100,:),2));

                    
                    chtype='target invisible';
%                     chtype='button press';
                    
             useBPNUM(1,:,:)= squeeze(nanmean(storeacrossPpant_onsetBP(:,51:100,:),2));
                    useBPNUM(2,:,:)= squeeze(nanmean(storeacrossPpant_onsetBP(:,1:50,:),2));
                    
                case 2
                    
                   
%                     useSNRNUM(1,:,:)=squeeze(storeacrossPpant_offsetSNR_PFI1);
%                     useSNRNUM(2,:,:)=squeeze(storeacrossPpant_offsetSNR_PFI2);

%                     useSNRNUM(1,:,:)=squeeze(nanmean(storeacrossPpant_offsetSNR_PFI_trialbytrial(:,1:25,:),2));
%                     useSNRNUM(2,:,:)=squeeze(nanmean(storeacrossPpant_offsetSNR_PFI_trialbytrial(:,26:50,:),2));
%                     useSNRNUM(3,:,:)=squeeze(nanmean(storeacrossPpant_offsetSNR_PFI_trialbytrial(:,51:75,:),2));
%                     useSNRNUM(4,:,:)=squeeze(nanmean(storeacrossPpant_offsetSNR_PFI_trialbytrial(:,76:100,:),2));


                    useSNRNUM(1,:,:)=squeeze(nanmean(storeacrossPpant_offsetSNR_PFI_trialbytrial(:,1:50,:),2));
                    useSNRNUM(2,:,:)=squeeze(nanmean(storeacrossPpant_offsetSNR_PFI_trialbytrial(:,51:100,:),2));
%                     useSNRNUM(3,:,:)=squeeze(nanmean(storeacrossPpant_offsetSNR_PFI_trialbytrial(:,67:100,:),2));


chtype='target visible';
%                     chtype='button release';




                    useBPNUM(1,:,:)= squeeze(nanmean(storeacrossPpant_offsetBP(:,51:100,:),2));
                    useBPNUM(2,:,:)= squeeze(nanmean(storeacrossPpant_offsetBP(:,1:50,:),2));
                    
                    
            end
            
            
            
%             
            %% SNR
            figure(1); 
            % 
           
            
            for inum=1:size(useSNRNUM,1)
           if inum==2
               lint='-';
           else 
%                lint=':';
           end
                if superimposeBP==1
                    
                    subplot(2,2,itimezero)
                    hold on
                    %need to resample to fit on same plot.
                    BPpl= squeeze(useBPNUM(inum,:,:));
                    %                     BPpl=resample(BPis', length(tbase),length(BPis))';
                    xt = [1:length(BPpl)]/60-3;
                    stE = std(BPpl)/sqrt(size(BPpl,1));
                    sh=shadedErrorBar(xt, squeeze(mean(BPpl,1)) ,stE,[lint 'r'],[1]);
                    
                    sh.mainLine.LineWidth = inum*2;
                    xlabel(['Time from ' chtype])
                    set(gca, 'fontsize', 25)
                    ylabel('Buttons pressed')
                    
                    xlim([-2.5 2.5])
                end
                
                
                ppantMeanSNR=squeeze(useSNRNUM(inum,:,:));    
            
            plis = 3+(itimezero-1);
            subplot(2,2,plis)
            hold on
            %plot across ppant trace
            
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            x = ppantMeanSNR;
            
            mXppant =squeeze( nanmean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = nanmean(nanmean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            stE = nanstd(NEWdata)/sqrt(size(x,1));
            %
            sh=shadedErrorBar(tbase, nanmean(ppantMeanSNR,1),stE,[lint],[1]);
                
                sh.mainLine.Color = col;
                sh.mainLine.LineWidth = 2;
                sh.patch.FaceColor = col;
                sh.edge(1).Color = col;
                sh.edge(2).Color = col;
                
                
                
%                 if superimposeBP==1 % add BP 
%                     hold on
%                 %need to resample to fit on same plot.
%                 BPis= squeeze(useBPNUM(inum,:,:));
%                 BPpl=resample(BPis', length(tbase),length(BPis))';
%                 stE = std(BPpl)/sqrt(size(BPpl,1));
%                 shadedErrorBar(tbase, mean(BPpl,1) ,stE,[lint 'b'],[1]);
%                 
%                 end
            
                
                
                
                
                
                set(gca, 'fontsize', 15)
            hold on;
            %                 %%
            %
            ylabel({['RESS log(SNR)']})
%             ylim([1.4 2.6])
%             title({[num2str(usehz) ' Hz SSVEP']})
            xlabel(['Time from ' chtype])
            set(gca, 'fontsize', 25)
           
            set(gcf, 'color', 'w')
            anovatestdata([plcount],1:size(ppantMeanSNR,1),:) = ppantMeanSNR;
            legendprint(inum)=sh.mainLine;
            plcount=plcount+1;
            end
        end
        
        %%
        %check for sig 
        %use RMANOVA.
        %also calculate and plot significance, using ANOVA across
                %conditions.
                factorarray = [ones(21,1);repmat(2,21,1);repmat(3,21,1)];
                ss=[1:21]';
                subjects=[ss;ss;ss];
                
                %store ANOVA results
                pvals=nan(1,size(anovatestdata,3));
                Fvals=pvals;
                %store planned comparisons; for AttL vs nAttL
%                 presultPC=nan(1,size(anovatestdata,3));
                
                for itime=2:size(anovatestdata,3) %at each timeponit
                    
                    %arrange for rmanova. single cols.
%                     rmanova(data,factorarray,subjects,varnames, btw_ss_col)
%                     datat= [squeeze(anovatestdata(1,:,itime))'; squeeze(anovatestdata(2,:,itime))';...
%                         squeeze(anovatestdata(3,:,itime))'];
%                     
%                     anvret=rmanova(datat,factorarray,subjects);
%         pvals(itime)=anvret.p(1);
%         Fvals(itime)=anvret.table{2,6};
        
        
        
                 [~, pvals(itime),~,stat]=ttest(squeeze(anovatestdata(plcount-2,:,itime)),squeeze(anovatestdata(plcount-1,:,itime)), .05);
        Fvals(itime)=stat.tstat;
                end
%         %%
       
        sigs=find(pvals<.05);
%         
%             %perform cluster based correction.
            if length(sigs)>1 && checkcluster==1 && itimezero==1
                % find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                %grab largest
%                 ignore bad points.
                
                
                  % find biggest cluster:
                %finds adjacent time points
                sigs = find(pvals<.05);
                
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
               
                %%
                for icl=1%:size(clusterSTandEND,1)

                    %start and end are now:
                    STC=sigs(clusterSTandEND(maxClust,1));
                    ENDC=sigs(clusterSTandEND(maxClust,2)+1);
                    checktimes =STC:ENDC;
                    observedCV = sum(abs(Fvals(checktimes)));
                    % now shuffle condition labels to see if this cluster is
                    % sig (compared to chance).
                    sumTestStatsShuff = zeros(1,2000);
                    
                    for irand = 1:2000
                        %testing the null that it isn't mismatched - matched at time 2
                        % which creates a diff. so select from either!
                        shD= zeros(size(anovatestdata,1),length(checktimes),size(anovatestdata,2));
                        
                        %change shuffle parameters based on test of
                        %interest. (ie between conditions, or temporal
                        %null).
                        
                        for ipartition = 1:size(anovatestdata,1)
                            for ippant = 1:size(anovatestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    if mod(randi(100),2)==0 %if random even number
                                        pdata = anovatestdata(plcount-2,randi(size(anovatestdata,2)), checktimes(itime)); %select both chans
                                    else
                                        pdata = anovatestdata(plcount-1,randi(size(anovatestdata,2)), checktimes(itime));
                                    end
                                    
%                                         pdata = anovatestdata(randi(3),randi(size(anovatestdata,2)), checktimes(itime)); %select both chans
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                       
                        %now compute difference between out hypothetical topoplots,
                        % and test for sig, checking the accumulated test statistic at our
                        % times of interest
                        tvalspertimepoint = zeros(1,length(checktimes));
                        
                        testdata = squeeze(shD(1,:,:)) - squeeze(shD(2,:,:));
                        pvalsNew=[];
                        FvalsNew=[];
                        for itest = 1:length(checktimes) %test each time point
                            
%                             
%                             datat= [squeeze(shD(1,:,itest))'; squeeze(shD(2,:,itest))';...
%                         squeeze(shD(3,:,itest))'];
%                     
%                     anvret=rmanova(datat,factorarray,subjects);
% %                     pvalsNew(itime)=anvret.p(1);
%                     FvalsNew(itime)=anvret.table{2,6};
                    
%                             
%                             tvalspertimepoint(1,itest) = anvret.table{2,6};
                            
                            
                             [~, p, ~,stat]= ttest(testdata(itest,:));
                            
                            tvalspertimepoint(1,itest) = stat.tstat;
                            
                        end
                        
                        sumTestStatsShuff(1,irand) = sum(abs(tvalspertimepoint));
                    end %repeat nshuff times
                    
                    
                    %is the observed greater than CV?
                    % plot histogram:
                    figure(2);
                    
                        clf
                    
%                     subplot(2,1, plcount-1)
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
                    
                    %observed pvalue in distribution=
                    [~, c2] = min(abs(H.Data-observedCV)); %closest value.
                    pvalis= 1-cdf(c2);
                    title(['sum tvals = ' num2str(observedCV), 'p=' num2str(pvalis)]);
                    
                    %              title('Spatial Cluster  Significant!')
                        for itime=checktimes
                            figure(1);
                            hold on
                            plot(timeidDYN(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', col)
                        end
                    end
            
%            
            end
            else
                try %just plotem
                for itime=sigs
                            figure(1);
                            hold on
%                             plot(timeidDYN(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', col)
                end
                catch
                end
                
            end
            
            
           %itimezero
%                 
            %%
            xlim([-2.5 2.5])
            end
    end
            %%
%             figure(1);
%             hold on;
%             
%             axis tight
% %             ylim([-4 4.5])
%             hold on
%             plot(xlim, [0 0], ['k:']);
%             plot([0 0] , ylim, ['k:'], 'linewidth', 1) 
%             %plot patch:
%             title([num2str(usehz) ' Hz SSVEP'])
%                xlim([-2.5 2.5])
%             shg
%             
%             hold on
%             
%             %
%             
%         
%         %
%         
%         lg=legend(legendprint, {'PFI<median', 'PFI>median', '3 Targets PFI'});
%         if itimezero==1
%             
%             set(lg, 'location', 'NorthWest')
%               
%         else
%             
%             set(lg, 'location', 'SouthWest')
%         end
%                     set(lg, 'location', 'SouthWest')
%         %
%          cd(basefol)
%             cd('newplots-MD')
%             
%             %
            cd('Figures')
            %%
        print('-dpng', ['PFI trace ' num2str(usehz) 'Hz SSVEP summary, during ' num2str(chtype) '_sepNumberPFI.png'])      
end


if job.BPandSSVEPtimecourseacrossppants_group==1


   getelocs
    
    rmvbase=0;
    checksigON=1;
    checkcluster=1;
    
    clf
    plcount=1;
    legendprint=[];
    cd(basefol)
    cd('newplots-MD')
    for hzis=1:2%
        
        switch hzis
            case 1
                load('GFX_PFIperformance_withSNR_20_min0_RESS')
                
                lint=':';
                usehz=20;
                sigheight=1.55;
                if rmvbase==1
%                     sigheight=-3.7;
                     sigheight=1.9;
                end
                col='b';
            case 2
                
                load('GFX_PFIperformance_withSNR_40_min0_RESS')
                
                usehz=40;
                lint=':';
                sigheight=1.6;
                if rmvbase==1
%                     sigheight=-4.1;
                        sigheight=1.5;
                end
                
                col=[0 .5 0];
                col='m';
        end
        tgrm=timeid;
        tbase = tgrm;
        %%
        
        ttestdata=[];
%         legendprint=[];
        for itimezero=1:2
            
            
            switch itimezero
                case 1
                    useBP=storeacrossPpant_onsetBP;
                    useSNR=storeacrossPpant_onsetSNR;
                    
                    
                    chtype='target invisible';
%                     chtype='button press';
                    
%                     col=[0 .5 0];
                    lint='-';
                case 2
                    
                    useBP=storeacrossPpant_offsetBP;
                    useSNR=storeacrossPpant_offsetSNR;
                    

                    
                    
                    chtype='target visible';
                    chtype='subjective report';
%                     chtype='button release';
                    
                    
%                     col='r';
                    bsrem=[-3 -1];
                    lint=':';
                    
            end
            
            
            %%
%           
            %SNR
            figure(1);
%           
            
            for ip=2%1:2
                switch ip
                    case 1
                        used=useBP;
                        basep=1;
                        col='b';
                        lint='-';
                        ylimsare = [0 2.5];
                    case 2
                        used=useSNR;
                        basep=3;
                        if hzis==1
                            
                            lgc=1;
                        else
                            lgc=1;
                        end
                        
                        ylimsare= [ 1.5 2.7];
                end
                
                placeme = basep+ (1*itimezero-1);
                
%                 subplot(2,2,3)
            
            
            hold on
            
%             
%             if rmvbase==1
%                 tmp = zeros(size(useSNR));
%                 tIND = dsearchn(tbase', [bsrem]');
%             for ippant=1:size(useSNR,1)
%                 for itrial = 1:size(useSNR,2)
%                     tbase = squeeze(nanmean(useSNR(ippant,itrial,tIND(1):tIND(2)),3));
%                     rbase = repmat(tbase, [1  size(useSNR,3)]);
%                     
%                     tmp(ippant,itrial,:) =squeeze(useSNR(ippant,itrial,:))' - rbase;
%                 end
%             end
%             useSNR=tmp;
%             end
            
            %plot across ppant trace
            ppantMeanSNR= squeeze(mean(used,2));
            
            %adjust standard error as per COusineau(2005)
            %confidence interval for within subj designs.
            % y = x - mXsub + mXGroup,
            if ip~=1
            x = ppantMeanSNR;
            
            mXppant =squeeze( mean(x,2)); %mean across conditions we are comparing (within ppant ie. time points).
            mXgroup = mean(mean(mXppant)); %mean overall (remove b/w sub differences
            
            %for each observation, subjtract the subj average, add
            %the group average.
            NEWdata = x - repmat(mXppant, 1, size(x,2)) + repmat(mXgroup, size(x,1),size(x,2));
            
            %             %compute new stErr %which version?
            
            
            else
                NEWdata=ppantMeanSNR;
            end
                stE = std(NEWdata)/sqrt(length(allppants));
                
            if ip==1
                sh=shadedErrorBar(-3:1/60:3, mean(ppantMeanSNR,1),stE,[lint],[1]);
                ylabel(['Buttons pressed'])
            else
            sh=shadedErrorBar(tgrm, mean(ppantMeanSNR,1),stE,[lint],[1]);
            ylabel({['RESS log(SNR)']})
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
            xlabel(['Time from ' chtype ' [s]'])
%             xlabel('Time from perceptual report')
            set(gca, 'fontsize', 25)
            
            timeidDYN=tgrm;
            xlim([timeidDYN(1) timeidDYN(end)])
            ylim( [ylimsare])
            
            
            
            
            
            set(gcf, 'color', 'w')
            
            if ip==2
            ttestdata(itimezero,:,:) = ppantMeanSNR;
            legendprint(lgc)=sh.mainLine;
            end
            
            plcount=plcount+1;
            
            
            %place legend
            if hzis==2 && ip==2
                if itimezero==1
%                 lg=legend([legendprint(1) legendprint(2)], {'20 Hz', '40 Hz'});
                end
%             set(lg, 'location', 'NorthEast')
        
            end
            
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
        sigs=find(pvals<.05);
%         
%             %perform cluster based correction.
            if length(sigs)>2 &&checkcluster==1
                % find biggest cluster:
                %finds adjacent time points
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                %grab largest
%                 ignore bad points.
                
                
                  % find biggest cluster:
                %finds adjacent time points
                sigs = find(pvals<.05);
                
                vect1 = diff(sigs);
                v1 = (vect1(:)==1);
                d = diff(v1);
                clusterSTandEND= [find([v1(1);d]==1) find([d;-v1(end)]==-1)];
                [~,maxClust] = max(clusterSTandEND(:,2)-clusterSTandEND(:,1));
               
                %%
                for icl=1:size(clusterSTandEND,1)

                    %start and end are now:
                    % change icl to maxClust if we only want the largest
                    % cluster.
                    STC=sigs(clusterSTandEND(icl,1));
                    ENDC=sigs(clusterSTandEND(icl,2)+1);
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
                            for ippant = 1:size(ttestdata,2)
                                for itime=1:length(checktimes)
                                    
                                    if mod(randi(100),2)==0 %if random even number
                                        pdata = ttestdata(1,randi(size(ttestdata,2)), checktimes(itime)); %select both chans
                                    else
                                        pdata = ttestdata(2,randi(size(ttestdata,2)), checktimes(itime));
                                    end
                                    
                                    shD(ipartition,itime,ippant) = pdata;
                                end
                            end
                        end
                        else %null is that there are no temporal coincident sig values.
                            for ipartition = 1:2
                            for ippant = 1:size(ttestdata,2)
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
                    
                        clf
                    
                    
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
                    timeidDYN=tgrm;
                        for itime=checktimes
                            figure(1);
                            hold on
                            plot(timeidDYN(itime), sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', sh.mainLine.Color)
%                             plot(tgrm(itime)-3, sigheight, ['*' ],'markersize', 15, 'linewidth', 3, 'color', 'm')
                        end
                    end
            
%            
            end
            end 
                
    end
            %%
            end

            hold on
%             plot(xlim, [0 0], ['k:']);
            plot([0 0] , ylim, ['k-'], 'linewidth', 1) 
            %plot patch:
            set(gca, 'box', 'on')
               
            shg
            
            hold on
            
            %%
           
          
        %%
        cd('Figures')
        %%
        print('-dpng', ['PFI trace Bground SSVEP summary, during ' num2str(chtype) '.png'])       
        
        %%
        
% print('-dpng', ['PFI trace Bground SSVEP summary, during both, 40Hz.png'])      
end