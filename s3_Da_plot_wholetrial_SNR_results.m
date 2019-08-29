% s3_Da_plot_wholetrial_SNR_results
% clear all
cd('/Users/MattDavidson/Desktop/FromIrenes Mac/Real_Experiment/Data Experiment copy/newplots-MD');

load('Across all mSNR for mtspectrum, whole trials.mat')
load('Across all for mtspectrum, whole trials.mat')
load('MD_AllBP_perppant.mat', 'PFIppants')
goodppants=PFIppants;
% first plot avg mtspectrum across good ppants.
getelocs;
%
gooddataSNR = squeeze(allSNR_mtspectrum(goodppants,:,:));
gooddatamtspect=squeeze(acrossallmtspect(goodppants,:,:));

mAllSNR= squeeze(mean(gooddataSNR,1));
mAllmtspect= squeeze(mean(gooddatamtspect,1));

% indx for topoplots.
[~,id20] = min(abs(mtfreq - 20));
[~,id40] = min(abs(mtfreq - 40));
fontsize=25;
%plot now.
%%
job.plotwholetrialtopos=0;
job.plotwholetrialspectrumSNRvsmtspect=1;


cd('Figures')
cd('whole trial SNR, across all')
if job.plotwholetrialtopos==1
for ifreq=1:2
    switch ifreq
        case 1
            usef = id20;
            hz= '20';
        case 2
            usef=id40;
            hz = '40';
    end
    
figure(1);
%%
clf
topoplot(squeeze(mAllSNR(usef,:)), elocs(1:64) )
%%
c=colorbar;
caxis([0 20])

set(gca, 'fontsize', fontsize*2)
title([ num2str(hz) ' Hz'], 'fontsize', fontsize*3)
ylabel(c, 'SNR', 'fontsize', fontsize*3)

%%
printfilename = ['Topo for ' num2str(hz) 'Hz SNR, whole trial'];
print('-dpng', printfilename)


%% also pvalues
pplot=zeros(1, 64);
for ichan = 1:64
    tmp = squeeze(gooddataSNR(:,usef,ichan));
    
    [h,p]=ttest(tmp,0);
    
    pplot(1,ichan)=p;
end

thresh=fdr(pplot, .05)    ;
%%
clf
% if ifreq==1
% topoplot(log10(pplot), elocs(1:64), 'emarker2', {[5,29:31,38,39, 61:63], 'O' 'b' 15 15});
% else
    topoplot(log10(pplot), elocs(1:64));
% end
colormap('hot')
colormap(flipud(colormap))
%
c=colorbar;
caxis([min(log10(pplot)) log10(thresh)])
set(gca, 'fontsize', fontsize*2)
title({['Significant after FDR correction']}, 'fontsize', fontsize*2)
ylabel(c, 'log10(pvalues)', 'fontsize', fontsize*3)

%%
printfilename = ['Topo for ' num2str(hz) 'Hz SNR, whole trial pvalues'];
print('-dpng', printfilename)


end
end
%%
if job.plotwholetrialspectrumSNRvsmtspect==1
%now spectrum at representative channel (POz)
clf;
set(gcf, 'color', 'w')

% plot(mtfreq, squeeze(mAllmtspect(:,62)), 'color', 'k', 'linewidth', 5);
% hold on
plot(mtfreq, squeeze(mAllSNR(:,62)), 'color', 'k', 'linewidth', 5);
xlim([0 45])
% title({['Mean whole trial SNR of EEG-spectrum at channel POz'];['(\itN \rm\bf= ' num2str(length(goodppants)) ')']}, 'fontweight', 'bold', 'fontsize', fontsize)
xlabel('Hz')
ylabel('log_1_0(SNR)')

set(gca, 'fontsize', fontsize*2)
%%
print('-dpng', ['mean whole trial logSNR at POz across n=' num2str(length(goodppants))])

end