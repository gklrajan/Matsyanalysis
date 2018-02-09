clearvars; clc;
load('Danionella5dpf_050218.mat');
for cc=1:size(freeSwim,2)
tmp5_freeSwim(cc).boutData=freeSwim(cc).boutData;
norm5_freeSwim(cc).boutData= StatisticalNormaliz(tmp5_freeSwim(cc).boutData,'standard');
pooledMedian5(cc,:) = nanmedian(norm5_freeSwim(cc).boutData);
pooledMAD5(cc,:) = mad(norm5_freeSwim(cc).boutData);
end
close();

pooled5dpf = norm5_freeSwim(1).boutData;
for dataSize=2:size(freeSwim,2)
pooled5dpf=[pooled5dpf;norm5_freeSwim(dataSize).boutData];
end
    

load('Danionella8dpf_040218.mat');
for cc=1:size(freeSwim,2)
tmp8_freeSwim(cc).boutData=freeSwim(cc).boutData;
norm8_freeSwim(cc).boutData= StatisticalNormaliz(tmp8_freeSwim(cc).boutData,'standard');
pooledMedian8(cc,:) = nanmedian(norm8_freeSwim(cc).boutData);
pooledMAD8(cc,:) = mad(norm8_freeSwim(cc).boutData);
end
close();

pooled8dpf = norm8_freeSwim(1).boutData;
for dataSize2=2:size(freeSwim,2)
pooled8dpf=[pooled8dpf;norm8_freeSwim(dataSize).boutData];
end

for ii=1:size(pooled8dpf,2)
[p(ii),h(ii),stats(ii)] = ranksum(pooled5dpf(:,ii),pooled8dpf(:,ii),'alpha',0.01);
end

% for ii=1:size(pooledRaw_MT,2)
% [h(ii),p(ii)] = ttest2(pooledRaw_MT(:,ii),pooledRaw_WT(:,ii));
% end

pooledRaw_Danio(:,12)=rad2deg(pooled5dpf(:,12));
pooledRaw_WT(:,12)=rad2deg(pooled8dpf(:,12));

jugaad = [3,6,7,8,9,12,13,16];
title_jugaad = ["bout duration in ms","turn angle in deg","mean velocity in mm/sec","peak velocity in mm/sec","total distance in mm","angular velocity in deg/sec","inter-bout interval in ms","number of tail beats per bout"];
xlabel_jugaad = ["duration (ms)","turn (deg)","mean velocity(mm/sec)","peak velocity (mm/sec)","total distance (mm)","angular velocity (deg/sec)","inter-bout interval (ms)","number of tail beats/bout"];
figure1 = figure;
for gg = 1:length(jugaad)
subplot(2,8,gg);
histogram(pooled5dpf(:,jugaad(gg)),70);        
title(title_jugaad(gg));
xlabel(xlabel_jugaad(gg));
ylabel({'# of bouts'});

subplot(2,8,gg+8);
histogram(pooled8dpf(:,jugaad(gg)),70);
title(title_jugaad(gg));
xlabel(xlabel_jugaad(gg));
ylabel({'# of bouts'});
end

% %to plot overlayed histograms of WT & mutant
% for kk = 1:length(jugaad)
% subplot(3,3,kk);
% histogram(pooledRaw_WT(:,jugaad(kk)),100,'FaceColor','g','FaceAlpha',1);
% hold on;
% histogram(pooledRaw_Danio(:,jugaad(kk)),100,'FaceColor','r','FaceAlpha',0.5)
% hold off;
% end