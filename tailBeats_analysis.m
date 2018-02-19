%%

clearvars; clc;

[FileName, PathName] = uigetfile('*.csv');
cd(PathName);

tailBeats = csvread(FileName);
tailBeats_corrected=[tailBeats(:,1),tailBeats(:,2)-90];

%% angle correction and filtering

tailAngles_filtered = tailBeats_corrected(:,2);
tailAngles_filtered(isnan(tailAngles_filtered))=0; % remoivng nans for filtering orientation
windowWidth = 31;
polynomialOrder = 3;
tailAngles_filtered = sgolayfilt(tailAngles_filtered, polynomialOrder, windowWidth);

%% frequency domain

y = tailAngles_filtered;
Fs = 300; % Sampling frequency
T = 1/Fs; % Sample time
L = length(tailAngles_filtered); % Length of signal
t = (0:L-1)*T; % Time vector

NFFT = 2^nextpow2(L); % Next power of 2 from signal length
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% %% some plotting
% 
% % plot single-sided amplitude spectrum.
% fig1 = figure; plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
% 
% % compare raw and filtered signal
% fig2 = figure;
% plot(tmp_delta_ori_filtered);
% hold on;
% plot(tailBeats_corrected(:,2));
% hold off;
% title('Compare filtered signal to the raw one - another sanity check')
% xlabel('Frame')
% ylabel('Tail Angle (deg)')
% 




% plot final filetered tail angles
xx= 1:3.3:96832;
fig4 = figure;
plot(xx,tailAngles_filtered);
title('Tail Angle Trace over Time')
xlabel('Time(ms)')
ylabel('Tail Angle (deg)')




% dcm1 = datacursormode(fig1);
% set(dcm1, 'UpdateFcn', @Data_Cursor_precision, 'Enable', 'on');

%% %%%%%%%%%%%%%%%%%%%% 

%plot a gif
xx = 1:3.3:96832;
xx = xx';
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = CatStr(PathName,'testAnimated.gif');
for n = 1:150:29343
    % Draw plot
    y = tailAngles_filtered(1:n);
    x = xx(1:n);
    plot(x,y,'linewidth', 0.25, 'color', 'red');
    ylim([-20 20]);
    xlim([0 96832]);
    xlabel('Time (ms)');
    ylabel('Tail Angle (deg)');
    grid on;
    set(gcf,'color','w'); % set figure background to white
    drawnow; 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
 %%%%%%%%%%%%%%%%%%%%