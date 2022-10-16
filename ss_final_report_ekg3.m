% step 1
% load data
Fs = 250;%sampling rate
rawEKG = load('tEKG3.mat');
nosEKG = load('fEKG3.mat');
rawEKG = rawEKG.tEKG;
nosEKG = nosEKG.fEKG;
% Filtered signal
FilteredSignal = conv(Hd, nosEKG);
FilteredSignal = FilteredSignal(floor(length(Hd)/2)+[1:length(nosEKG)]);
FTFilteredSignal = fft(FilteredSignal);
% plot in freq domain
t = (0: length(rawEKG)-1)/Fs;
plot(t, rawEKG, 'b');
hold on;
plot(t, nosEKG,'r');
hold off;
axis([15 17 -150 500]);
xlabel('sec');
ylabel('amplitude');
legend('rawEKG', 'nosEKG');
% fast fourier transform
FTrawEKG = fft(rawEKG);
FTnosEKG = fft(nosEKG);
f = linspace(0, Fs, length(rawEKG));
figure
subplot(3, 1, 1)
plot(f, abs(FTrawEKG))
title('raw EKG in freq domain'); xlim([0, Fs]);
xlabel('freq(Hz)'); ylabel('magnitude')
subplot(3, 1, 2)
plot(f, abs(FTnosEKG))
title('noisy signal in freq domain'); xlim([0, Fs]);
xlabel('freq(Hz)'); ylabel('magnitude')
subplot(3, 1, 3)
plot(f, abs(FTFilteredSignal))
title('filtered signal in freq domain'); xlim([0, Fs]);
xlabel('freq(Hz)'); ylabel('magnitude')
% plot versus figure
figure
subplot(2, 1, 1)
plot(t, rawEKG, 'b');
hold on;
plot(t, nosEKG,'r');
hold off;
title('noise EKG vs raw EKG');
axis([15 17 -150 500]);
xlabel('sec'); 
ylabel('mV');
legend('rawEKG', 'nosEKG');
subplot(2, 1, 2)
plot(t, rawEKG, 'b');
hold on;
plot(t, FilteredSignal,'r');
hold off;
title('Filtered EKG vs raw EKG');
axis([15 17 -150 500]);
xlabel('sec');
ylabel('mV');
legend('rawEKG', 'nosEKG');

% step2: signal to noise ratio
snr1 = snr(rawEKG, nosEKG-rawEKG);
snr2 = snr(rawEKG, FilteredSignal-rawEKG);

% step3: baseline removal
baseline = msbackadj(t.', FilteredSignal.', 'WindowSize', 1, 'StepSize', 1, 'ShowPlot', 'yes');
FilteredSignal = FilteredSignal+baseline';% plus or minus
figure;
plot(t, FilteredSignal);
axis([15 17 -1000 1000]);
title('Baseline removed FilteredSignal')

% step 4: EKG features
%derivatived signal
h1=[ -1 -2 0 2 1]/8;
ConvSignal = conv(h1, FilteredSignal);
ConvSignal = ConvSignal(floor(length(h1)/2)+[1:length(FilteredSignal)]);%shift 2 delayed point y2=ConvSignal
figure;
plot((1: length(FilteredSignal))/Fs, FilteredSignal, 'k');
hold on;
plot((1: length(ConvSignal))/Fs, ConvSignal, 'b');
hold on;
% square signal
ConvSignal = ConvSignal/max(abs(ConvSignal));% y2=y2/max(abs(y2));
y3 = ConvSignal.^2;% y2 square
plot((1: length(ConvSignal))/Fs, y3*500,'r');% plot y3
hold on;
% % moving window
h2 = ones(1 ,31)/31;
y4 = conv(h2, y3);% conv
y4 = y4(floor(length(h2)/2)+[1:length(y3)]); % shift delay back
plot((1: length(y4))/Fs, y4*500, 'g');% plot y4
hold off;
axis([15 17 -300 1500]);
legend('filtered signal', 'derivative signal', 'y3*500', 'y4*500');
% setting threshold
max_h=max(y4);
thresh=mean(y4);
threshold(1:length(y4))=(thresh*3); % threshold 
MWI_set=(y4>(thresh*max_h)*3)';
figure;
plot([1:length(y4)]/Fs,threshold,'-k');
hold on;
plot([1:length(MWI_set)]/Fs,50*MWI_set,'-r');
hold on;
plot(t, FilteredSignal,'b');
hold off;
axis([15 17 -300 1500]);
legend('threshold', 'moving window', 'Filtered Signal');
MWI_set = MWI_set.';
left = find(diff([0 MWI_set])==1);
right = find(diff([MWI_set 0])==-1);
% find R
for i=1:length(left)
    [R_peak_value(i) R_location(i)] = max(FilteredSignal(left(i):right(i))); % record max value and its location betwwen left(i) and right(i)->thus R-peak
    R_location(i) = R_location(i)-1+left(i);  %  redefine R_location 
end
% find S
for i=1:length(left)
    [S_peak_value(i) S_location(i)] = min(FilteredSignal(left(i):right(i))); % record min value and its location betwwen left(i) and right(i)->thus S-peak
    S_location(i) = S_location(i)-1+left(i);  %  redefine R_location 
end
% find Q
for i=1:length(left)
    [Q_peak_value(i) Q_location(i)] = min(FilteredSignal(left(i):R_location(i))); % record min value and its location betwwen left(i) and right(i)->thus S-peak
    Q_location(i) = Q_location(i)-1+left(i);  %  redefine R_location 
end
figure;
plot(t,FilteredSignal , t(R_location) ,R_peak_value , 'ro', t(S_location), S_peak_value, 'g*', t(Q_location), Q_peak_value, 'r+');
legend('Filtered Signal', 'R points', 'S point', 'Q point');
axis([15 17 -250 1300]);
title('sample no.3 EKG');

HeartBeat = numel(R_location);
RR = diff(R_location)/Fs; %ms
AveRR = sum(RR)/(HeartBeat-1);






