clear;
clc;

%Load .wav files
[intro1, ~] = audioread("data/Intro/c7Mic3Intro.wav");
[intro2, ~] = audioread("data/Intro/c8Mic3Intro.wav");

% %Spectrogram with default settings
% figure('Name','Defualt Spectrogram of Guy Introduction Sample')
% subplot(1,2,1);
% spectrogram(intro1,'yaxis');
% title("Mic7 (Left Ear) Intro Spectrogram");
% subplot(1,2,2);
% spectrogram(intro2,'yaxis');
% title("Mic8 (Right Ear) Intro Spectrogram");

%Spectrogram with custom settings
nsc = 256;
nov = 0;
nff = 256;
figure('Name','Custom Spectrogram of Guy Introduction Sample')
subplot(1,2,1);
spectrogram(intro1,ones(nsc,1),nov,nff,'yaxis');
title("Mic7 (Left Ear) Intro Spectrogram")
subplot(1,2,2);
spectrogram(intro2,ones(nsc,1),nov,nff,'yaxis');
title("Mic8 (Right Ear) Intro Spectrogram");
clear nff nov nsc;

%Design filterbank
[y1,y2,y3,y4] = AnalysisFilterBank(intro1);
[z1,z2,z3,z4] = AnalysisFilterBank(intro2);

% %Plot spectrogram of each sub-channel
% figure('Name','Spectrogram of Sub-channels')
% subplot(2,2,1);
% spectrogram(y1,ones(256,1),0,256,'yaxis');
% title("Low (First) Frequency Channel")
% subplot(2,2,2);
% spectrogram(y2,ones(256,1),0,256,'yaxis');
% title("Second Frequency Channel")
% subplot(2,2,3);
% spectrogram(y3,ones(256,1),0,256,'yaxis');
% title("Third Frequency Channel")
% subplot(2,2,4);
% spectrogram(y4,ones(256,1),0,256,'yaxis');
% title("High (Fourth) Frequency Channel")

%Chunking
y1 = chunk(y1);
y2 = chunk(y2);
y3 = chunk(y3);
y4 = chunk(y4);
z1 = chunk(z1);
z2 = chunk(z2);
z3 = chunk(z3);
z4 = chunk(z4);
y = {y1, y2, y3, y4};
z = {z1, z2, z3, z4};

%Calculate lag index
lagindex = zeros(size(y1,2),4);
for j=1:4
    for i=1:size(y{j},2)
        [cor, lag] = xcorr(y{j}(:,i),z{j}(:,i)); 
        [~, ind] = max(abs(cor));
        lagindex(i,j) = lag(ind);
    end
end
%Plot
figure('Name','Histogram of Lag Indexes for Guy Introductoin Sample');
subplot(4,1,1);
h1 = histogram(lagindex(:,1));
ylabel("Low Frequency Channel");
title("Histogram of Lag Indexes for Guy Introduction Sample");
subplot(4,1,2);
h2 = histogram(lagindex(:,2));
subplot(4,1,3);
h3 = histogram(lagindex(:,3));
subplot(4,1,4);
h4 = histogram(lagindex(:,4));
ylabel("High Frequency Channel");
xlabel("Index of Maximum Correlation")
%Calculate mu and sigma
h = {h1;h2;h3;h4};
mu = zeros(4,1);
sigma2 = zeros(4,1);
for i=1:4
    c = h{i}.Values;
    d = h{i}.BinEdges(1:end-1) + h{i}.BinWidth/2;
    mu(i) = sum(c.*d)/sum(c);
    sigma2(i) = sum(((d - mu(i)*ones(size(d))).^2 ).*c)/sum(c);
end
sigma = sqrt(sigma2);

clear c d
clear h h1 h2 h3 h4
clear intro1 intro2


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%REPEAT STEPS 1 TO 4 FOR MAIN FILES
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load .wav files
[main1, Fs] = audioread("data/main/c7m19.wav");
[main2, ~] = audioread("data/main/c8m19.wav");

% %Spectrogram with default settings
% figure('Name','Defualt Spectrogram of Conversation Samples')
% subplot(1,2,1);
% spectrogram(main1,'yaxis');
% title("Mic7 (Left Ear) Conversation Signal Spectrogram");
% subplot(1,2,2);
% spectrogram(main2,'yaxis');
% title("Mic8 (Right Ear) Conversation Signal Spectrogram");

%Spectrogram with custom settings
nsc = 256;
nov = 0;
nff = 256;
figure('Name','Custom Spectrogram of Conversation Samples')
subplot(1,2,1);
spectrogram(main1,ones(nsc,1),nov,nff,'yaxis');
title("Mic7 (Left Ear) Conversation Signal Spectrogram")
subplot(1,2,2);
spectrogram(main2,ones(nsc,1),nov,nff,'yaxis');
title("Mic8 (Right Ear) Conversation Signal Spectrogram");
clear nff nov nsc;

%Design filterbank
[y1,y2,y3,y4] = AnalysisFilterBank(main1);
[z1,z2,z3,z4] = AnalysisFilterBank(main2);

%Chunking
y1 = chunk(y1);
y2 = chunk(y2);
y3 = chunk(y3);
y4 = chunk(y4);
z1 = chunk(z1);
z2 = chunk(z2);
z3 = chunk(z3);
z4 = chunk(z4);
y = {y1, y2, y3, y4};
z = {z1, z2, z3, z4};

%Calculate lag index
lagindex = zeros(size(y1,2),4);
for j=1:4
    for i=1:size(y{j},2)
        [cor, lag] = xcorr(y{j}(:,i),z{j}(:,i)); 
        [~, ind] = max(abs(cor));
        lagindex(i,j) = lag(ind);
    end
end
%Plot
figure('Name','Histogram of Lag Indexes for Conversation Samples');
subplot(4,1,1);
histogram(lagindex(:,1));
ylabel("Low Frequency Channel");
title("Histogram of Lag Indexes for Conversation Samples");
subplot(4,1,2);
histogram(lagindex(:,2));
subplot(4,1,3);
histogram(lagindex(:,3));
subplot(4,1,4);
histogram(lagindex(:,4));
ylabel("High Frequency Channel");
xlabel("Index of Maximum Correlation")

%Calculate cell weights
Mu = repmat(mu',size(lagindex,1),1);
Sigma2 = repmat(sigma2', size(lagindex,1),1);
w = exp( (-1)*((lagindex - Mu).^2)./(2*Sigma2) );
clear Mu Sigma2

clear y1 y2 y3 y4
clear z1 z2 z3 z4
clear cor lag ind lagindex

%Rebuild signals
yp = y;
zp = z;
for j=1:size(w,2)
    for i=1:size(w,1)
        yp{j}(:,i) = w(i,j) * y{j}(:,i);
        zp{j}(:,i) = w(i,j) * z{j}(:,i);
    end
end
N = length(main1);
prosig1 = SynthesisFilterBank(yp, N);
prosig2 = SynthesisFilterBank(zp, N);

%Save processed signals
audiowrite('results/processed_sig1.wav',prosig1,Fs);
audiowrite('results/processed_sig2.wav',prosig1,Fs);

%Filter weights
load('filters/f0.mat','f0');
fw = filter(f0, 1, w);

%Plot weights before and after filtering
n = 1:length(w);
figure('Name','Filtering of the Weight Matrix');
subplot(4,1,1);
plot(n,w(:,1),n,fw(:,1));
ylabel('Low Frequency Channel');
title("Filtering of the Weight Matrix");
legend("Original Weight Seq.","Filtered Weight Seq.");
subplot(4,1,2);
plot(n,w(:,2),n,fw(:,2));
subplot(4,1,3);
plot(n,w(:,3),n,fw(:,3));
subplot(4,1,4);
plot(n,w(:,4),n,fw(:,4));
ylabel('High Frequency Channel');

%Rebuild signals with filtered weights
for j=1:size(fw,2)
    for i=1:size(fw,1)
        yp{j}(:,i) = fw(i,j) * y{j}(:,i);
        zp{j}(:,i) = fw(i,j) * z{j}(:,i);
    end
end
fprosig1 = SynthesisFilterBank(yp, N);
fprosig2 = SynthesisFilterBank(zp, N);

%Save processed signals
audiowrite('results/fil_processed_sig1.wav',prosig1,Fs);
audiowrite('results/fil_processed_sig2.wav',prosig1,Fs);

clear y yp z zp N i j




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%FUNCTIONS
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y1, y2, y3, y4] = AnalysisFilterBank(x)
    %Build filters
    load('filters/f1.mat', 'f1');
    load('filters/f2.mat', 'f2');
    f3 = real( exp(-1i*(3/8)*pi*(0:(length(f2)-1))).*f2 );
    f4 = real( exp(-1i*(5/8)*pi*(0:(length(f2)-1))).*f2 );
    f5 = real( exp(-1i*(1)*pi*(0:(length(f1)-1))).*f1 );
    %Filter input signal
    y1 = filter(f1, 1, x);
    y2 = filter(f3, 1, x);
    y3 = filter(f4, 1, x);
    y4 = filter(f5, 1, x);
end

function prosig = SynthesisFilterBank(y, N)
    prosig = zeros(1,N);
    for j=1:4
        %De-chunk
        temp = y{j}(1:N);
        %Add
        prosig = prosig + temp;
    end
end

function prosig = SynthesisFilterBankwithFiltering(y, N)
    prosig = zeros(1,N);
    %Build filters
    load('filters/f1.mat', 'f1');
    load('filters/f2.mat', 'f2');
    f = {f1, f2, f2, f2};
    f{2} = real( exp(-1i*(3/8)*pi*(0:(length(f2)-1))).*f2 );
    f{3} = real( exp(-1i*(5/8)*pi*(0:(length(f2)-1))).*f2 );
    f{4}= real( exp(-1i*(1)*pi*(0:(length(f1)-1))).*f1 );
    for j=1:4
        %De-chunk
        temp = y{j}(1:N);
        %Filter
        temp = filter(f{j}, 1, temp); 
        %Add
        prosig = prosig + temp;
    end
end

function r = chunk(x)
    r = zeros(256, ceil(length(x)/256));
    r(1:length(x)) = x(:);
end
