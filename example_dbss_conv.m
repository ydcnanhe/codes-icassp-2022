function example_dbss_conv()

%Description:
%
% example_ubss_conv: example of determined blind source separation 
% algorithms, including iva and ilrma. The mixing 
% data is live-recording data from SiSEC 2011 evaluation campaign 
% (http://sisec.wiki.irisa.fr/) convolutive mixtures of "Under-determined 
% speech and music mixtures" task
%
% Input:
% -----
%
% Output:
% ------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% (yhebh@connect.ust.hk)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
data_dir = './data/sisec2011 dev1/2mix2src130ms/';
[x1,fs]=audioread([data_dir 'mix_ch1.wav']);
[x2,fs]=audioread([data_dir 'mix_ch2.wav']);
nsampl=length(x1);
x=zeros(2,nsampl);
x(1,:)=x1;
x(2,:)=x2;
nsrc=2;
s=zeros(nsrc,nsampl);
for n=1:nsrc
    tmp=audioread([data_dir 'src' num2str(n) '.wav']);
    s(n,:)=tmp;
end
%% IVA
% addpath('./iva')
% epochs = 200;
% nfft = 1024;
% % performing separation
% fprintf('Performing IVA...\n');
% tic
% se_iva=ivabss(x, nfft, epochs);
% time_duet=toc;
% fprintf('Separation done in %.2f seconds.\n',time_duet);
% % evaluation
% fprintf('Evaluation of estimated source signals...\n');
% [sdr_iva,sir_iva,sar_iva,~]=bss_eval_sources(se_iva,s);
% fprintf('SDRs are %.2f dB and %.2f dB, respectively. Average is %.2f dB\n',...
%     sdr_iva(1),sdr_iva(2), mean(sdr_iva));
% % writeout estimated source signals
% result_dir = './sep/iva/';
% if ~exist(result_dir,'dir')
%     mkdir(result_dir);
% end
% fprintf(['Writing out separated signals in ' result_dir ' ...\n']);
% for n=1:nsrc
%     tmp=se_iva(n,:);
%     audiowrite([result_dir 'se' num2str(n) '.wav'],tmp,fs);
% end
% fprintf('Done.\n');
%% ILRMA
addpath('ilrma');
seed = 1; % pseudo random seed
refMic = 1; % reference microphone for back projection
ns = 2; % number of sources
fftSize = 4096; % window length in STFT [points]
shiftSize = 2048; % shift length in STFT [points]
nb = 2; % number of bases (for type=1, nb is # of bases for "each" source. for type=2, nb is # of bases for "all" sources)
it = 200; % number of iterations (define by checking convergence behavior with drawConv=true)
type = 1; % 1 or 2 (1: ILRMA w/o partitioning function, 2: ILRMA with partitioning function)
drawConv = false; % true or false (true: plot cost function values in each iteration and show convergence behavior, false: faster and do not plot cost function values)
normalize = true; % true or false (true: apply normalization in each iteration of ILRMA to improve numerical stability, but the monotonic decrease of the cost function may be lost. false: do not apply normalization)
% Fix random seed
RandStream.setGlobalStream(RandStream('mt19937ar','Seed',seed))
% % performing separation
% fprintf('Performing DUET_LSF...\n');
tic
[se_ilrma, ~] = bss_ILRMA(x',ns,nb,fftSize,shiftSize,it,type,refMic,drawConv,normalize);
time_ilrma=toc;
fprintf('Separation done in %.2f seconds.\n',time_ilrma);
% evaluation
fprintf('Evaluation of estimated source signals...\n');
se_ilrma=se_ilrma';
[sdr_ilrma,sir_ilrma,sar_ilrma,~]=bss_eval_sources(se_ilrma,s);
fprintf('SDRs are %.2f dB and %.2f dB, respectively. Average is %.2f dB\n',...
    sdr_ilrma(1),sdr_ilrma(2),mean(sdr_ilrma));
% writeout estimated source signals
result_dir = './sep/ilrma/';
if ~exist(result_dir,'dir')
    mkdir(result_dir);
end
fprintf(['Writing out separated signals in ' result_dir ' ...\n']);
for n=1:nsrc
    tmp=se_ilrma(n,:);
    audiowrite([result_dir 'se' num2str(n) '.wav'],tmp,fs);
end
fprintf('Done.\n');