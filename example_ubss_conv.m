function example_ubss_conv()

%Description:
%
% example_ubss_conv: example of underdetermined blind source separation 
% algorithms, including duet, duet_lsf, multi_nmf, fullrank. The mixing 
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
data_dir = './data/sisec2011 dev1/2mix4src130ms/';
[x1,fs]=audioread([data_dir 'mix_ch1.wav']);
[x2,fs]=audioread([data_dir 'mix_ch2.wav']);
nsampl=length(x1);
x=zeros(2,nsampl);
x(1,:)=x1;
x(2,:)=x2;
nsrc=4;
s=zeros(nsrc,nsampl);
for n=1:nsrc
    tmp=audioread([data_dir 'src' num2str(n) '.wav']);
    s(n,:)=tmp;
end
%% DUET
% addpath('./duet')
% addpath('./duet/utils')
% stft_win_len = 1024;
% stft_win_overlap=0.75*stft_win_len;
% isPlot_hist=0;
% % performing separation
% fprintf('Performing DUET...\n');
% tic
% se_duet=duet(x,nsrc,fs,stft_win_len,stft_win_overlap,isPlot_hist);
% time_duet=toc;
% fprintf('Separation done in %.2f seconds.\n',time_duet);
% % evaluation
% fprintf('Evaluation of estimated source signals...\n');
% [sdr_duet,sir_duet,sar_duet,~]=bss_eval_sources(se_duet,s);
% fprintf('SDRs are %.2f dB, %.2f dB, %.2f dB, and %.2f dB, respectively. Average is %.2f dB\n',...
%     sdr_duet(1),sdr_duet(2),sdr_duet(3),sdr_duet(4),mean(sdr_duet));
% % writeout estimated source signals
% result_dir = './sep/duet/';
% if ~exist(result_dir,'dir')
%     mkdir(result_dir);
% end
% fprintf(['Writing out separated signals in ' result_dir ' ...\n']);
% for n=1:nsrc
%     tmp=se_duet(n,:);
%     audiowrite([result_dir 'se' num2str(n) '.wav'],tmp,fs);
% end
% fprintf('Done.\n');
%% DUET_ISR
addpath('duet_isr');
addpath('./duet_isr/utils');
stft_win_len = 1024;
stft_win_overlap=0.75*stft_win_len;
mu=0.05;
lsf_na='isr';
% performing separation
fprintf('Performing DUET_LSF...\n');
tic
se_duetlsf = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu,lsf_na);
time_duetisr=toc;
fprintf('Separation done in %.2f seconds.\n',time_duetisr);
% evaluation
fprintf('Evaluation of estimated source signals...\n');
[sdr_duetlsf,sir_duetlsf,sar_duetlsf,~]=bss_eval_sources(se_duetlsf,s);
fprintf('SDRs are %.2f dB, %.2f dB, %.2f dB, and %.2f dB, respectively. Average is %.2f dB\n',...
    sdr_duetlsf(1),sdr_duetlsf(2),sdr_duetlsf(3),sdr_duetlsf(4),mean(sdr_duetlsf));
% writeout estimated source signals
result_dir = './sep/duet_lsf/';
if ~exist(result_dir,'dir')
    mkdir(result_dir);
end
fprintf(['Writing out separated signals in ' result_dir ' ...\n']);
for n=1:nsrc
    tmp=se_duetlsf(n,:);
    audiowrite([result_dir 'se' num2str(n) '.wav'],tmp,fs);
end
fprintf('Done.\n');
%% Multi_NMF
% addpath('multi_nmf');
% addpath('./multi_nmf/utils');
% addpath('./multi_nmf/ltfat-2.4.0-src/ltfat')
% ltfatstart
% stft_win_len = 2048;
% stft_win_overlap=0.75*stft_win_len;
% cpsn=10;
% % performing separation
% fprintf('Performing MULTI_NMF...\n');
% tic
% [se_nmf,ie_nmf]=bss_multinmf(x,nsrc,fs,stft_win_len,stft_win_overlap,cpsn);
% time_nmf=toc;
% fprintf('Separation done in %.2f seconds.\n',time_nmf);
% % evaluation
% ie_nmf=squeeze(ie_nmf(:,:,1)); % only consider the source image of the first channel
% fprintf('Evaluation of estimated source signals...\n');
% [sdr_nmf,sir_nmf,sar_nmf,~]=bss_eval_sources(ie_nmf,s);
% fprintf('SDRs are %.2f dB, %.2f dB, %.2f dB, and %.2f dB, respectively. Average is %.2f dB\n',...
%     sdr_nmf(1),sdr_nmf(2),sdr_nmf(3),sdr_nmf(4),mean(sdr_nmf));
% % writeout estimated source signals
% result_dir = './sep/multi_nmf/';
% if ~exist(result_dir,'dir')
%     mkdir(result_dir);
% end
% fprintf(['Writing out separated signals in ' result_dir ' ...\n']);
% for n=1:nsrc
%     tmp=ie_nmf(n,:);
%     audiowrite([result_dir 'se' num2str(n) '.wav'],tmp,fs);
% end
% fprintf('Done.\n');
%% full_rank
% addpath('./fullrank');
% addpath('./fullrank/utils')
% addpath('./fullrank/ltfat-2.4.0-src/ltfat')
% ltfatstart
% stft_win_len = 2048; % STFT window length
% stft_win_overlap=0.75*stft_win_len;
% K=30;        % No. of clusters for hierarchical clustering algorithm
% d=0.05;      % microphone spacing for test mixture
% EMIter=20;   % No. of EM iterations
% tic
% [ie_fullrank,~]=bss_fullrank(x,nsrc,d,fs,EMIter,stft_win_len,stft_win_overlap,K);
% time_fullrank=toc;
% fprintf('Separation done in %.2f seconds.\n',time_fullrank);
% % evaluation
% ie_fullrank=squeeze(ie_fullrank(:,:,1));% only consider source image of the first channel
% fprintf('Evaluation of estimated source signals...\n');
% [sdr_fullrank,sir_fullrank,sar_fullrank,~]=bss_eval_sources(ie_fullrank,s);
% fprintf('SDRs are %.2f dB, %.2f dB, %.2f dB, and %.2f dB, respectively. Average is %.2f dB\n',...
%     sdr_fullrank(1),sdr_fullrank(2),sdr_fullrank(3),sdr_fullrank(4),mean(sdr_fullrank));
% % writeout estimated source signals
% result_dir = './sep/fullrank/';
% if ~exist(result_dir,'dir')
%     mkdir(result_dir);
% end
% fprintf(['Writing out separated signals in ' result_dir ' ...\n']);
% for n=1:nsrc
%     tmp=ie_fullrank(n,:);
%     audiowrite([result_dir 'se' num2str(n) '.wav'],tmp,fs);
% end
% fprintf('Done.\n');