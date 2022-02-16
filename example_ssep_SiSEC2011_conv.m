function example_ssep_SiSEC2011_conv()

%
% example_ssep_SiSEC2011_conv();
%
% DUET_LSF ISR algorithm for SiSEC 2011 evaluation campaign (http://sisec.wiki.irisa.fr/)
%   convolutive mixtures of "Under-determined speech and music mixtures" task
%
% Input:
% -----
%
% Output:
% ------
%
% estimated sources are written in the results_dir
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% If you use this code, please cite this paper
%
% Yudong He, He Wang, Qifeng Chen, and Richard H.Y. So
% "HARVESTING PARTIALLY-DISJOINT TIME-FREQUENCY INFORMATION FOR IMPROVING DEGENERATE UNMIXING ESTIMATION TECHNIQUE,"
% in ICASSP 2022 - 2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), in press. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
data_dir = './data/sisec2011 dev1/2mix3src130ms/';
results_dir = './';
addpath('duet_isr');
addpath('./duet_isr/utils');
% Load mixture signals
[mix,fs]=audioread([data_dir 'mixtures.wav']);
% Parameters
nsrc=3;
stft_win_len = 1024;
stft_win_overlap=0.75*stft_win_len;
mu=0.05;
resp_na='isr';
% Performing BSS
fprintf('Performing blind source separation...\n');
x=mix';
tic
se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu,resp_na);
fprintf('Separation done in %.2f seconds.\n',toc);
% Evaluation
nsampl=size(se,2);
s=zeros(nsrc,nsampl);
fprintf('Evaluation of estimated source signals...\n');
for n=1:nsrc
    tmp=audioread([data_dir 's' num2str(n) '.wav']);
    s(n,:)=tmp;
end
[sdr,sir,sar,perm]=bss_eval_sources(se,s);
fprintf('SDRs are %.2f dB, %.2f dB, and %.2f dB, respectively. Average is %.2f dB\n',sdr(1),sdr(2),sdr(3),mean(sdr));
% Writeout estimated source signals
fprintf('Writing out separated signals in result_dir...\n');
for n=1:nsrc
    tmp=s(n,:);
    audiowrite([results_dir 'se' num2str(n) '.wav'],tmp,fs);
end
fprintf('Done.\n');