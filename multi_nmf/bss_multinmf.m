function [se,ie]=bss_multinmf(x,nsrc,fs,stft_win_len,stft_win_overlap,cpsn)
%
% Description:
%
% bss_multinmf: underdetermined blind source separation by multichannel
% non-negative factorization (optimized by EM algorithm). 
%
% Syntax:
% 
% [se,ie] = bss_multinmf(x,nsrc)
% [se,ie]=bss_multinmf(x,nsrc,fs)
% [se,ie]=bss_multinmf(x,nsrc,fs,stft_win_len)
% [se,ie]=bss_multinmf(x,nsrc,fs,stft_win_len,stft_win_overlap)
% [se,ie]=bss_multinmf(x,nsrc,fs,stft_win_len,stft_win_overlap,cpsn)
%
% Inputs:
%
% x - nchan x nsampl matrix containg two chanels of time domain mixture signals with
% nsampl samples.
%
% nsrc - number of source signals.
%
% fs - sampling frequency (defalut: 16000Hz);
%
% stft_win_len - window length of stft (default: 1024 samples or 64ms at 16 kHz.
%
% stft_win_overlap - overlap length of window of stft (default 
% 3/4 * stft_win_len).
%
% cpsn - component per source number (default: 10)
%
% Output:
%
% se - nsrc x nsampl matrix containing the estimated time domain source
% signals.
% 
% ie - nsrc x nsampl x nchan time domain spatial source images
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% (yhebh@connect.ust.hk)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% This code is adapted from the code written by Alexey Ozerov in 2010.
%
% If you use this code, please cite the paper 
%
% A. Ozerov and C. Fevotte,
% "Multichannel nonnegative matrix factorization in convolutive mixtures for audio source separation,"
% IEEE Trans. on Audio, Speech and Lang. Proc. special issue on Signal Models and Representations
% of Musical and Environmental Sounds, vol. 18, no. 3, pp. 550-563, March 2010.
% Available: http://www.irisa.fr/metiss/ozerov/Publications/OzerovFevotte_IEEE_TASLP10.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
[nchan, nsampl] = size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if nchan>2, x=x(1:2,:); warning('Only two channels of signals will be used'); end
if nargin<3, fs=16000; stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; cpsn=10; end
if nargin<4, stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; cpsn=10; end
if nargin<5, stft_win_overlap=0.75*stft_win_len; cpsn=10; end
if nargin<6, cpsn=10; end
%%% Main %%%
nsampl = size(x,2);
X=multi_stft_ltfat(x,stft_win_len,stft_win_overlap);

fprintf('INITIALIZATION (using convolutive filters estimation and source separation via binary masking)\n\n');

% Estimation of the frequency-dependent mixing matrix
fprintf('Estimation of the frequency-dependent mixing matrix\n');
Ae=estmix_conv(X,nsrc,fs);

% Source separation via binary masking
fprintf('Source separation via binary masking\n');
Se_bm=sep_binmask(X,Ae);

fprintf('Source separation via multichannel NMF EM algorithm\n\n');

fprintf('Parameters initialization\n\n');

A_init = Ae;

% initialize W and H from bm separated sources
[W_init, H_init, source_NMF_ind] = Init_KL_NMF_fr_sep_sources(Se_bm, cpsn);

% initialize additive noise variances as mixture PSD / 100 
mix_psd = 0.5 * (mean(abs(X(:,:,1)).^2 + abs(X(:,:,2)).^2, 2));
Sigma_b_init = mix_psd / 100;

% run 200 iterations of multichannel NMF EM algorithm (with annealing and noise injection)
A_init = permute(A_init, [3 1 2]);

[W_EM, H_EM, Ae_EM, Sigma_b_EM, Se_EM, log_like_arr] = ...
    multinmf_conv_em(X, W_init, H_init, A_init, Sigma_b_init, source_NMF_ind, 200, 2);

Ae_EM = permute(Ae_EM, [2 3 1]);
% Computation of the source signals
se_EM=istft_ltfat(Se_EM,nsampl,stft_win_len,stft_win_overlap);
% Computation of the spatial source images
fprintf('Computation of the spatial source images\n');
Ie_EM=src_image(Se_EM,Ae_EM);
ie_EM=multi_istft_ltfat(Ie_EM,nsampl,stft_win_len,stft_win_overlap);

%Return
se=se_EM;
ie=ie_EM;
end
