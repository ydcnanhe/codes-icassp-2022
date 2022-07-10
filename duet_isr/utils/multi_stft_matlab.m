function Y=multi_stft_matlab(x,stft_win_len,stft_win_overlap)
%
% Description:
%
% mulich_stft_matlab: Multichannel short-time Fourier transform (STFT) using
% hann windows. STFT function is defined by MATLAB.
%
% Syntax:
%
% Y=multi_stft_matlab(x)
% Y=multi_stft_matlab(x,stft_win_len)
% Y=multi_stft_matlab(x,stft_win_len,stft_win_overlap)
%
% Input:
%
% x - nchan x nsampl matrix containing nchan channels of signals, each with
% nsampl samples.
%
% stft_win_len - window length of stft (default: 1024 samples or 64ms at 
% 16 kHz.
%
% stft_win_overlap - overlap length of window of stft (default 
% 3/4f-overlapping).
%
% Output:
%
% y - nbin x nfram x nchan tensor containing stft coefficients of x with nbin
% frequency bins and nfram time frames.
%
% Remark: In matlab, mapping from frequency to discrete point is:
% -(N/2)(fs/(N-2)):fs/(N-2):fs/2 --> 1:1:N, where N is the number of
% descrete Fourier point (Implicitly assume N is an even number)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Errors and warnings %%%
if nargin<1, error('Not enough input arguments.'); end
if nargin<2, stft_win_len=1024; end
[nchan,nsampl]=size(x); nfft=stft_win_len;
if nchan>nsampl, error('The signals must be within rows.'); end
if nargin<3, stft_win_overlap=0.75*stft_win_len; end

%%% Computing STFT coefficients %%%
tmp= stft(x(1,:),'Window',hann(stft_win_len,'periodic'),...
        'FFTLength',nfft,'OverlapLength',stft_win_overlap);
[~,nfram]=size(tmp);
Y=zeros(nfft,nfram,nchan);
Y(:,:,1)=tmp;
for m = 2:nchan
    Y(:,:,m) = stft(x(m,:),'Window',hann(stft_win_len,'periodic'),...
        'FFTLength',nfft,'OverlapLength',stft_win_overlap);
end
nbin=nfft/2+1;
Y = Y(nbin-1:end,:,:);

return;