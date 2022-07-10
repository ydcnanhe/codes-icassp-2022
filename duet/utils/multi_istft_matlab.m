function y=multi_istft_matlab(X,nsampl,stft_win_len,stft_win_overlap)
%
% Description:
%
% multich_istft_matlab: Multichannel inverse short-time Fourier transform (ISTFT)
% using hann windows.
%
% Syntax:
%
% y=multi_istft_matlab(X,nsampl)
% y=multi_istft_matlab(X,nsampl,stft_win_len)
% y=multi_istft_matlab(X,nsampl,stft_win_len,stft_win_overlap)
% 
% Input:
%
% X - nbin x nfram x nsrc tensor containing STFT coefficients for nsrc 
% sources with nbin frequency bins and nfram time frames or nbin x nfram x
% nsrc x nchan matrix containing the STFT coefficients for nsrc spatial
% source images over nchan channels.
%
% stft_win_len - window length of stft (default: 1024 samples or 64ms at 
% 16 kHz.
%
% stft_win_overlap - overlap length of window of stft. (default 
% 3/4-overlapping).
%
% nsampl - # of samples of the corresponding time domain signals.
%
% Output:
%
% y - nsrc x nsampl matrix or nsrc x nsampl x nchan tensor containing the
% corresponding time-domain signals
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
if nargin<2, error('Not enough input arguments.'); end
[nbin,nfram,nsrc,nchan]=size(X);nfft=(nbin-1)*2;
if nargin<3, stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; end
if nargin<4,  stft_win_overlap=0.75*stft_win_len; end

%%% Computing inverse STFT signal %%%
y=zeros(nsrc,nsampl,nchan);
for m=1:nchan
    Y=zeros(nfft,nfram,nsrc);
    Y(nfft/2:end,:,:) = X(:,:,:,m);
    Y(nfft/2-1:-1:1,:,:) = conj(X(2:nfft/2,:,:,m));
    for n = 1:nsrc
        tmp = istft(Y(:,:,n),'Window',hann(stft_win_len,'periodic'),...
        'FFTLength',nfft,'OverlapLength',stft_win_overlap);
        tmp(1:50)=0;% sometimes there are abnormal numbers at first 50 samples.
        tmp(end:-1:end-50)=0;% sometimes there are abnormal numbers at last 50 samples.
        % normalization
        tmp=real(tmp);
        tmp = tmp - mean(tmp);
        tmp = tmp ./ max(abs(tmp));
        if length(tmp)<nsampl
            y(n,1:length(tmp),m)=tmp;
        else
            y(n,1:nsampl,m)=tmp(1:nsampl);
        end
    end
end

return;