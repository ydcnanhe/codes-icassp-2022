function y=istft_ltfat(X,nsampl,stft_win_len,stft_win_overlap)
%
% Description:
% istft - inverse stft coefficients. ltfat has to be started first
%
% Syntax:
%
% y=stft_ltfat(x,nsampl)
% y=stft_ltfat(x,nsampl,stft_win_len)
% y=stft_ltfat(x,nsampl,stft_win_len,stft_win_overlap)
% 
% Input:
%
% x - input stft coefficients. Format: Nf X Nt X N, where Nf is the number
% of frequency bins, Nt is the number of time frames, and N is the number
% of signals.

% nsampl - # of samples of time domain singal.
%
% stft_win_len - window length of stft (default: 1024 samples or 64ms at 16 kHz.
%
% stft_win_overlap - overlap length of window of stft (default 
% 3/4-overlapping).
%
% Output:
%
% y - time domain signals. Format: N X T where T is the number of time
% domain samples per signal.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% (yhebh@connect.ust.hk)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
if nargin<3, stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; end
if nargin<4, stft_win_overlap=0.75*stft_win_len; end
%%% Main
    % Set the STFT and iSTFT parameters 
    a = stft_win_len-stft_win_overlap;
    g = gabwin({'tight', 'hann'}, a, stft_win_len);
    synthesis=@(x) idgtreal(x,g,a,stft_win_len,nsampl);  
    N=size(X,3);
    tr=synthesis(X(:,:,1));
    nsampl=size(tr,1);
    y=zeros(N,nsampl);
    for n=1:N
        tmp=synthesis(X(:,:,n));
        y(n,:)=tmp';
    end
end