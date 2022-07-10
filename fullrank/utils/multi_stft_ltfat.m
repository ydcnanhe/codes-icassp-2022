function y=multi_stft_ltfat(x,stft_win_len,stft_win_overlap)
%
% Description:
%
% multi_stft_ltfat - stft multichannel signals. ltfat has to be started first

% Syntax:
%
% y=multi_stft_ltfat(x)
% y=multi_stft_ltfat(x,stft_win_len)
% y=multi_stft_ltfat(x,stft_win_len,stft_win_overlap)
%
% Input:
%
% x - input signals. Format: M X T numerical matrix, where M is the number
% of signals and T is the number of samples per signal.
%
% stft_win_len - window length of stft (default: 1024 samples or 64ms at 16 kHz.
%
% stft_win_overlap - overlap length of window of stft (default 
% 3/4-overlapping).

% Output:
% y - stft coefficients of x. Format: Nf X Nt X M, where Nf is the number
% of frequency bins and Nt is the number of time frames. Nf=L_win/2+1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% (yhebh@connect.ust.hk)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Errors and warnings %%%
if nargin<1, error('Not enough input arguments.'); end
[nchan, nsampl] = size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if nargin<2, stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; end
if nargin<3, stft_win_overlap=0.75*stft_win_len; end
%%% Main %%%
    % Set the STFT and iSTFT parameters
    a = stft_win_len-stft_win_overlap;
    g = gabwin({'tight', 'hann'}, a, stft_win_len);
    analysis=@(x) dgtreal(x,g,a,stft_win_len);
    M=size(x,1);
    tr=analysis(x(1,:));
    [Nf,Nt]=size(tr);
    y=zeros(Nf,Nt,M);
    for m=1:M
        y(:,:,m)=analysis(x(m,:));
    end
end