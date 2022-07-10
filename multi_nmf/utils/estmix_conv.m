function Af_est = estmix_conv(Y,nsrc,fs)
%
% Description:
%
% estmix_conv - estimate the mixing matrix by duet clustering
% assuming convolutive mixing

% The following problem is solved
% Given STFT coefficients matrix Y, Y(f,t)=Af(f)X(f,t), estimate the mixing 
% matrix Af.
%
% Syntax:
% 
% Af_est = estmix_conv(Y,nsrc)
% Af_est = estmix_conv(Y,nsrc,fs)
%
% Input:
% Y - nbin x nfram x nchan matrix containing stft coefficients of mixed
% signals

% nsrc - the number of estimated signals.

% fs - sampling frequency (default: 16000Hz).

% Output:
% Af_est -  nchan X nsrc X nbin matrix containing frequency dependent 
% mixing parameters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% (yhebh@connect.ust.hk)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
if nargin<3, fs=16000; end
%%% Main %%%
[nbin,~,nchan]=size(Y);
Af_est=zeros(nchan,nsrc,nbin); % initialization of A_est
D=dist(Y,nsrc,fs,0); % obtain distance by duet clustering
for f=1:nbin
    Yf=reshape(Y(f,:,:),[],2);
    Yf=Yf./vecnorm(Yf,2,2).*exp(-1i*angle(Yf(:,1))); % phase and amplitude normalization
    Df=squeeze(D(f,:,:));
    Df=Df./mean(Df(:));
    threshold=inf;
    [val,loc]=min(Df,[],2);
    T = loc .* (val<threshold); % 0 represents MSP;
    Aff_est=zeros(nchan,nsrc);
    for n=1:nsrc
            Aff_est(:,n)=mean(Yf(T==n,:)); % normal mean
    end
    % normalize mixing vector
    Aff_est(isnan(Aff_est))=1;
    Aff_est=Aff_est./vecnorm(Aff_est,2,1);
    Af_est(:,:,f)=Aff_est;
end