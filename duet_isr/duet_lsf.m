function se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu,lsf_na)
%
% Description:
%
% duet_lsf: Multiple sources blind separation with two microphones.
% Preliminary separation by DUET then post-filtering by linear spatial
% filter (beamformer).
% 
% Syntax:
% 
% se = duet_lsf(x,nsrc)
% se = duet_lsf(x,nsrc,fs)
% se = duet_lsf(x,nsrc,fs,stft_win_len)
% se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap)
% se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu)
% se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu,lsf_na)
%
% Input:
%
% x - 2 x nsampl matrix containg two chanels of time domain mixture signals with
% nsampl samples.
%
% nsrc - number of sources.
%
% fs - sampling frequency (defalut: 16000Hz);
%
% stft_win_len - stft window length (default: 1024 samples or 64ms at 16 kHz).
%
% stft_win_overlap - stft window overlap length (default:
% 3/4 * stft_win_len).
%
% mu - the closeness threshold to determine an SSP (default: 0.05).
%
% lsf_na - name of linear spatial filter response, 'mvdr' or 'isr' (defalut: 'isr').
%
% Output:
%
% se - nsrc x nsampl matrix containing the estimated time domain source
% signals.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% If you use this code, please cite this paper
%
% @inproceedings{he2022harvesting,
%   title={Harvesting Partially-Disjoint Time-Frequency Information for Improving Degenerate Unmixing Estimation Technique},
%   author={He, Yudong and Wang, He and Chen, Qifeng and So, Richard HY},
%   booktitle={ICASSP 2022-2022 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},
%   pages={506--510},
%   year={2022},
%   organization={IEEE}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
[nchan, nsampl] = size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if nchan>2, x=x(1:2,:); warning('Only two channels of signals will be used'); end
if nargin<3, fs=16000; stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; mu=0.05; lsf_na='isr'; end
if nargin<4, stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; mu=0.05; lsf_na='isr';  end
if nargin<5, stft_win_overlap=0.75*stft_win_len; mu=0.05; lsf_na='isr'; end
if nargin<6, mu=0.05; lsf_na='isr'; end
if nargin<7, lsf_na='isr'; end

%%% Main %%%
%% Short-time Fourier Transform
X=multi_stft_matlab(x,stft_win_len,stft_win_overlap); % STFT
[nbin,nfram,~]=size(X);
Se = zeros(nbin,nfram,nsrc);
%% find single source points (SSPs)
[ssp_masks,rho] = find_ssp(X,nsrc,fs,mu);
ssp_stft = zeros(nbin,nfram,nchan,nsrc); % stft coeff at the ssp of nsrcs sources
for n = 1:nsrc
    for m=1:nchan
        ssp_stft(:,:,m,n)=X(:,:,m).*ssp_masks(:,:,n);
    end
end
%% Calculate linear spatial filter
lsf_coeff = zeros(nsrc,nbin,nchan);% lsf_coeff: coefficients of the linear spatial filters corresponding to nsrc sources
for n = 1:nsrc
    tgt_stft = squeeze(ssp_stft(:,:,:,n)); % ssp_tgt_stft: stft coeff of the target signal
    ind_intrf =1:nsrc; % ind_inter: index of interference signals
    ind_intrf(n)=[]; % n-th signal is the target signal and the rest are interference signals
    intrf_stft = squeeze(ssp_stft(:,:,:,ind_intrf(1)));% intrf_stft: stft coeff of the interference signals
    for i = 2:nsrc-1
        intrf_stft = intrf_stft+squeeze(ssp_stft(:,:,:,ind_intrf(i)));
    end
    lsf_coeff(n,:,:) = Calculate_lsf(intrf_stft,tgt_stft,lsf_na);
end
%% Masking and Post-filtering
for f = 1:nbin
    for t = 1:nfram
        [~,loc] = sort(rho(f,t,:));
        n = loc(1);
        Se(f,t,n) = lsf_coeff(n,f,1)*X(f,t,1)+...
            lsf_coeff(n,f,2)*X(f,t,2);  % combing masking and post-filtering
    end
end
%% Back to the time domain
se=multi_istft_matlab(Se,nsampl,stft_win_len,stft_win_overlap);
% normalization
for n=1:nsrc
    tmp=se(n,:);
    tmp = tmp - mean(tmp);
    se(n,:) = tmp ./ max(abs(tmp));
end
end
%%% Internally defined function %%%
%% calculate coefficients of linear spatial filter
function lsf_coeff = Calculate_lsf(intrf_stft,tgt_stft,lsf_na)
%
% Description:
%
% Calculate_lsf: Given stft coeff of the target signal and the interference
% signals, calculate the coefficients of the linear spatial filter to
% extract the target signal
%
% Input:
% 
% intrf_stft - stft coefficients of the interference signals.
% tgt_stft - stft coefficients of the target signal.
% lsf_na - type of linear spatial filter, e.g., minimum variance
% distortionless response (MVDR) or interference suppression response
% (ISR).
%
% Output:
%
% lsf_coeff - coeffficients of linear spatial filters.
% 
% Main:
    nbin = size(intrf_stft,1);
    nchan = size(intrf_stft,3);
    lsf_coeff=zeros(nbin,nchan);
    for  f = 1:nbin
        w = zeros(nchan,1);
        if sum(abs(intrf_stft(f,:,1)))<eps 
            % in the case of inadequate information of the interference
            % signals, no filtering effect
            w(1)=1; 
        end
        
        if sum(abs(tgt_stft(f,:,1)))<eps && sum(abs(intrf_stft(f,:,1)))>eps 
            % in the case of inadequate information of the target signals,
            % using ISR
            a_1 = transpose(intrf_stft(f,:,1));
            a_2 = transpose(intrf_stft(f,:,2));
            Arg_x = -1*angle(a_1'*a_2);
            Amp_x = -1*abs(a_1'*a_2)/(a_2'*a_2);
            w(1) = 1;
            w(2) = Amp_x.*exp(1i.*Arg_x);
            w=w';
        end
        
        if sum(abs(tgt_stft(f,:,1)))>eps && sum(abs(intrf_stft(f,:,1)))>eps  
            % in the case of adequate information of the target and the interference signals
            
            % choose one type response
                % ISR
            if lsf_na == "isr"
                a_1 = transpose(intrf_stft(f,:,1));
                a_2 = transpose(intrf_stft(f,:,2));
                Arg_x = -1*angle(a_1'*a_2);
                Amp_x = -1*abs(a_1'*a_2)/(a_2'*a_2);
                w(1) = 1;
                w(2) = Amp_x.*exp(1i.*Arg_x);
                w=w';
            end
                % MVDR
            if lsf_na == "mvdr"
                af = ones(nchan,1);
                x2 = squeeze(tgt_stft(f,:,2));
                x1 = squeeze(tgt_stft(f,:,1));
                x2(x2==0)=[];
                x1(x1==0)=[];
                af(2) = mean(x2./x1);
                af=af./norm(af);
                uf = squeeze(intrf_stft(f,:,:)); % shape: nfram X nchan
                uf = uf(1:find(uf==0,1),:); 
                Nt = size(uf,2); % # of time frames
                uf = transpose(uf); % shape: nchan X nfram
                Sigma_u = uf*uf'/Nt; % sample-covariance matrix
                inv_Sigma_u = (Sigma_u+eps*eye(nchan))^(-1);% inverse sample-covariance matrix
                mu = 0;
                w = inv_Sigma_u * af / (mu+af'*inv_Sigma_u*af);
            end
        end
        lsf_coeff(f,:)=w';
    end
end