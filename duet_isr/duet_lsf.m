function se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu,resp_na)

% DUET_LSF Blind source separation by DUET_LSF algorithm.
% 
% se = duet_lsf(x,nsrc)
% se = duet_lsf(x,nsrc,fs)
% se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap)
% se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu)
% se = duet_lsf(x,nsrc,fs,stft_win_len,stft_win_overlap,mu,resp_na)
%
% Input:
% x: nchan x nsampl matrix containg nchan time domain mixture signals with
% nsampl samples.
% nsrc: number of source signals.
% fs: sampling frequency (defalut: 16000Hz);
% stft_win_len: window length of stft (default: 1024 samples or 64ms at 16 kHz.
% stft_win_overlap: overlap length of window of stft (default 
% half-overlapping).
% mu: the closeness threshold to determine a SSP (default: 0.05).
% resp_na: name of spatial linear filter response (defalut: 'isr').
%
% Output:
% se: nsrc x nsampl matrix containing the corresponding time domain source
% signals.
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


%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
[nchan, nsampl] = size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if nargin<3, fs=16000; stft_win_len=1024; stft_win_overlap=0.5*stft_win_len; mu=0.05; resp_na='isr'; end
if nargin<4, error('Lack of setting of STFT window'); end
if nargin<5, mu=0.05; resp_na='isr'; end
if nargin<6, resp_na='isr'; end

%%% Blind source separation %%%
%% STFT
X=mystft_matlab(x,stft_win_len,stft_win_overlap); % STFT
[nbin,nfram,~]=size(X);
Se = zeros(nbin,nfram,nsrc);
%% find ssp
[ssp_masks,rho] = find_ssp(X,nsrc,fs,mu);
ssp_stft = zeros(nbin,nfram,nchan,nsrc);
for n = 1:nsrc
    for m=1:nchan
        ssp_stft(:,:,m,n)=X(:,:,m).*ssp_masks(:,:,n);
    end
end
%% Calculate responses
responses =zeros(nsrc,nbin,nchan);
for n = 1:nsrc
    Cn = squeeze(ssp_stft(:,:,:,n));
    ind_intf =1:nsrc;
    ind_intf(n)=[];
    intf_stft = squeeze(ssp_stft(:,:,:,ind_intf(1)));
    for i = 2:nsrc-1
        intf_stft = intf_stft+squeeze(ssp_stft(:,:,:,ind_intf(i)));
    end
    responses(n,:,:) = Calculate_responses(intf_stft,Cn,resp_na);
end
%% Linear filtering
for f = 1:nbin
    for t = 1:nfram
        [~,loc] = sort(rho(f,t,:));
        n = loc(1);
        Se(f,t,n) = responses(n,f,1)*X(f,t,1)+...
            responses(n,f,2)*X(f,t,2);
    end
end
%% Back to the time domain
se=myistft_matlab(Se,stft_win_len,stft_win_overlap,nsampl);
% normalization
for n=1:nsrc
    tmp=se(n,:);
    tmp = tmp - mean(tmp);
    se(n,:) = tmp ./ max(abs(tmp));
end
end

function responses = Calculate_responses(intf_stft,Cn,resp_na)
    
    nbin = size(intf_stft,1);
    nchan = size(intf_stft,3);
    responses=zeros(nbin,nchan);
    for  f = 1:nbin
        w = zeros(nchan,1);
        if sum(abs(intf_stft(f,:,1)))<eps
            w(1)=1;
        end
        
        if sum(abs(Cn(f,:,1)))<eps && sum(abs(intf_stft(f,:,1)))>eps
            a_1 = transpose(intf_stft(f,:,1));
            a_2 = transpose(intf_stft(f,:,2));
            Arg_x = -1*angle(a_1'*a_2);
            Amp_x = -1*abs(a_1'*a_2)/(a_2'*a_2);
            w(1) = 1;
            w(2) = Amp_x.*exp(1i.*Arg_x);
            w=w';
        end
        
        if sum(abs(Cn(f,:,1)))>eps && sum(abs(intf_stft(f,:,1)))>eps
            % choose one type response
                % ISR
            if resp_na == "isr"
                a_1 = transpose(intf_stft(f,:,1));
                a_2 = transpose(intf_stft(f,:,2));
                Arg_x = -1*angle(a_1'*a_2);
                Amp_x = -1*abs(a_1'*a_2)/(a_2'*a_2);
                w(1) = 1;
                w(2) = Amp_x.*exp(1i.*Arg_x);
                w=w';
            end
                % MVDR
            if resp_na == "mvdr"
                af = ones(nchan,1);
                x2 = squeeze(Cn(f,:,2));
                x1 = squeeze(Cn(f,:,1));
                x2(x2==0)=[];
                x1(x1==0)=[];
                af(2) = mean(x2./x1);
                af=af./norm(af);
                uf = squeeze(intf_stft(f,:,:)); % Nt X Ch
                uf = uf(1:find(uf==0,1),:); 
                Nt = size(uf,2); % eff # of time frames
                uf = transpose(uf); % Ch X Nt
                Sigma_u = uf*uf'/Nt; % sample-covariance matrix
                inv_Sigma_u = (Sigma_u+eps*eye(nchan))^(-1);% inverse sample-covariance matrix
                mu = 0;
                w = inv_Sigma_u * af / (mu+af'*inv_Sigma_u*af);
            end
        end
        responses(f,:)=w';
    end
end