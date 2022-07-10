
function [ie,crit]=bss_fullrank(x,nsrc,d,fs,EMIter,stft_win_len,stft_win_overlap,K)    
%
% Description:
% 
% bss_fullrank: Blind reverberant source separation using full-rank spatial
% covariance model
%
% Syntax:
%
% [ie,crit]=bss_fullrank(x,nsrc,d)  
% [ie,crit]=bss_fullrank(x,nsrc,d)  
% [ie,crit]=bss_fullrank(x,nsrc,d,fs,)  
% [ie,crit]=bss_fullrank(x,nsrc,d,fs,EMIter)
% [ie,crit]=bss_fullrank(x,nsrc,d,fs,EMIter,stft_win_len)
% [ie,crit]=bss_fullrank(x,nsrc,d,fs,EMIter,stft_win_len,stft_win_overlap)
% [ie,crit]=bss_fullrank(x,nsrc,d,fs,EMIter,stft_win_len,stft_win_overlap,K)   
%
% Inputs:
%
% x - nchan x nsampl mixture signal.
%
% nsrc - number of sources.
%
% d - microphone distance (unit: meter).
%
% fs: sampling rate (default: 16000Hz).
%
% EMIter - number of EM iterations (default: 20).
%
% stft_win_len - stft window length (default: 1024 samples or 64ms at 16 kHz.
%
% stft_win_overlap - overlap length of window of stft (default 
% 3/4 * stft_win_len).
%
% K - number of cluster for the hierarchical clustering (default:30).
% 
% Outputs:
%
% ie - nsrc x nsampl x nchan estimated source images.
%
% crit - average log-likelihood criterion.
%
%**************************************************************************
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
% N. Q. K. Duong, E. Vincent and R. Gribonval, "Under-determined reverberant
% audio source separation using a full-rank spatial covariance model",
% IEEE Transactions on Audio, Speech and Language Processing, 2010
%**************************************************************************

% ----- Errors and warnings
if nargin<4, error('Not enough input arguments.'); end
[nchan, nsampl] = size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if nargin<5, fs=16000; EMIter=20; stft_win_len=1024; stft_win_overlap=0.75* stft_win_len; K=30;  end
if nargin<6, EMIter=20; stft_win_len=1024; stft_win_overlap=0.75* stft_win_len; K=30;  end
if nargin<7, stft_win_len=1024; stft_win_overlap=0.75* stft_win_len; K=30;  end
if nargin<8, stft_win_overlap=0.75* stft_win_len; K=30;  end
if nargin<9, K=30;  end
% ----- Main
[nchan,nsampl]=size(x);
X=multi_stft_ltfat(x,stft_win_len,stft_win_overlap); clear x;
[nbin,nfram,~]=size(X);

Cx=zeros(nchan,nchan,nbin,nfram);
R=zeros(nchan,nchan,nbin,nsrc);     v=ones(nbin,nfram,nsrc);        

% ----- Compute empirical covariance matrix from the mixture
lf=1; lt=1;
winf=hanning(2*lf-1);   wint=hanning(2*lt-1).';
for f=1:nbin
    for t=1:nfram
        indf=max(1,f-lf+1):min(nbin,f+lf-1);
        indt=max(1,t-lt+1):min(nfram,t+lt-1);
        nind=length(indf)*length(indt);
        wei=reshape(winf(indf-f+lf)*wint(indt-t+lt),1,nind); 
        X1=reshape(X(indf,indt,:),nind,nchan).';                
        Cx(:,:,f,t)=(X1.*(ones(nchan,1)*wei))*X1'/sum(wei);
    end
end

% ----- Initialization by hierarchical clustering
fprintf('\nParameter initialization...');
[A,R]=para_init(X,fs,nsrc,d,343,K);

% ----- EM algorithm
fprintf('\nEM iterations...');
Af=[];  Cc=zeros(nchan*nsrc,nchan*nsrc);
for j=1:nsrc, Af=[Af,diag(ones(nchan,1))]; end
s=ones(nbin,nfram,nsrc);           
c=ones(nchan*nsrc,nbin,nfram);
CovX=zeros(nchan,nchan,nbin,nfram);    
       
for iter=1:EMIter
    for f=1:(nbin)
        tempRj=zeros(nchan,nchan,nsrc);
        for t=1:nfram                   
            % E-step
            for j=1:nsrc
                Cc((j-1)*nchan+1:(j*nchan),(j-1)*nchan+1:(j*nchan))=v(f,t,j)*R(:,:,f,j);
            end
            Rx=Af*Cc*Af';          
            G=Cc*Af'*inv(Rx);         
            Rc=G*Cx(:,:,f,t)*G'+(diag(ones(nchan*nsrc,1))-G*Af)*Cc;

            % M-step
            for j=1:nsrc
                subind=(j-1)*nchan+1:(j*nchan);
                v(f,t,j)=real((1/nchan)*trace(inv(R(:,:,f,j))*Rc(subind,subind)));
                tempRj(:,:,j)=tempRj(:,:,j)+(1/v(f,t,j))*Rc(subind,subind);
            end

            % log-likelihood
            MLcrit(f,t)=abs(-log(det(Rx))-reshape(X(f,t,:),nchan,1)'*inv(Rx)*reshape(X(f,t,:),nchan,1));
        end
        for j=1:nsrc, R(:,:,f,j)=tempRj(:,:,j)/nfram;  end
    end
    crit(iter)=sum(sum(MLcrit))/((nbin)*nfram);
end

fprintf('\nSolving the permutation problem...');
% ----- Perform PCA on Rj
for f=1:nbin
    tempW=[];
    for j=1:nsrc
        [E,D]=eig(R(:,:,f,j));
        [maxd,ind]=max(diag(abs(D)));
        W(:,j,f)=E(:,ind);
    end
end
% ----- Align source orders
[W,permu] = permutation(W,fs,d,343);
for f=1:nbin
    tempR(:,:,f,:)=R(:,:,f,permu(:,f));
    tempv(f,:,:)=v(f,:,permu(:,f));
end
R=tempR; 
v=tempv;

% ----- Wiener filtering        
for f=1:(nbin)
    for t=1:nfram
        vR=zeros(nchan,nchan); vR1=zeros(nchan,nchan);
        for j=1:nsrc,  vR=vR+v(f,t,j)*R(:,:,f,j); end
        for j=1:nsrc,  Ie(:,f,t,j)=v(f,t,j)*R(:,:,f,j)*(vR^-1)*reshape(X(f,t,:),nchan,1); end
    end
end

% ----- Inverse time-frequency transformation
Ie=permute(Ie,[2,3,4,1]);    
ie=multi_istft_ltfat(Ie,nsampl,stft_win_len,stft_win_overlap); 

    
