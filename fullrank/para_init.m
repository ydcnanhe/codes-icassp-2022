
function [A,R]=para_init(X,fs,J,d,c,K)

% Mixing matrix/spatial covariance matrix estimation by hierarchical
% clustering in the T-F domain
%
% Inputs:
% X: F x N x I STFT mixture signal
% fs: sampling rate
% J: number of sources
% d: microphone distance
% c: sound velocity
% K: number of clusters
% 
% Outputs:
% A: I x J x F estimated mixing matrix
% R: I x I x F x J estimated spatial covariance matrix
%
%**************************************************************************
% Copyright 2010
% Author Ngoc Q. K. Duong
% This source code is distributed under the terms of the GNU General Public
% License version 3 (http://www.gnu.org/licenses/gpl.txt)
% If you find it useful, please cite the following reference:
% N. Q. K. Duong, E. Vincent and R. Gribonval, "Under-determined reverberant
% audio source separation using a full-rank spatial covariance model",
% IEEE Transactions on Audio, Speech and Language Processing, 2010
%**************************************************************************

% ----- Errors and warnings
if nargin<6, error('Not enough input arguments.'); end
if J<0, error('The number of sources must be larger 0.'); end

[F,N,I]=size(X);
XX=permute(X,[3,2,1]); X0=XX;

% ----- normalization
for f=1:F,
    for t=1:N
        clear i;  
        XX(:,t,f)=XX(:,t,f)*exp(-i*angle(XX(1,t,f)))/norm(XX(:,t,f));
        X0(:,t,f)=X0(:,t,f)*exp(-i*angle(X0(1,t,f)));
    end
end
XX=permute(XX,[2,1,3]); X0=permute(X0,[2,1,3]); 
      
% ----- hierachical clustering   
R=zeros(I,I,F,J);
Rt=zeros(I,I,F,N,J);
for f=1:F
    ind=0;
    for t=1:N-1
        for tt=t+1:N
            ind=ind+1;  Y(ind)=norm(XX(t,:,f)-XX(tt,:,f));
        end
    end            
    Z = linkage(Y,'average'); 
    T = cluster(Z,'maxclust',K);
        
    for k=1:K,  C(k)=length(find(T==k));  end
    for j=1:J
        [maxj,ind]=max(C);   C(ind)=0;  
        Gj=find(T==ind);     NN(f,j)=maxj;
        A(:,j,f)=sum(X0(Gj,:,f))/maxj;
        temp=zeros(I,I);             
        for k=1:maxj, 
            temp=temp+reshape(X0(Gj(k),:,f),I,1)*reshape(X0(Gj(k),:,f),I,1)'; 
            Rt(:,:,f,k,j)=reshape(X0(Gj(k),:,f),I,1)*reshape(X0(Gj(k),:,f),I,1)';
        end
        R(:,:,f,j)=temp/maxj;                     
    end
end
     

 