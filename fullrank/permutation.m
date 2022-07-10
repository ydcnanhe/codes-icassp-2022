
function [A,permu] = permutation(A,fs,d,c)

% Solve the permutation problem

% Reference: H. Sawada et.al., "Grouping Separated Frequency Components by
% Estimating Propagation Model Parameters in Frequency-domain Blind Source
% Separation", IEEE Trans. on Audio, Speech, and Language Processing, 2007, pp 1592-1604

% Inputs/Outputs:
%
% A: I x J x F estimated mixing matrix 
% fs: sampling rate
% J: number of sources
% d: microphone distance
% c: sound velocity
%
% permu: J x F aligned order of sources across all frequencies

%**************************************************************************
% Copyright 2010
% Author Ngoc Q. K. Duong
% This source code is distributed under the terms of the GNU General Public
% License version 3 (http://www.gnu.org/licenses/gpl.txt)
%**************************************************************************
[I,J,Nh]=size(A);
N = 2*(Nh-1);
order=zeros(J,Nh);

A1 = zeros(I,J,Nh);
A2 = zeros(I,J,Nh);   
Jr = 1;                % reference sensor for normalization
d = d*(I-Jr);          % maximum distance between microphones
fL = floor(c/(2*d)); 
KL = floor((fL/fs)*N); % low-frequency bin without spatial aliasing
K = min(KL,Nh);   
permu_list=perms(1:J);

% ----- Normalization 
clear i;
for k=1:Nh    
    for j=1:J, 
        A1(:,j,k) = A(:,j,k)*exp(-i*angle(A(Jr,j,k))) / norm(A(:,j,k));
    end
    A2(:,:,k) = abs(A1(:,:,k)).*exp( i*angle(A1(:,:,k))/(4*(k*fs/N)*d/c) ); 
end

% take the farest mic to maximize angle difference
for j=1:J,  arg_org(j,:) = angle(A1(I,j,:));  end

% ----- For low frequencies without spatial aliasing
% Initialize h
sumh = zeros(I,J);
for j=1:J
    for k=1:min(KL,Nh),  sumh(:,j) = sumh(:,j) + A2(:,j,k);   end
    h(:,j) = sumh(:,j)/norm(sumh(:,j));
end    
        
% ----- k-mean algorithm
NITER = 5; 
for ITER=1:NITER      
    sumh = zeros(I,J);
    for k = 1:min(KL,Nh)         
        for m = 1:length(permu_list)
            order = permu_list(m,:);  JM(m) = 0;
            % check all the permutation of a2
            for j=1:J,  JM(m) = JM(m) + norm(A2(:,order(j),k)-h(:,j));  end
        end
        [minJM,index] = min(JM);
        trueorder = permu_list(index,:);     
        % exchange columns according to true order
        for j=1:J,  sumh(:,j) = sumh(:,j) + A2(:,trueorder(j),k);  end

        if (ITER == NITER)
            temp1 = A1(:,trueorder,k);   A1(:,:,k) = temp1;
            temp  = A(:,trueorder,k);    A(:,:,k) = temp;
            permu(:,k) = trueorder;
        end
    end
    for j=1:J,  h(:,j)=sumh(:,j)/norm(sumh(:,j));  end
end
    
% Compute average value for lamda & tau [see reference paper]
lamda = abs(h);      tau = -2*d*angle(h)/(pi*c);
sum1 = zeros(I,J);   sum2 = zeros(I,J);
clear i;

% ----- For high frequencies with spatial aliasing
if KL<Nh  % if aliasing exist     
    mu = 10^-5;      
    for k = (KL+1):Nh        
        h1 = lamda.*exp(-i*2*pi*(k*fs/N)*tau);
        
        for m = 1:length(permu_list)
            order = permu_list(m,:);  JM(m) = 0;        
            % check all the permutation
            for j=1:J, 
                JM(m) = JM(m) + norm(A1(:,order(j),k)-h1(:,j));
                %JM(m) = JM(m) + sum(abs(angle(A1(:,order(j),k)-h1(:,j)))); 
            end
        end
        [minJM,index] = min(JM);
        trueorder = permu_list(index,:);
        temp1 = A1(:,trueorder,k);   A1(:,:,k)= temp1;
        temp  = A(:,trueorder,k);    A(:,:,k) = temp;
        permu(:,k) = trueorder;
            
        % update lamda(i,l) and tau(i,l)
        sum1 = ((k*fs/N)/fs) * imag(A1(:,:,k).*exp(i*2*pi*(k*fs/N)*tau));
        sum2 = lamda - real(A1(:,:,k).*exp(i*2*pi*(k*fs/N)*tau));       
        tau = tau - mu*lamda.*sum1;   lamda = lamda - mu*sum2;
    end
end
  
% plot angle BEFORE permutation alignment

figure(1);
subplot(2,1,1);
for j=1:J, color(j,:) = rand(1,3); plot(arg_org(j,:),'.','Color',color(j,:)); hold on; end
title('Angle of A_{Ij} before the permutation'); 
xlabel('Frequency bin'); ylabel('Argument');  axis([1 Nh -pi pi]);    

% plot angle AFTER permutation alignment
for j=1:J, arg_per(j,:) = angle(A1(I,j,:)); end
subplot(2,1,2); 
for j=1:J, plot(arg_per(j,:),'.','Color',color(j,:)); hold on;  end
title('Angle of A_{Ij} after the permutation'); 
xlabel('Frequency bin'); ylabel('Argument');  axis([1 Nh -pi pi]);  


