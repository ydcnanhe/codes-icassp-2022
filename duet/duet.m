function se = duet(x,nsrc,fs,stft_win_len,stft_win_overlap,isPlot_hist)
%
% Description:
%
% duet: Implementation of the degenerate unmixing estimation technique
% (Rickard, 2007), a multiple sources blind separation with two microphones. 
% Separation by binary masks that are constructed by spatial cues clustering.
%
% Syntax:
% 
% se = duet(x,nsrc)
% se = duet(x,nsrc,fs)
% se = duet(x,nsrc,fs,stft_win_len)
% se = duet(x,nsrc,fs,stft_win_len,stft_win_overlap)
% se = duet(x,nsrc,fs,stft_win_len,stft_win_overlap,isPlot_hist)
%
% Inputs:
%
% x - 2 x nsampl matrix containg two chanels of time domain mixture signals with
% nsampl samples.
%
% nsrc - number of sources.
%
% fs - sampling frequency (defalut: 16000Hz);
%
% stft_win_len - stft window length (default: 1024 samples or 64ms at 16 kHz.
%
% stft_win_overlap - overlap length of window of stft (default 
% 3/4 * stft_win_len).
%
% isPlot_hist: true/false, whether to plot the delay-attenuation histogram
% (default: false).
%
% Output:
%
% se - nsrc x nsampl matrix containing the estimated time domain source
% signals.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% (yhebh@connect.ust.hk)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% This code is rewritten from a python version provided in "nussl" package
%
% If you use this code, please cite these papers  
%
% @incollection{rickard2007duet,
%   title={The DUET blind source separation algorithm},
%   author={Rickard, Scott},
%   booktitle={Blind speech separation},
%   pages={217--241},
%   year={2007},
%   publisher={Springer}
% }
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
%
%%% Errors and warnings %%%
if nargin<2, error('Not enough input arguments.'); end
[nchan, nsampl] = size(x);
if nchan>nsampl, error('The signals must be within rows.'); end
if nchan>2, x=x(1:2,:); warning('Only two channels of signals will be used'); end
if nargin<3, fs=16000; stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; isPlot_hist=0; end
if nargin<4, stft_win_len=1024; stft_win_overlap=0.75*stft_win_len; isPlot_hist=0; end
if nargin<5, stft_win_overlap=0.75*stft_win_len; isPlot_hist=0; end
if nargin<6, isPlot_hist=0; end
%%% Main %%%
    X=multi_stft_matlab(x,stft_win_len,stft_win_overlap); % STFT
    [nbin,nfram,~]=size(X); Wd2=nbin-1;
    X1=X(:,:,1);
    X2=X(:,:,2);
    vf=0:fs/(2*Wd2):fs/2;
    vf=vf';
    wmat = vf * (ones(1,nfram)) * (2 * pi / fs);
    wmat = wmat + eps;
    % compute_atn_delay
    inter_channel_ratio = (X2 + eps) ./ (X1 + eps);
    attenuation = abs(inter_channel_ratio);
    symmetric_atn = attenuation - 1 ./ attenuation; % symmetric_attenuation
    delay = -1*imag(log(inter_channel_ratio)) ./ wmat;% relative delay. unit: sample
    % parameters for histogram
    p = 1; q=0;
    attenuation_min = -3; attenuation_max = 3; 
    delay_min = -3; delay_max = 3;
    num_attenuation_bins = 50; num_delay_bins=50;
    attenuation_min_distance = 5; delay_min_distance = 5;
    peak_threshold = 0.0;
    % calculate the weighted histogram
    time_frequency_weights = (abs(X1) .* abs(X2)) .^p .* ...
                                (abs(wmat)) .^ q;
    % only consider time-freq. points yielding estimates in bounds                 
    attenuation_premask = (attenuation_min < symmetric_atn) & (symmetric_atn < attenuation_max);
    delay_premask = (delay_min < delay) &(delay < delay_max);
    attenuation_delay_premask = attenuation_premask & delay_premask;
    nonzero_premask = (attenuation_delay_premask~=0);
    symmetric_attenuation_vector = symmetric_atn(nonzero_premask);
    delay_vector = delay(nonzero_premask);
    time_frequency_weights_vector = time_frequency_weights(nonzero_premask);
    bins_array = [num_attenuation_bins,num_delay_bins];
    range_array = [attenuation_min,attenuation_max;delay_min,delay_max];
    % compute the histogram
    [histogram, atn_bins, delay_bins] = histogram2d(symmetric_attenuation_vector, ...
        delay_vector,bins_array,range_array,time_frequency_weights_vector,isPlot_hist);
    % find peak indices
    min_dist=[attenuation_min_distance, delay_min_distance];% only be a 1x2 vector
    [pks,peak_indices] = find_peak_indices(histogram,nsrc,peak_threshold,min_dist);
    % compute delay_peak, attenuation peak, and attenuation/delay estimates
    atn_indices = peak_indices(1,:);
    delay_indices = peak_indices(2,:);
    symmetric_atn_peak = atn_bins(atn_indices);
    delay_peak = delay_bins(delay_indices);
    atn_peak = (symmetric_atn_peak + sqrt(symmetric_atn_peak.^ 2 + 4)) ./ 2;
    % compute masks
    best_so_far = ones(size(X1)) * inf;
    result_masks = true([size(X1), nsrc]);
    for i = 1:nsrc
        mask_array = zeros(size(X1));
        phase = exp(-1j *wmat * delay_peak(i));
        score = abs(atn_peak(i) * phase .* X1 - X2) .^ 2 / (1 + atn_peak(i).^ 2);
        mask = (score < best_so_far);
        mask_array(mask) = true;
        background_mask = mask_array;
        result_masks(:,:,i) = background_mask;
        result_masks(:,:,1) = xor(result_masks(:,:,i), result_masks(:,:,1));
        best_so_far(mask) = score(mask);
    end
    result_masks(:,:,1) = not(result_masks(:,:,1));
    % applying mask
    Se = apply_mask(X1,result_masks);
    % convert to time domain
    se=multi_istft_matlab(Se,nsampl,stft_win_len,stft_win_overlap);
    % normalization
    for n=1:nsrc
        tmp=se(n,:);
        tmp = tmp - mean(tmp);
        se(n,:) = tmp ./ max(abs(tmp));
    end
end

function  [histogram, atn_bins, delay_bins] = histogram2d(alpha_vec, delta_vec, bins_array, range_array, tfweight,isPlot)
    abins = bins_array(1);
    dbins = bins_array(2);
    mina = range_array(1,1);
    maxa = range_array(1,2);
    mind = range_array(2,1);
    maxd = range_array(2,2);
    alpha_ind = round(1+(abins-1)*(alpha_vec-mina)/(maxa-mina)); 	
    delta_ind = round(1+(dbins-1)*(delta_vec-mind)/(maxd-mind));
    histogram=accumarray([alpha_ind delta_ind],tfweight,[abins dbins]);
    % scale histogram from 0 to 1
    histogram = histogram ./ max(histogram(:));
    % smooth the normalized histogram - local average 3-by-3 neighboring bins
    Kernel = ones(3,3);
    Kernel = Kernel ./ (3*3);
    copy_row = floor(size(Kernel,1) / 2);  % number of rows to copy on top and bottom
    copy_col = floor(size(Kernel,2) / 2);  % number of columns to copy on either side
    % augment matrix to be smoothed
    aug_histogram = zeros(2*copy_row+size(histogram,1),2*copy_col+size(histogram,2));
    aug_histogram(1+copy_row:copy_row+size(histogram,1),1+copy_col:copy_col+size(histogram,2)) = histogram;
    histogram = conv2(aug_histogram, Kernel,'valid');
    atn_bins = mina + (maxa-mina)/(abins-1)*(0:1:abins-1);
    delay_bins = mind + (maxd-mind)/(dbins-1)*(0:1:dbins-1);
    if isPlot==true
        aaxis=linspace(mina,maxa,abins);
        daxis=linspace(mind,maxd,dbins);
        figure()
        surf(daxis, aaxis,histogram);
        ylabel('symmetric attenuation');
        xlabel('relative delay');
    end
end

function [pks,locs] = find_peak_indices(input_array,n_peaks,threshold,min_dist)
    % throw out everything below threshold
    input_array = input_array.* (input_array >= threshold);
    % check to make sure we didn't throw everything out
    if sum(input_array~=0) == 0
        error('Threshold set incorrectly. No peaks above threshold.')
    end
    if sum(input_array~=0) < n_peaks
        warning('Threshold set such that there will be less peaks than n_peaks.')
    end
    pks = zeros(1,n_peaks);
    locs = zeros(2,n_peaks);
    for i = 1:n_peaks
        [pk,loc] = max(reshape(input_array,1,[]));
        pks(i) = pk;
        [row,col] = ind2sub(size(input_array),loc);
        locs(1,i)=row; locs(2,i)=col;
        lower1 = (locs(1,i)-min_dist(1));
        upper1 = locs(1,i)+min_dist(1)+1;
        lower2 = locs(2,i)-min_dist(2);
        upper2 = locs(2,i)+min_dist(2)+1;
        if lower1 < 1 
            lower1 = 1;
        end
        if lower2 < 1 
            lower2 = 1;
        end
        if  upper1 >= size(input_array,1)
            upper1 = size(input_array,1); 
        end
        if  upper2 >= size(input_array,2)
            upper2 = size(input_array,2); 
        end
        input_array(lower1:upper1,lower2:upper2)=0;
    end
end

function masked_stft = apply_mask(specg_1,masks)
    nscr = size(masks,3);
    magnitude = abs(specg_1);
    phase =  angle(specg_1);
    masked_stft = zeros(size(masks));
    for i = 1: nscr
        masked_abs = magnitude .* squeeze(masks(:,:,i));
        masked_stft(:,:,i) = masked_abs .* exp(1j * phase);
    end
end
    