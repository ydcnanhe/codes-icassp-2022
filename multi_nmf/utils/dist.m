function D = dist(X,nsrc,fs,isPlot_hist)
%
% Description:
%
% dist - calculate the distance between local mixing parameters and 
% estimated global mixing parameters.

% Syntax:
% 
% D = dist(X,nsrc)
% D = dist(X,nsrc,fs)
% D = dist(X,nsrc,fs,isPlot_hist)

% Inputs:
%
% X - stft coefficients of mixed signals. Format: Nf X Nt X M complex 
% matrix, where Nf is the number of frequency bins, Nt is the number of 
% frames, and M is the number of channels/sensors.
%
% nsrc: the number of sources.
%
% fs: frequency of sampling (default: 16000Hz).
%
% isPlot_hist: true/false, whether to plot the delay-attenuation histogram
% (default: false).

% Output:

% D: the distance between local mixing parameters and 
% estimated global mixing parameters. Format: Nf X Nt X N tensor.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2022 Yudong He
% (yhebh@connect.ust.hk)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Errors and warnings %%%
if nargin<3, error('Not enough input arguments.'); end
if nargin<4, fs=16000; isPlot_hist=0; end
if nargin<5, isPlot_hist=0; end
%%% Main %%%
    X1=X(:,:,1);
    X2=X(:,:,2);
    Wd2=size(X1,1)-1;
    vf=0:fs/(2*Wd2):fs/2;
    vf=vf';
    n_time_bins = size(X1,2);
    wmat = vf * (ones(1,n_time_bins)) * (2 * pi / fs);
    wmat = wmat + eps;
    % compute_atn_delay
    inter_channel_ratio = (X2 + eps) ./ (X1 + eps);
    attenuation = abs(inter_channel_ratio);
    symmetric_atn = attenuation - 1 ./ attenuation; % symmetric_attenuation
    delay = -1*imag(log(inter_channel_ratio)) ./ wmat;% relative delay. unit: sample
    % parameters for duet
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
    % compute distance
    D=ones([size(X1), nsrc]);
    for i = 1:nsrc
        phase = exp(-1j *wmat * delay_peak(i));
        D(:,:,i) = abs(atn_peak(i) * phase .* X1 - X2) .^ 2 / (1 + atn_peak(i).^ 2);
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
        colorbar;
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