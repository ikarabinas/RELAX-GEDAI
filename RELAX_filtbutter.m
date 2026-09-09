% RELAX_filtbutter() - filters the data using a zero-phase butterworth
%                   filter. Either a band pass, band stop, high pass, or low pass 
%                   filter can be implemented. The filter order is defined by the user.
%                   This function uses the matlab butter and filtfilt functions.
%
% Usage:
%   >>  EEG = RELAX_filtbutter( EEG, high, low, ord, type );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
%   high            - integer (required). Highpass frequency (allow frequencies
%                   above this value with bandpass or high pass)
%                   Example: 1
%   low             - integer (required). Lowpass frequency (allow
%                   frequencies below this value with bandpass or low pass)
%                   Example: 50
%   ord             - integer (required). Filter order.
%                   Example: 4 (designs a fourth order butterworth filter)
%   type            - 'str' (required). 'bandpass' | 'bandstop' | 'highpass' | 'lowpass'. Designs either
%                   a zero-phase band pass, band stop, high pass, or low pass butterworth filter
%    
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = RELAX_filtbutter( EEG, 1, 100, 4, 'bandpass' ); %zero-phase, 4th-order bandpass butterworth filter between 1-100 Hz
%   EEG = RELAX_filtbutter( EEG, 48, 52, 4, 'bandstop' ); %zero-phase, 4th-order bandstop butterworth filter between 48-52 Hz
%   EEG = RELAX_filtbutter( EEG, 1, [], 4, 'highpass' ); %zero-phase, 4th-order high pass butterworth filter allowing frequencies above 1 Hz
%   EEG = RELAX_filtbutter( EEG, [], 45, 4, 'lowpass' ); %zero-phase, 4th-order low pass butterworth filter allowing frequencies below 45 Hz
% 
% See also:
%   eegfiltnew

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Change log
% 20.6.2018 - included low and high pass filter options
% 20.6.2018 - changed to filtering across concatenated trials instead of
% individual trials to try and minimise edge artifacts and speed up script
% 28.5.2024 - NWB changed to include option to apply causal filters instead
% of just acausal filters

function EEG = RELAX_filtbutter(EEG, high, low, ord, type, causal_or_acausal_filter, respect_boundaries)

if nargin < 5
	error('Not enough input arguments.');
end

if exist('causal_or_acausal_filter', 'var')==0
        causal_or_acausal_filter='acausal';
end

if exist('respect_boundaries', 'var')==0
        respect_boundaries=false;
end

%Check inputs
if ~isnumeric(high)
    error('Input ''high'' needs to be a number, not a string. e.g. 1, not ''1''.')
elseif ~isnumeric(low)
    error('Input ''low'' needs to be a number, not a string. e.g. 100, not ''100''.')
elseif high > low
    error('Input ''high'' needs to be less than input ''low''.')
elseif ~isnumeric(ord)
    error('Input ''ord'' needs to be a number, not a string. e.g. 4, not ''4''.')
elseif ~(strcmp(type,'bandpass') || strcmp(type,'bandstop') || strcmp(type,'highpass') || strcmp(type,'lowpass'))
    error('Input ''type'' needs to be either ''bandpass'', ''bandstop'', ''highpass'', or ''lowpass''.')
elseif ~(strcmp(causal_or_acausal_filter,'causal') || strcmp(causal_or_acausal_filter,'acausal'))
    error('Input ''causal_or_acausal_filter'' needs to be either ''causal'' or ''acausal'' '.')
end

if strcmp(type,'highpass') && ~isempty(low)
    warning('Performing high pass filter - ignoring low pass input value.');
end
if strcmp(type,'lowpass') && ~isempty(high)
    warning('Performing low pass filter - ignoring high pass input value.');
end

Fs = EEG.srate;
ordIn = ord/2;

if strcmpi(type,'bandstop')
    type = 'stop';
end
               
if strcmp(type,'bandpass') || strcmp(type,'stop')
    [z1 p1] = butter(ordIn, [high low]./(Fs/2), type);
elseif strcmp(type,'highpass')
    [z1 p1] = butter(ordIn, high./(Fs/2), 'high');
elseif strcmp(type,'lowpass')
    [z1 p1] = butter(ordIn, low./(Fs/2), 'low');
end

data = double(EEG.data);

if respect_boundaries && size(data,3)==1
    % Filter each stretch of data between boundary events separately, with odd
    % reflection padding at each segment edge. The padding width is set from the
    % frequency that governs the filter's transient length: the low edge for a high
    % pass or band pass, the cutoff for a low pass, and the stop-band width for a
    % band stop (a notch rings at its bandwidth, not at its centre frequency).
    segments = RELAX_get_segments(EEG);

    if strcmp(type,'lowpass')
        f_governing = low;
    elseif strcmp(type,'stop')
        f_governing = abs(low - high);
    else % 'bandpass' or 'highpass'
        f_governing = high;
    end
    if isempty(f_governing) || f_governing <= 0
        f_governing = Fs/2;
    end
    padSamples = round(3*Fs/f_governing);

    if strcmpi(causal_or_acausal_filter,'causal')
        % NOTE: reflect-padding a causal filter breaks its causality guarantee, because
        % the reflected samples appended after the segment end are mirrored copies of
        % samples from BEFORE that point - so filter() output near a segment's start
        % would depend on samples that occur later in time. causal_or_acausal_filter is
        % 'acausal' throughout this codebase today (see RELAX_SET_PARAMETERS_AND_RUN*),
        % so this is latent, but causal filtering exists specifically to prevent
        % post-stimulus activity from leaking backwards - reflect-padding it would
        % reintroduce exactly that leak. If this branch is ever exercised with
        % respect_boundaries=true, it needs one-sided (not reflection) padding instead.
        filterFn = @(d) filter(z1,p1,d,[],2);
    else
        filterFn = @(d) filtfilt(z1,p1,d')';
    end

    % minSegmentLength=16: RELAX_apply_by_segments defaults to applying fn to every
    % segment, however short - fine for a generic/pointwise fn, but filtfilt is not
    % generic. RELAX_reflect_pad clamps the pad to segLen-1, so a segment shorter than
    % the requested pad gets a padded length of 3*segLen-2 (worst case). filtfilt on an
    % order-4 filter needs a signal longer than 3*(order)=12 samples, so segLen>=5 is
    % the true floor (3*5-2=13>12); 16 gives >3x that margin without being arbitrary.
    temp = RELAX_apply_by_segments(data, segments, filterFn, padSamples, 16);
    fprintf('Filtering applied separately to %d segment(s) between boundary events (pad = %d samples).\n', ...
        size(segments,1), padSamples);
else
    temp = NaN(size(data,1),EEG.pnts,size(data,3));
    for x = 1:size(data,1) 
        dataIn = data(x,:);
        if strcmpi(causal_or_acausal_filter,'causal')
            dataFilt1 = filter(z1,p1,dataIn);
        end
        if strcmpi(causal_or_acausal_filter,'acausal')
            dataFilt1 = filtfilt(z1,p1,dataIn);
        end
        temp(x,:,:) = reshape(dataFilt1,1,EEG.pnts,[]);
    end 
end

EEG.data = temp;

%display message
if strcmp(type,'bandpass') || strcmp(type,'stop')
    fprintf('Data filtered using a %s zero-phase %s butterworth filter (order = %d) between %d and %d Hz.\n',type,causal_or_acausal_filter,ord,high,low);
elseif strcmp(type,'highpass')
    fprintf('Data filtered using a high pass zero-phase %s butterworth filter (order = %d) below %d Hz.\n',causal_or_acausal_filter,ord,high);
elseif strcmp(type,'lowpass')
    fprintf('Data filtered using a low pass zero-phase %s butterworth filter (order = %d) above %d Hz.\n',causal_or_acausal_filter,ord,low);
end
    

end
