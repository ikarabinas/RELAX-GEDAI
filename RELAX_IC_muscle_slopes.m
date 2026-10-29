% RELAX_IC_muscle_slopes() - Log-frequency log-power slope of each ICA component over
%                           12-70 Hz (excluding +/-6 Hz around line noise), from one
%                           whole-recording Welch spectrum per component. Slopes at or
%                           above RELAX_cfg.MuscleSlopeThreshold are what targeted wICA
%                           treats as muscle.
%
% Usage:
%   >> slopes = RELAX_IC_muscle_slopes(icaact, srate, lineNoiseFrequency);
%
% Inputs:
%   icaact             - components x samples (x trials) activations
%   srate              - sampling rate in Hz
%   lineNoiseFrequency - RELAX_cfg.LineNoiseFrequency
%
% Outputs:
%   slopes - 1 x components

function slopes = RELAX_IC_muscle_slopes(icaact, srate, lineNoiseFrequency)

    options.muscleFreqIn=[12,70];
    options.Freq_to_compute = [1,100];

    eegData = reshape(icaact,size(icaact,1),[]);
    [pxx,fp] = pwelch(eegData',size(eegData,2),[],size(eegData,2),srate);
    FFTout = pxx';
    fp = fp';

    % Calculate FFT bins
    freq=options.Freq_to_compute(1,1):0.5:options.Freq_to_compute(1,2);
    fftBins = zeros(size(FFTout,1),size(freq,2)); %preallocate
    for a=1:size(freq,2)
        [~, index1]=min(abs(fp-((freq(1,a)-0.25))));
        [~, index2]=min(abs(fp-((freq(1,a)+0.25))));
        fftBins(:,a)=mean(FFTout(:,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
    end

    comps=size(icaact,1);
    options.muscleFreqEx=[lineNoiseFrequency-6 lineNoiseFrequency+6];
    slopes = zeros(1, comps);
    for compNum =1:comps
        % Define frequencies to include in the analysis
        [~,fin1] = min(abs(options.muscleFreqIn(1) - freq));
        [~,fin2] = min(abs(options.muscleFreqIn(2) - freq));
        freqHz = freq(1,fin1:fin2);
        freqPow = fftBins(compNum,fin1:fin2);
        % Define frequencies to exclude from fit
        [~,fex1] = min(abs(options.muscleFreqEx(1) - freqHz));
        [~,fex2] = min(abs(options.muscleFreqEx(2) - freqHz));
        freqHz(fex1:fex2) = [];
        freqPow(fex1:fex2) = [];
        % Fit linear regression to log-log data
        p = polyfit(log(freqHz),log(freqPow),1);
        slopes(compNum) = p(1);
    end
end
