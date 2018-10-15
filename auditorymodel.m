function [ outsig, step, cfg ] = auditorymodel( varargin )
% AUDITORYMODEL 
% Use as [ outsig, steps, cfg ] = auditorymodel([],fs,cfg) 
% or [ outsig, steps, cfg ] = auditorymodel(S, fs, cfg)
% Processes input signal S (sampling rate fs) through the model
% defined in structure cfg
%
% Outputs:
% - outsig is a Ntime X Nchan X Nmod matrix with Ntime number of samples,
%  Nchan number of frequency channels, Nmod number of modulation channels.
% - step contains the temporary representations of the stimulus at the
% different stages of the model.
% - cfg config structure as input with all unspecified parameters set to
% default.
%
% List of cfg parameters, step-by-step:
% - Stimulus creation: cfg.stim_gen (in case S is empty)
%  cfg.stim_fc, cfg.stim_fs, cfg.stim_fm, cfg.stim_Mdepth,
%  cfg.stim_duration, cfg.stim_A, cfg.stim_carrier_waveform,
%  cfg.stim_M_waveform, cfg.stim_fun  
% - Peripheral filtering: cfg.gtone_filterbank, cfg.gtone_O_filterbank or
% cfg.DRNL_filterbank 
%  cfg.gtone_fmin, cfg.gtone_fmax (only for gtone_filterbank and
%  gtone_O_filterbank)
%  gtone_order, cfg.gtone_ERB (only for gtone_O_filterbank)
% - Hilbert envelope extraction: cfg.envelope_extract
% - Compression: cfg.compression_power or compression_brokenstick
%  cfg.compression_n, cfg.compression_knee (only for brokenstick),
%  cfg.compression_smooth (only for brokenstick) 
% - Hair-cell transduction: cfg.HC_trans
%  cfg.HC_fcut
% - Adaptation: cfg.adapt_FBloops or cfg.adapt_HP
%  cfg.adapt_FBloops_timeconst (only for FBloops), cfg.adapt_HP_fc and
%  cfg.adapt_HP_order (only for HP)
% - Modulation filterbank: cfg.mod_filterbank or cfg.modS_filterbank
%  cfg.modbank_fmin, cfg.modbank_fmax, cfg.LP_filter
%  cfg.modbank_Qfactor (only for mod_filterbank)
% - Phase insensitivity: cfg.phase_insens_hilbert or cfg.phase_insens_filt
%  cfg.phase_insens_cut, cfg.phase_insens_order
% - Downsampling: cfg.downsampling
%  cfg.downsampling_factor
% - Memory decay: cfg.memorydecay
%  cfg.memorydecay_tau
% - Internal noise: cfg.intnoise
%  cfg.intnoise_addstd, cfg.intnoise_multratio, intnoise_memstd
% - General display options: cfg.verbose, cfg.display_step, cfg.display_out
%  cfg.channels2plot
%
% Leo Varnet 2016 - last update : january 2018

if nargin == 1
    cfg = varargin{1};
elseif nargin ==3
    input = varargin{1};
    fs = varargin{2};
    cfg = varargin{3};
end
step = {};

%% Model structure / load arguments

% display options
cfg = set_default_cfg(cfg, 'verbose', 'yes', 'display_step', 'no', 'display_out', 'no', 'display_Freq', [0 8000], 'display_Overlap', 0.5, 'display_Nwindow', 512, 'display_NFFT', 1024);
        
% Creating the stimulus
if isfield(cfg, 'stim_gen')
    if isyes(cfg.stim_gen)
        cfg = set_default_cfg(cfg, 'stim_fc', 1000, 'stim_fs', 44100, 'stim_fm', 20, 'stim_Mdepth', 0.5, 'stim_duration', 1, 'stim_A', 1, 'stim_carrier_waveform', 'sin', 'stim_M_waveform', 'sin', 'stim_fun', @create_AM);
        fs = cfg.stim_fs;
    end
else
    cfg.stim_gen = 'no';
end

% Simple Hilbert envelope extraction
if isfield(cfg, 'envelope_extract')
    
else
    cfg.envelope_extract = 'no';
end

% Gammatone filterbank
if isfield(cfg, 'gtone_filterbank')
    if isyes(cfg.gtone_filterbank)
        cfg = set_default_cfg(cfg, 'gtone_fmin', 70, 'gtone_fmax', 6700);
    end
else
    cfg.gtone_filterbank = 'no';
end

% Oldenburg's Gammatone filterbank
if isfield(cfg, 'gtone_O_filterbank')
    if isyes(cfg.gtone_O_filterbank)
        cfg = set_default_cfg(cfg, 'gtone_fmin', 70, 'gtone_fmax', 6700, 'gtone_order', 4, 'gtone_ERB', 1);
    end
else
    cfg.gtone_O_filterbank = 'no';
end

% DRNL filterbank
if isfield(cfg, 'DRNL_filterbank')
    if isyes(cfg.DRNL_filterbank)
        cfg = set_default_cfg(cfg, 'drnl_fmin', 70, 'drnl_fmax', 6700);
    end
else
    cfg.DRNL_filterbank = 'no';
end

if isyes(cfg.gtone_filterbank)+isyes(cfg.gtone_O_filterbank)+isyes(cfg.DRNL_filterbank)>1
    error('too many peripheral filtering stages. You have to choose between gtone_filterbank, gtone_O_filterbank or DRNL_filterbank');
end

% "Power" compression
if isfield(cfg, 'compression_power')
    if isyes(cfg.compression_power)
        cfg = set_default_cfg(cfg, 'compression_n', 0.6);
    end
else
    cfg.compression_power = 'no';
end

% "Brokenstick" compression
if isfield(cfg, 'compression_brokenstick')
    if isyes(cfg.compression_brokenstick)
        cfg = set_default_cfg(cfg, 'compression_n', 0.25);%, 'compression_a', 30274, 'compression_b', 0.06);
    end
else
    cfg.compression_brokenstick = 'no';
end

% Stefan's "Brokenstick" compression
if isfield(cfg, 'compression_sbrokenstick')
    if isyes(cfg.compression_sbrokenstick)
        cfg = set_default_cfg(cfg, 'compression_n', 0.6, 'compression_knee', 1, 'compression_smooth', 0.5);
    end
else
    cfg.compression_sbrokenstick = 'no';
end

if isyes(cfg.compression_power)+isyes(cfg.compression_sbrokenstick)>1
    error('too many compression stages. You have to choose between compression_power or compression_brokenstick');
end

% Hair-cell transduction
if isfield(cfg, 'HC_trans')
    if isyes(cfg.HC_trans)
        cfg = set_default_cfg(cfg, 'HC_fcut', 1500);
    end
else
    cfg.HC_trans = 'no';
end

% Feedback loops
if isfield(cfg, 'adapt_FBloops')
    if isyes(cfg.adapt_FBloops)
        cfg = set_default_cfg(cfg, 'adapt_FBloops_timeconst', [0.005 0.05 0.129 0.253 0.5]);
    end
else
    cfg.adapt_FBloops = 'no';
end

% Adaptation by high-pass filtering
if isfield(cfg, 'adapt_HP')
    if isyes(cfg.adapt_HP)
        cfg = set_default_cfg(cfg, 'adapt_HP_fc', 3, 'adapt_HP_order', 1);
    end
else
    cfg.adapt_HP = 'no';
end

if isyes(cfg.adapt_FBloops)+isyes(cfg.adapt_HP)>1
    error('too many adaptation stages. You have to choose between adapt_FBloops or adapt_FBloops');
end

% Modulation filterbank
if isfield(cfg, 'mod_filterbank')
    if isyes(cfg.mod_filterbank)
        cfg = set_default_cfg(cfg, 'modbank_fmin', 2, 'modbank_fmax', 150, 'modbank_LPfilter', 'no', 'modbank_Nmod', [], 'modbank_Qfactor', 1);
    end
else
    cfg.mod_filterbank = 'no';
end 

% Stefan's modulation filterbank
if isfield(cfg, 'modS_filterbank')
    if isyes(cfg.modS_filterbank)
        cfg = set_default_cfg(cfg, 'modbank_fmin', 2, 'modbank_fmax', 120, 'modbank_LPfilter', 'no', 'modbank_Nmod', []);
    end
else
    cfg.modS_filterbank = 'no';
end 

if isyes(cfg.mod_filterbank)+isyes(cfg.modS_filterbank)>1
    error('too many modulation filtering stages. You have to choose between mod_filterbank or modS_filterbank');
end

% Phase insensitivity (hilbert)
if isfield(cfg, 'phase_insens_hilbert')
    if isyes(cfg.phase_insens_hilbert)
        cfg = set_default_cfg(cfg, 'phase_insens_cut',10); 
    end
else
    cfg.phase_insens_hilbert = 'no';
end

% Phase insensitivity (filtering)
if isfield(cfg, 'phase_insens_filt')
    if isyes(cfg.phase_insens_filt)
        cfg = set_default_cfg(cfg, 'phase_insens_cut',6, 'phase_insens_order', 2); 
    end
else
    cfg.phase_insens_filt = 'no';
end

if isyes(cfg.phase_insens_hilbert)+isyes(cfg.phase_insens_filt)>1
    error('too many phase insensitivity stages. You have to choose between phase_insens_hilbert or phase_insens_filt');
end

% Downsampling
if isfield(cfg, 'downsampling')
    if isyes(cfg.downsampling)
        cfg = set_default_cfg(cfg, 'downsampling_factor',10); 
    end
else
    cfg.downsampling = 'no';
end

% Memory decay
if isfield(cfg, 'memorydecay')
    if isyes(cfg.memorydecay)
        cfg = set_default_cfg(cfg, 'memorydecay_tau', 1.2);
    end
else
    cfg.memorydecay = 'no';
end 

% Internal noises
if isfield(cfg, 'intnoise')
    if isyes(cfg.intnoise)
        cfg = set_default_cfg(cfg, 'intnoise_addstd', 0, 'intnoise_multratio', 0, 'intnoise_memstd', 0, 'memorydecay_tau', 0);
    end
else
    cfg.intnoise = 'no';
end 

% general options
cfg = set_default_cfg(cfg, 'verbose', 'yes', 'display_step', 'no', 'display_out', 'yes');
if isyes(cfg.display_step) || isyes(cfg.display_out)
    cfg = set_default_cfg(cfg, 'channels2plot', 1);
end

if isyes(cfg.verbose)
    fprintf('\nStarting the model\n')
end

%% Creating the stimulus

if isyes(cfg.stim_gen)
    if isyes(cfg.verbose)
        tic
        fprintf('Creating the stimulus\n')
        display_cfg(cfg,'stim_fc','stim_fs','stim_fm','stim_Mdepth','stim_duration','stim_A','stim_carrier_waveform','stim_M_waveform','stim_fun');
    end
    
    if isequal(cfg.stim_fun, @create_AM) || isequal(cfg.stim_fun, @create_FM)
        [input,fs] = cfg.stim_fun(cfg.stim_fc, cfg.stim_fm, cfg.stim_Mdepth, cfg.stim_duration, fs, cfg.stim_A, cfg.stim_carrier_waveform, cfg.stim_M_waveform);
    elseif isequal(cfg.stim_fun, @create_MM)
        [input,fs] = cfg.stim_fun(cfg.stim_fc, cfg.stim_fm, cfg.stim_Mdepth, cfg.stim_fm, cfg.stim_Mdepth, 0);
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
end

%% Prepare input

input = input(:);

%input = setdbspl(input,60);

if isyes(cfg.display_step)
    plot_sound(input, fs, cfg.display_Freq, cfg.display_Overlap, cfg.display_Nwindow, cfg.display_NFFT);
    set(gcf, 'Name', 'input sound')
end

t = (1:length(input))'/fs;
Nsamples = length(input);

if nargout>1
    step.t = t;
    step.Nsamples = Nsamples;
end

%% Gammatone filterbank

if isyes(cfg.gtone_filterbank)
    if isyes(cfg.verbose)
        tic
        fprintf('Gammatone filterbank\n')
        display_cfg(cfg,'gtone_fmin','gtone_fmax')
    end
    
    [gtone_responses, ~, fc, fb] = gammatone_filtering( input, cfg.gtone_fmin, cfg.gtone_fmax, 1, 0, fs);
    Nsamples = length(gtone_responses);
    Nchannels = length(fc);
    gtone_responses = real(gtone_responses);
    outsig = gtone_responses;
    if isyes(cfg.display_step)
        plot_channels(t, gtone_responses, cfg.channels2plot);        
        set(gcf, 'Name', 'Gammatones responses');
    end
    if nargout>1
        step.fc = fc;
        step.fb = fb;
        step.Nchannels = Nchannels;
        step.gtone_responses = gtone_responses;
    end
    
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear gtone_responses
end

%% Oldenburg gammatone filterbank

if isyes(cfg.gtone_O_filterbank)
    if isyes(cfg.verbose)
        tic
        fprintf('Gammatone filterbank\n')
        display_cfg(cfg, 'gtone_fmin', 'gtone_fmax', 'gtone_order', 'gtone_ERB')
    end
    
    analyzer = Gfb_Analyzer_new(fs, cfg.gtone_fmin, cfg.gtone_fmin, cfg.gtone_fmax, 1, cfg.gtone_order, cfg.gtone_ERB);
    fc = analyzer.center_frequencies_hz;
    Nchannels = length(fc);
    [gtone_responses, ~] = Gfb_Analyzer_process(analyzer, input);
        
    gtone_responses = real(gtone_responses)';
    outsig = gtone_responses;
    if isyes(cfg.display_step)
        plot_channels(t, gtone_responses, cfg.channels2plot);        
        set(gcf, 'Name', 'Gammatones responses');
    end
    if nargout>1
        step.fc = fc;
        step.Nchannels = Nchannels;
        step.gtone_responses = gtone_responses;
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear gtone_responses
end

%% DRNL filterbank

if isyes(cfg.DRNL_filterbank)
    if isyes(cfg.verbose)
        tic
        fprintf('DRNL filterbank\n')
        display_cfg(cfg,'drnl_fmin','drnl_fmax')
    end

    [DRNL_responses, fc] = drnl(input, fs, 'flow', cfg.drnl_fmin, 'fhigh', cfg.drnl_fmax, 'nlinonly');%drnl(input, fs); %
    DRNL_responses=DRNL_responses;
    Nchannels = length(fc);
    
    outsig = DRNL_responses;
    
    if nargout>1
        step.DRNL_responses = DRNL_responses;
        step.fc = fc;
        step.Nchannels = Nchannels;
    end
    if isyes(cfg.display_step)
        plot_channels(t, DRNL_responses, cfg.channels2plot);        
        set(gcf, 'Name', 'DRNL responses');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear DRNL_responses
end

%% Simple Hilbert envelope extraction

if isyes(cfg.envelope_extract)
    if isyes(cfg.verbose)
        tic
        fprintf('Hilbert envelope extraction\n')
    end
    if ~exist('output')
        outsig = input;
        fc = 1;
        Nchannels = 1;
    end
    hilbert_envelope = hilbert_extraction(input, fs);
    outsig = hilbert_envelope;
    if isyes(cfg.display_step)
        plot_sound(hilbert_envelope, fs, cfg.display_Freq, cfg.display_Overlap, cfg.display_Nwindow, cfg.display_NFFT);
        set(gcf, 'Name', 'Hilbert envelope');
    end
    if nargout>1
        step.fc = fc;
        step.Nchannels = Nchannels;
        step.hilbert_envelope = hilbert_envelope;
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear hilbert_envelope
end

%% "Power" compression (relative to frequency fa)

if isyes(cfg.compression_power)
    if isyes(cfg.verbose)
        tic
        fprintf('Power compression\n')
        display_cfg(cfg, 'compression_n');
    end
    
    if size(cfg.compression_n)==1
        cfg.compression_n = cfg.compression_n*ones(1,Nchannels);
    elseif size(cfg.compression_n,2) ~= Nchannels
        error('simwork.cfg.compression_n does not have the same number of channels as the output of the Gammatone or DRNL Filterbank');
    end
    
    compressed_response = sign(outsig).*abs(outsig).^(ones(Nsamples,1)*cfg.compression_n);
    outsig = compressed_response;
    
    if nargout>1
        step.compressed_response = compressed_response;
    end
    if isyes(cfg.display_step)
        plot_channels(t, compressed_response, cfg.channels2plot);        
        set(gcf, 'Name', 'Power-compressed responses responses');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear compressed_response
end

%% "Brockenstick" compression (relative to frequency fa)

if isyes(cfg.compression_brokenstick)
    if isyes(cfg.verbose)
        tic
        fprintf('Brokenstick compression\n')
        display_cfg(cfg, 'compression_n', 'compression_a', 'compression_k');
    end
    
    if length(cfg.compression_n)==1
        cfg.compression_n = cfg.compression_n*ones(1,Nchannels);
    elseif size(cfg.compression_n,2) ~= Nchannels
        error('simwork.cfg.compression_n does not have the same number of channels as the output of the Gammatone or DRNL Filterbank');
    end
    cfg.compression_b = cfg.compression_a.*(cfg.compression_k.^(1-cfg.compression_n));
    
    compressed_response = sign(outsig).*min(cfg.compression_a*abs(outsig),(ones(Nsamples,1)*cfg.compression_b).*(abs(outsig).^(ones(Nsamples,1)*cfg.compression_n)));
    outsig = compressed_response;
    
    if nargout>1
        step.compressed_response = compressed_response;
    end
    if isyes(cfg.display_step)
        plot_channels(t, compressed_response, cfg.channels2plot);        
        set(gcf, 'Name', 'Brockenstick-compressed responses');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear compressed_response
end

%% Stefan's "Brokenstick" compression (relative to frequency fa)

if isyes(cfg.compression_sbrokenstick)
    if isyes(cfg.verbose)
        tic
        fprintf('Stefan''s brokenstick compression\n')
        display_cfg(cfg, 'compression_n','compression_knee','compression_smooth');
    end
      
    if size(cfg.compression_n)==1
        cfg.compression_n = cfg.compression_n*ones(1,Nchannels);
    elseif size(cfg.compression_n,2) ~= Nchannels
        error('simwork.cfg.compression_n does not have the same number of channels as the output of the Gammatone or DRNL Filterbank');
    end
    
    for ichan=1:Nchannels
        compressed_response(:,ichan) = nltrans3(outsig(:,ichan), cfg.compression_knee, cfg.compression_n(ichan), cfg.compression_smooth);
    end
    
%     idx_ch_above = find(fc>cfg.BScomp_fa);
%     idx_ch_below = find(fc<=cfg.BScomp_fa);
%     compressed_response = zeros(Nchannels, Nsamples);
%     compressed_response(idx_ch_above,:) = nltrans3(outsig(idx_ch_above,:), cfg.BScomp_knee, cfg.BScomp_n_above, cfg.BScomp_smooth);
% %     compressed_response(idx_ch_above(1),:) = nltrans3(outsig(idx_ch_above(1),:), cfg.BScomp_knee, cfg.BScomp_n_above, cfg.BScomp_smooth);
% %     compressed_response(idx_ch_above(2),:) = nltrans3(outsig(idx_ch_above(2),:), cfg.BScomp_knee, 1, cfg.BScomp_smooth);
%     compressed_response(idx_ch_below,:) = nltrans3(outsig(idx_ch_below,:), cfg.BScomp_knee, cfg.BScomp_n_below, cfg.BScomp_smooth);
    outsig = compressed_response;
    
    if nargout>1
        step.compressed_response = compressed_response;
    end
    if isyes(cfg.display_step)
        plot_channels(t, compressed_response, cfg.channels2plot);        
        set(gcf, 'Name', 'Brokenstick-compressed responses responses');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear compressed_response
end
clear idx_ch_above idx_ch_below

%% Hair-cell transduction (envelope extraction by half-wave rectification + low-pass filter with high cutoff frequency)

if isyes(cfg.HC_trans)
    if isyes(cfg.verbose)
        tic
        fprintf('Hair-cell transduction\n')
        display_cfg(cfg, 'HC_fcut');
    end

    HC_response = envelope_extraction(outsig, fs, cfg.HC_fcut);
    outsig = HC_response;
    if nargout>1
        step.HC = HC_response;
    end
    if isyes(cfg.display_step)
        plot_channels(t, HC_response, cfg.channels2plot);        
        set(gcf, 'Name', 'Hair-cell responses');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear HC
end

%% Expansion

 
%% Feedback loops

if isyes(cfg.adapt_FBloops)
    if isyes(cfg.verbose)
        tic
        fprintf('Feedback loops\n')
        display_cfg(cfg, 'adapt_FBloops_timeconst');
    end
    adapted_response = outsig;
    for ichan = 1:Nchannels
        adapted_response(:,ichan) = adaptloop(adapted_response(:,ichan),fs,10,0, cfg.adapt_FBloops_timeconst);
    end
    outsig = adapted_response;
    if nargout>1
        step.FB = adapted_response;
    end
    if isyes(cfg.display_step)
        plot_channels(t, adapted_response, cfg.channels2plot);
        set(gcf, 'Name', 'Feedback loop responses');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear ichan adapted_response
end

%% Adaptation by high-pass filtering

if isyes(cfg.adapt_HP)
    if isyes(cfg.verbose)
        tic
        fprintf('Adaptation by high-pass filtering\n')
        display_cfg(cfg, 'adapt_HP_fc', 'adapt_HP_order');
    end
    [B,A] = butter(cfg.adapt_HP_order,2*(cfg.adapt_HP_fc/fs),'high');
    %adapted_response=filtfilt(B,A,outsig')';
    adapted_response=filter(B,A,outsig);
    if isyes(cfg.display_step)
        plot_channels(t, adapted_response, cfg.channels2plot);
        set(gcf, 'Name', 'summed HP-adapted responses');
    end
    outsig = adapted_response;
    if nargout>1
    step.adapted_response = adapted_response;
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear B A adapted_response
end

%% Excitation Patterns

if isyes(cfg.display_step) || ((~isyes(cfg.mod_filterbank) && ~isyes(cfg.modS_filterbank)) && isyes(cfg.display_out))

    % parameters
    decision_statistics_fun = @(x) std(x, 0, 1);
    time_EP = 1:length(t);%floor(length(t)/5):length(t);
    
    excitation_pattern = decision_statistics_fun(outsig(time_EP,:));
    if isyes(cfg.display_step)
        figure; plot(fc, excitation_pattern, 'k-o')
        ylabel('Decision statistics (std(E))')
        xlabel('center frequency (Hz)');
        title('excitation pattern');
    end
    clear excitation_pattern time_EP
end

%% Modulation filterbank

if isyes(cfg.mod_filterbank)
    if isyes(cfg.verbose)
        tic
        fprintf('Modulation filterbank\n')
        display_cfg(cfg, 'modbank_fmin', 'modbank_fmax', 'modbank_LPfilter', 'modbank_Nmod', 'modbank_Qfactor');
    end
    
    if isyes(cfg.modbank_LPfilter)
        [B,A] = butter(6,2*(100/fs)); % [B,A] = butter(1,2*(150/fs));%
        outsig = filter(B,A,outsig);
    end
    
    [BB, AA, fmc] = modulation_filterbank(cfg.modbank_fmin, cfg.modbank_fmax, fs, cfg.modbank_Nmod, cfg.modbank_Qfactor);
    Nmodchannels = length(fmc);
    E_mod = apply_filterbank(BB, AA, outsig);
    
%     for i=1:length(fmc)
%         if fmc(i)>10
%             for j=1:length(fc)
% %                 rectif = max((squeeze(E_mod(:,j,i))),0);
% %                 [B,A] = butter(2,2*(6/fs));
% %                 envenv = filtfilt(B,A,rectif);
%                 envenv = hilbert_extraction(squeeze(E_mod(:,j,i)),fs);
%                 E_mod(:,j,i) = (envenv / rms(envenv)) * rms(E_mod(:,j,i));
%             end
%         end
%     end
%     
    outsig = E_mod;
        
    if nargout>1
        step.E_mod = E_mod;
        step.fmc = fmc;
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    if isyes(cfg.display_step)
        plot_channels(t, E_mod, cfg.channels2plot);        
        set(gcf, 'Name', 'Modulation filterbank responses');    
        legend(num2str(fmc'))
    end
    clear E_mod AA BB A B
end

%% Stefan's modulation filterbank

if isyes(cfg.modS_filterbank)
    if isyes(cfg.verbose)
        tic
        fprintf('Modulation filterbank\n')
        display_cfg(cfg, 'modbank_fmin', 'modbank_fmax', 'modbank_LPfilter');
    end
    
    style=2-isyes(cfg.LP_filter);
    [fmc, ~] = mFB(outsig(1,:),cfg.modbank_fmin,cfg.modbank_fmax,2,style,fs);
    Nmodchannels = length(fmc);
    %E_mod = nan(Nsamples, Nmodchannels, Nchannels);
    E_mod = nan(Nchannels,Nmodchannels,Nsamples);
    for ichan = 1:Nchannels
        %[fmc, E_mod(:,:,ichan)] = mFB(outsig(ichan,:),cfg.modbank_fmin,cfg.modbank_fmax,2,style,fs);
        [~, A] = mFB(outsig(ichan,:),cfg.modbank_fmin,cfg.modbank_fmax,2,style,fs);
        E_mod(ichan,:,:) = permute(A,[3,2,1]);
    end
    outsig = E_mod;
    if nargout>1
        step.E_mod = E_mod;
        step.fmc = fmc;
    end
    if isyes(cfg.display_step)
        plot_channels(t, E_mod, cfg.channels2plot);            
        legend(num2str(fmc'))
        set(gcf, 'Name', 'Modulation filterbank responses');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear E_mod style A
end

%% Phase insensitivity (hilbert)

if isyes(cfg.phase_insens_hilbert)
    if isyes(cfg.verbose)
        tic
        fprintf('Phase insensitivity (hilbert)\n')
        display_cfg(cfg, 'phase_insens_cut');
    end
    
    E_phase_ins = outsig;
    
    for i=1:Nmodchannels
        if fmc(i)>cfg.phase_insens_cut
            for j=1:length(fc)
                % hilbert
                envenv = hilbert_extraction(squeeze(outsig(:,j,i)),fs);
                E_phase_ins(:,j,i) = (envenv / rms(envenv)) * rms(outsig(:,j,i));
            end
       end
    end
    
    outsig = E_phase_ins;
    if isyes(cfg.display_step)
        plot_channels(t, E_phase_ins, cfg.channels2plot);            
        legend(num2str(fmc'))
        set(gcf, 'Name', 'Phase insensitive responses');
    end
    if nargout>1
        step.E_phase_ins = E_phase_ins;
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
end

%% Phase insensitivity (filtering)

if isyes(cfg.phase_insens_filt)
    if isyes(cfg.verbose)
        tic
        fprintf('Phase insensitivity (filtering)\n')
        display_cfg(cfg, 'phase_insens_cut', 'phase_insens_order');
    end
    
    E_phase_ins = outsig;
    
    for i=1:Nmodchannels
            for j=1:length(fc)
%                 % filtre
                rectif = max((squeeze(outsig(:,j,i))),0);
                [B,A] = butter(cfg.phase_insens_order,2*(cfg.phase_insens_cut/fs));%butter(1,2*(20/fs));%
                envenv = filtfilt(B,A,rectif);
%                 
                E_phase_ins(:,j,i) = (envenv / rms(envenv)) * rms(outsig(:,j,i));
            end
    end
    
    outsig = E_phase_ins;
    if isyes(cfg.display_step)
        plot_channels(t, E_phase_ins, cfg.channels2plot);            
        legend(num2str(fmc'))
        set(gcf, 'Name', 'Phase insensitive responses');
    end
    if nargout>1
        step.E_phase_ins = E_phase_ins;
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
end

%% Downsampling

if isyes(cfg.downsampling)
    if isyes(cfg.verbose)
        tic
        fprintf('downsampling\n')
        display_cfg(cfg, 'downsampling_factor'); 
    end
    fs=fs/cfg.downsampling_factor;
    step.fs_outsig = fs;
    if (isyes(cfg.mod_filterbank) || isyes(cfg.modS_filterbank))
        outsig = outsig(1:cfg.downsampling_factor:end,:,:);
    else
        outsig = outsig(1:cfg.downsampling_factor:end,:);        
    end    
    t_initial_sampling = t;
    t = (1:length(outsig))'/fs;
    Nsamples = length(t);
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    if nargout>1
        step.t_initial_sampling = t_initial_sampling;
        step.t = t;
        step.Nsamples = Nsamples;
    end
else
    step.fs_outsig = fs;
end

%% Plotting the Modulation Excitation Pattern

if isyes(cfg.display_step) && (isyes(cfg.mod_filterbank) || isyes(cfg.modS_filterbank))
    plot_modep(fc, fmc, outsig);
    plot_channels(t, outsig, cfg.channels2plot);
    set(gcf, 'Name', 'Modulation filterbank responses');
    legend(num2str(fmc'))
end

%% Internal noise

if isyes(cfg.intnoise) & (cfg.intnoise_addstd>0 | cfg.intnoise_multratio>0 | cfg.intnoise_memstd>0)
    if isyes(cfg.verbose)
        tic
        fprintf('Internal noise\n');
        display_cfg(cfg, 'intnoise_addstd', 'intnoise_multratio', 'memorydecay_tau', 'intnoise_memstd');
    end
    rms_Emod = sqrt(mean(outsig.^2,1));
    noisy_Emod = outsig;
    % additive noise
    if cfg.intnoise_addstd>0
        noisy_Emod = noisy_Emod+randn(size(noisy_Emod))*cfg.intnoise_addstd;
    end
    % multiplicative noise
    if cfg.intnoise_multratio>0
        multnoise = randn(size(outsig));
        for i_chan = 1:Nchannels
            for i_mod = 1:Nmodchannels
                multnoise(:,i_chan,i_mod) = cfg.intnoise_multratio*multnoise(:,i_chan,i_mod)*rms_Emod(1,i_chan,i_mod)/sqrt(mean(multnoise(:,i_chan,i_mod).^2,1));
            end
        end
        noisy_Emod = noisy_Emod+multnoise;
    end
    % memory noise
    if cfg.intnoise_memstd>0
        [sizex, sizey, sizez] = size(noisy_Emod);
        [idx_value, ~, ~] = ndgrid(1:sizex,1:sizey,1:sizez);
        noisedecay = exp(t(flip(idx_value,1))/cfg.memorydecay_tau);
        %noise = randn(size(noisy_Emod)).*(exp(-t(idx_value)/cfg.memorydecay_tau));
        noise = randn(size(noisy_Emod)).*noisedecay;
        noisy_Emod = noisy_Emod+noise*cfg.intnoise_memstd;
    end
    outsig = noisy_Emod;
    if nargout>1
        step.noisy_Emod = noisy_Emod;
    end
    if isyes(cfg.display_step)
        plot_channels(t, noisy_Emod, cfg.channels2plot);
        set(gcf, 'Name', 'internal representation with internal noise');
    end
    if isyes(cfg.verbose)
        fprintf(['  elapsed time: ' num2str(toc) '\n'])
    end
    clear noisy_Emod rms_Emod multnoise i_chan i_mod
end

%% Memory decay
% 
% if isyes(cfg.memorydecay)
%     
%     if isyes(cfg.verbose)
%         tic
%         fprintf('memory decay\n')
%         display_cfg(cfg, 'memorydecay_tau');
%     end
%     memory_Emod = outsig;
%     
%     [sizex, sizey, sizez] = size(memory_Emod);
%     [idx_value, ~, ~] = ndgrid(1:sizex,1:sizey,1:sizez);
%     memory_Emod = memory_Emod.*(exp(t(idx_value)/cfg.memorydecay_tau));
%     outsig = memory_Emod;
%     if nargout>1
%         step.memory_Emod = memory_Emod;
%     end
%     if isyes(cfg.display_step)
%         plot_channels(t, memory_Emod, cfg.channels2plot);
%         set(gcf, 'Name', 'internal representation with memory decay');
%     end
%     if isyes(cfg.verbose)
%         fprintf(['  elapsed time: ' num2str(toc) '\n'])
%     end
%     clear memory_Emod sizex sizey sizez idx_value
% end

%% Additive memory noise

% if isyes(cfg.intnoise) & cfg.intnoise_memstd>0
%     if isyes(cfg.verbose)
%         tic
%         fprintf('Memory noise\n');
%         display_cfg(cfg, 'intnoise_memstd');
%     end
%     noisy_Emod = outsig;
%     if cfg.intnoise_memstd>0
%         noisy_Emod = noisy_Emod+randn(size(noisy_Emod))*cfg.intnoise_memstd;
%     end
%     outsig = noisy_Emod;
%     if nargout>1
%         step.noisy_Emod = noisy_Emod;
%     end
%     if isyes(cfg.display_step)
%         plot_channels(t, noisy_Emod, cfg.channels2plot);
%         set(gcf, 'Name', 'internal representation with memory noise');
%     end
%     if isyes(cfg.verbose)
%         fprintf(['  elapsed time: ' num2str(toc) '\n'])
%     end
%     clear noisy_Emod
% end

%% Plotting the Modulation Excitation Pattern

if (isyes(cfg.display_step) || isyes(cfg.display_out)) && (isyes(cfg.mod_filterbank) || isyes(cfg.modS_filterbank))
    plot_modep(fc, fmc, outsig);
    title('Final internal representation');
    plot_channels(t, outsig, cfg.channels2plot, @plot, 1,5);
    set(gcf, 'Name', 'Final internal representation');
    legend(num2str(fmc'))
end

%% Suppr onset/offset

%outsig = outsig(0.05*fs:end-0.02*fs,:,:);

end

