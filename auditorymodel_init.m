%% Example of auditory model initialisation
% Leo Varnet - 2018
    
close all
clear cfg
fs = 44100;

%% display settings

cfg.verbose = 'yes'; % details of the computation with process time for each step
cfg.display_step = 'yes'; % automatically plot the results for each step
cfg.display_out = 'no'; % plot the model output
cfg.channels2plot = 1:5;
 
%% generate stim
% create the input within the model function

cfg.stim_gen = 'yes';
cfg.stim_fc = 5000;
cfg.stim_fm = 10;
cfg.stim_Mdepth = 0.5;
cfg.stim_A = 1;
cfg.stim_fun=@create_AM;

%% filterbank
% the options DRNL_filterbank and gtone_O_filterbank require the
% corresponding toolboxes. gtone_filterbank is my implementation of
% Hohmann (2002)'s gammatone definition

cfg.gtone_filterbank = 'yes';
cfg.gtone_fmin=ERB2f(f2ERB(cfg.stim_fc)-2);
cfg.gtone_fmax=ERB2f(f2ERB(cfg.stim_fc)+2)+1;

%% compression
% WARNING: this module requires a prliminary calibration (set kneepoint
% level at 40 dB SPL)

cfg.compression_power = 'yes'; 
cfg.compression_n = [1 1 0.3 1 1];
cfg.compression_knee = 1;

%% hair-cell transduction

cfg.HC_trans = 'yes';

%% adaptation

cfg.adapt_HP = 'yes';

%%
cfg.phase_insens_hilbert = 'yes';

%% modulation filterbank

cfg.mod_filterbank = 'yes';
cfg.modbank_fmin = 2;
cfg.modbank_fmax = 100;
cfg.modbank_Nmod = 10;
cfg.LP_filter = 'yes';

%% downsampling

cfg.downsampling = 'yes';
cfg.downsampling_factor = 10;

%% internal noise
% several internal noise types are implemented, I recommend using additive
% noise only though. Internal noise level is expressed in arbitrary units -
% needs to be calibrated.

cfg.intnoise = 'yes';
cfg.intnoise_addstd = 1;

%% detection
% WARNING : In development. All options are not fully available. Please use
% detect_PEMOxcorr_limit (detection by template matching)

cfg.detect_PEMOxcorr = 'yes';
%cfg.detect_PEMOxcorr_limit = ((1/cfg.stim_fm)*(fs/cfg.downsampling_factor))/2;

%% Generating the template (non-noisy representation of the target)

cfg_temp = cfg;
cfg_temp.intnoise = 'no';
[ template, step ] = auditorymodel([], fs, cfg_temp);

%% Generating internal representations for the target and non-target intervals

[ target, step ] = auditorymodel([], fs, cfg);

cfg.stim_Mdepth=0;
[ nontarget, step ] = auditorymodel([], fs, cfg);

%% Detecting the target

[ interval, ir ] = auditorymodel_detection( {target, nontarget}, template, cfg );

answer = {'correct', 'incorrect'};
fprintf(['The interval selected is ' num2str(interval) ' (' answer{interval} ' answer)\n']);