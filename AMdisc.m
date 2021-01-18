%% Example of AM-depth discrimination experiment
% Leo Varnet - 2020
    
close all
clear cfg
fs = 44100;

%% display settings

cfg.verbose = 'no'; % details of the computation with process time for each step
cfg.display_step = 'no'; % automatically plot the results for each step
cfg.display_out = 'no'; % plot the model output
cfg.channels2plot = 1:5;
 
%% generate stim
% create the input within the model function

cfg.stim_gen = 'no';
cfg.stim_fc = 4000;
cfg.stim_fm = 10;

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

cfg.compression_power = 'no'; 
% cfg.compression_n = [1 1 0.3 1 1];
% cfg.compression_knee = 1;

%% hair-cell transduction

cfg.HC_trans = 'yes';

%% adaptation

cfg.adapt_HP = 'yes';

%%
cfg.phase_insens_hilbert = 'no';

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
cfg.intnoise_addstd = 0.01;

%% detection

%cfg.detect_PEMOxcorr = 'yes';
%cfg.detect_PEMOxcorr_limit = ((1/cfg.stim_fm)*(fs/cfg.downsampling_factor))/2;
cfg.detect_maxlag = 10000;
cfg.detect_threshold = 0;%0.95;%0.1;
cfg.detect_weberf = 1;%1.005;

%% experiment

clear m_thres

standard_mdB = -30:2:-4;

for i_cond=1:length(standard_mdB)
     figure
     for i_rep = 1:100
         m_inc = 1;
         cfg.stim_Mdepth = 10.^((standard_mdB(i_cond))/20)*sqrt(1+m_inc);%m_inc;%
         cfg.intnoise = 'no';
         [ Stemplate, fs] = create_AM(cfg.stim_fc, cfg.stim_fm, cfg.stim_Mdepth, 0.5, fs); Stemplate = Stemplate/rms(Stemplate);
         [ template ] = auditorymodel(Stemplate, fs, cfg); template = template/rms(template(:));
         for i_trial = 1:50
             cfg.intnoise = 'yes';
             fprintf([num2str(m_inc) '\n'])
             cfg.stim_Mdepth = 10.^(standard_mdB(i_cond)/20);%0;%
             [ Sreference, fs] = create_AM(cfg.stim_fc, cfg.stim_fm, cfg.stim_Mdepth, 0.5, fs);  Sreference = Sreference/rms(Sreference);
             fprintf([num2str(cfg.stim_Mdepth) '\n'])
             [ reference ] = auditorymodel(Sreference, fs, cfg);
             cfg.stim_Mdepth = 10.^((standard_mdB(i_cond))/20)*sqrt(1+m_inc);%m_inc;%
             [ Starget, fs] = create_AM(cfg.stim_fc, cfg.stim_fm, cfg.stim_Mdepth, 0.5, fs); Starget = Starget/rms(Starget);
             fprintf([num2str(cfg.stim_Mdepth) '\n'])
             [ target ] = auditorymodel(Starget, fs, cfg);
             
             cfg.display_out = 'yes';
             [ interval, ir ] = auditorymodel_TMdetect({target, reference}, template, cfg );
             cfg.display_out = 'no';
             
             if interval == 1
                 m_inc = m_inc/2;
             else
                 m_inc = m_inc*2;
             end
             if max(ir)<cfg.detect_threshold
                 semilogy(i_trial,m_inc,'ro');break
             elseif (ir(1)/ir(2))<cfg.detect_weberf & (ir(2)/ir(1))<cfg.detect_weberf
                 semilogy(i_trial,m_inc,'bo');break
             else
                 semilogy(i_trial,m_inc,'go');
             end
             drawnow;hold on;
         end
         m_thres(i_cond,i_rep) = m_inc;
     end
end