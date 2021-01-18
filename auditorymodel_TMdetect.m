function [ response, ir, cfg ] = auditorymodel_TMdetect( in, template, cfg )
%AUDITORYMODEL_TMDETECT PEMO-like detection device for the auditory model
%   use as : [ response, ir, cfg ] = auditorymodel_detection( in, template, cfg )
% inputs: 
%   in: structure of n internal representations to be analyzed (matrices of
%   size M x N x L) 
%   template: (matrix of size M x N x L)
%   cfg: config sructure with field cfg.detect_PEMOcorrcoef
%                                or cfg.detect_PEMOxcorr 
%                                or cfg.detect_EPSM. 
%   Additional fields:
%   cfg.detect_PEMOxcorr_limit (default : no limit) and
%   cfg.detect_EPSM_threshold (default = 1 dB)
% output:
%   response: number between 1 and n indicating which response has been
%   selected by the model (interval or template number)
%   ir: n X 1 vector of internal response
%   cfg: config structure as input with all unspecified parameters set to
% default.
%
% Leo Varnet 2016 - last update : 2020

%% Model structure / load arguments

global lagxcorr

if ~iscell(template) & ~iscell(in) & isfield(cfg, 'detect_threshold')
    detect_mode = 1; % detect mode: threshold
    Nalter = 2;
    Ninterval = 1;
    Nsamples = size(in,1);
    Nchan = size(in,2);
    Nmod = size(in,3);
elseif iscell(template) & ~iscell(in)
    detect_mode = 2; % detect mode: single interval (yes/no)
    Nalter = length(template);
    Ninterval = 1;
    Nsamples = size(in,1);
    Nchan = size(in,2);
    Nmod = size(in,3);
elseif iscell(in)
    detect_mode = 3; % detect mode: several intervals (Forced choice)
    Nalter = length(in);
    Ninterval = Nalter;
    Nsamples = size(template,1);
    Nchan = size(template,2);
    Nmod = size(template,3);
else
    error('there must be either (1) several templates and a single interval, (2) a single template and several intervals, or (3) a single template, a single interval, and a threshold value')
end

cfg = set_default_cfg(cfg, 'detect_maxlag', 0, 'detect_norm', 'no', 'detect_threshold', 0, 'detect_weberf', 1);

if cfg.detect_maxlag>Nsamples
    cfg.detect_maxlag = Nsamples;
end

if isyes(cfg.verbose)
    tic
    if detect_mode == 1
        fprintf('detection (threshold)\n')
    elseif detect_mode == 2
        fprintf('detection (single interval - yes/no task)\n')
    elseif detect_mode == 3
        fprintf('detection (several intervals - Forced choice)\n')
    end
    display_cfg(cfg,'detect_maxlag','detect_norm','detect_threshold','detect_weberf');
end

%% Normalization of signals and template

if 1%isyes(cfg.detect_norm)
    if iscell(template)
        for ialter = 1:Nalter
            template{ialter} = template{ialter}/sqrt(cov(template{ialter}(:)));
        end
    else
        template = template/sqrt(cov(template(:),1));
    end
end

if isyes(cfg.detect_norm)
    if iscell(in)
        for ialter = 1:Nalter
            in{ialter} = in{ialter}/sqrt(cov(in{ialter}(:)));
        end
    else
        in = in/sqrt(cov(in(:),1));
    end
end

%% Template-matching

switch detect_mode
    case 1 % yes/no with 1 template (threshold comparison)
        for i_chan = 1:Nchan
            for i_mod = 1:Nmod
                if cfg.detect_maxlag == 0
                    Rmatrix(:,:,i_chan,i_mod) = cov([in(:,i_chan,i_mod) template(:,i_chan,i_mod)]);
                else
                    [Rmatrix(:,:,i_chan,i_mod) lags] = xcov([in(:,i_chan,i_mod) template(:,i_chan,i_mod)],cfg.detect_maxlag, 'biased'); % 'unbiased' so that large lags are not penalized
                end
            end
        end
        
        Rmatrix(isnan(Rmatrix)) = 0;
        Rmatrix = sum(sum(Rmatrix,4),3);
        
        if cfg.detect_maxlag == 0
            ir = Rmatrix(end,1:end-1);
            [~,response] = max(ir);
        else
            xcorr_ir = Rmatrix(:,Nalter);%Rmatrix(:,(1:Nalter)*(Nalter+1));
            
            %[ir, max_xcorr] = max(xcorr_ir);
            
             [ir, max_xcorr] = max(xcorr_ir);
             lagxcorr = [lagxcorr lags(max_xcorr)];
            [~,response] = max(ir);
        end
       
    case 2 % yes/no with 2 templates
        for i_chan = 1:Nchan
            for i_mod = 1:Nmod
                if cfg.detect_maxlag == 0
                    Rmatrix(:,:,i_chan,i_mod) = cov([in(:,i_chan,i_mod) cell2mat(cellfun(@(x)(x(:,i_chan,i_mod)),template,'UniformOutput',0))]);
                else
                    [Rmatrix(:,:,i_chan,i_mod) lags] = xcov([in(:,i_chan,i_mod) cell2mat(cellfun(@(x)(x(:,i_chan,i_mod)),template,'UniformOutput',0))],cfg.detect_maxlag,'unbiased');
                end
            end
        end
        
        Rmatrix(isnan(Rmatrix)) = 0;
        Rmatrix = sum(sum(Rmatrix,4),3);
        
        if cfg.detect_maxlag == 0
            ir = Rmatrix(1,2:end);%Rmatrix(end,1:end-1);
            [~,response] = max(ir);
        else
            xcorr_ir = Rmatrix(:,(1:Nalter)*(Nalter+1));
            ir = max(xcorr_ir);
            [~,response] = max(ir);
        end
    case 3 % Forced choice
        for i_chan = 1:Nchan
            for i_mod = 1:Nmod
                if cfg.detect_maxlag == 0
                    Rmatrix(:,:,i_chan,i_mod) = cov([cell2mat(cellfun(@(x)(x(:,i_chan,i_mod)),in,'UniformOutput',0)) template(:,i_chan,i_mod)]);
                else
                    [Rmatrix(:,:,i_chan,i_mod) lags] = xcov([cell2mat(cellfun(@(x)(x(:,i_chan,i_mod)),in,'UniformOutput',0)) template(:,i_chan,i_mod)],cfg.detect_maxlag,'unbiased');
                end
            end
        end
        
        Rmatrix(isnan(Rmatrix)) = 0;
        Rmatrix = sum(sum(Rmatrix,4),3);
        
        if cfg.detect_maxlag == 0
            ir = Rmatrix(end,1:end-1);
            %[~,response] = max(ir);
        else
            xcorr_ir = Rmatrix(:,(1:Nalter)*(Nalter+1));
            ir = max(xcorr_ir);
            %[~,response] = max(ir);
        end
end

%% Decision

switch detect_mode
    case 1
        %%ir = Rmatrix(1,2:end); %%%TODO
        if ir >= cfg.detect_threshold
            response = 2;
        else
            response = 1; %randi(2)-1;
        end
        
    otherwise
%         if max(ir) < cfg.detect_threshold %%% DOES IT WORK?
%             response = randi(Nalter);
%         else
            if ir(1)/ir(2) >= cfg.detect_weberf
                response = 1;
            elseif ir(2)/ir(1) >= cfg.detect_weberf
                response = 2;
            else
                response = randi(Nalter);
            end
%        end
end

%% Plot

if isyes(cfg.display_out) | isyes(cfg.display_step)
    switch detect_mode
        case 1
            figure; bar(ir); xlabel('interval'); ylabel('matching value'); title('xcorr with template')
        case 2
            
        case 3
            if cfg.detect_maxlag == 0
                figure; bar(ir); xlabel('interval'); ylabel('matching value'); title('xcorr with template')
            else
                figure; plot(lags,xcorr_ir); xlabel('lags (sample)'); ylabel('matching value'); legend({'interval 1','interval 2'});title('xcorr with template')
            end
    end
end

