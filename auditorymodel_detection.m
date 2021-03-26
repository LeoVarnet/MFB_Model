function [ interval, ir, cfg ] = auditorymodel_detection( in, template, cfg )
%AUDITORYMODEL_DETECTION Detection device for the auditory model
%   use as : [ interval, ir, cfg ] = auditorymodel_detection( in, template, cfg )
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
%   interval: number between 1 and n indicating which interval has been
%   selected by the model 
%   ir: n X 1 vector of internal response
%   cfg: config structure as input with all unspecified parameters set to
% default.
%
% Leo Varnet 2016 - last update : 2020

%% Model structure / load arguments

Nalter = length(in);
Nsamples = size(template,1);
Nchan = size(template,2);
Nmod = size(template,3);

if isfield(cfg, 'detect_PEMOcorrcoef')
    
else
    cfg.detect_PEMOcorrcoef = 'no';
end

if isfield(cfg, 'detect_PEMOxcorr')
    if isyes(cfg.detect_PEMOxcorr)
        cfg = set_default_cfg(cfg, 'detect_PEMOxcorr_limit', Nsamples);
    end
else
    cfg.detect_PEMOxcorr = 'no';
end

if isfield(cfg, 'detect_EPSM')
    if isyes(cfg.detect_EPSM)
        cfg = set_default_cfg(cfg, 'detect_EPSM_threshold', 0);
    end
else
    cfg.detect_EPSM = 'no';
end

if isyes(cfg.detect_PEMOcorrcoef)+isyes(cfg.detect_PEMOxcorr)+isyes(cfg.detect_EPSM)>1
    error('too many decision devices selected. You have to choose between detect_PEMOcorrcoef, detect_PEMOxcorr or detect_EPSM');
end

%% Decision device

if isfield(cfg,'detect_PEMOcorrcoef') && isyes(cfg.detect_PEMOcorrcoef)
    for ialter = 1:Nalter
         cctemp = corrcoef(in{ialter}(:),template(:));
         %cctemp = cov(in{ialter}(:),template(:));
         ir(ialter)=cctemp(1,2);
    end
   
%     % xcorr of internal representation and template
%     cc = sum(in.*repmat(template,1,size(in,2)));
%     % now select the interval with the maximum standard deviation
    [~,interval] = max(ir);	% select cross power

elseif isfield(cfg,'detect_PEMOeuclid') && isyes(cfg.detect_PEMOeuclid)
    for ialter = 1:Nalter
        ir(ialter)=sqrt(sum(sum(sum(in{ialter}(:)-template(:))^2)));
%          cctemp = corrcoef(in{ialter}(:),template(:));
%          ir(ialter)=cctemp(1,2);
    end
   
%     % xcorr of internal representation and template
%     cc = sum(in.*repmat(template,1,size(in,2)));
%     % now select the interval with the maximum standard deviation
    [~,interval] = max(ir);	% select cross power
    % if another interval than the first has the max cross power, the interval is wrong, since work.signal always carries the
    % signal interval in the first column
%     if interval ~= 1
%         interval = 0;
%     end
elseif isfield(cfg,'detect_PEMOxcorr') && isyes(cfg.detect_PEMOxcorr)
    for ialter = 1:Nalter
        %%% NORMALIZED XCORR
        for ichan = 1:Nchan
            for imod = 1:Nmod
                [crosscorr(:,ichan,imod),lags]=xcorr(in{ialter}(:, ichan, imod),template(:, ichan, imod), ceil(cfg.detect_PEMOxcorr_limit)); %, period_sample/2
                [autocorr1(:,ichan,imod)]=xcorr( in{ialter}(:, ichan, imod), in{ialter}(:, ichan, imod), 0);
                [autocorr2(:,ichan,imod)]=xcorr( template(:, ichan, imod), template(:, ichan, imod), 0);
            end
        end
        xcorrsum(ialter,:) = sum(sum(crosscorr,2),3);
        autocorr1sum(ialter,:) = sum(sum(autocorr1,2),3);
        autocorr2sum(ialter,:) = sum(sum(autocorr2,2),3);
        ir(ialter)= max(xcorrsum(ialter,:))/sqrt(autocorr1sum(ialter,1)*autocorr2sum(ialter,1));
    end
    
    %     temp = reshape(template,[simwork.Nsample, simwork.Nchan, simwork.Nmod]);
    %     for i=1:size(in,2)
    %         stim = reshape(in(:,i),[simwork.Nsample, simwork.Nchan, simwork.Nmod]);
    %         period_sample = (1/work.exppar1)*(def.samplerate/simwork.cfg.downsampling_factor);
    %         for ichan = 1:simwork.Nchan
    %             for imod = 1:simwork.Nmod
    %                 [crosscorr(:,ichan,imod),lags]=xcorr(stim(:, ichan, imod),temp(:, ichan, imod)); %, period_sample/2
    %             end
    %         end
    %         xcorrsum(i,:) = sum(sum(crosscorr,2),3);
    %         %lags_of_interest = abs(lags)<period_sample/2;1:length(lags);%
    %         cc(i)= max(xcorrsum(i,:));%max(xcorrsum(i,lags_of_interest));
    %     end
    %cfg.display_step = 'yes';
    if isyes(cfg.display_out) | isyes(cfg.display_step)
        figure;
        for ialter = 1:Nalter
            plot(lags, xcorrsum(ialter,:)/sqrt(autocorr1sum(ialter,1)*autocorr2sum(ialter,1))); hold on
        end
        title('xcorr with template')
        legend({'target', 'nontarget'})
    end
    
    % now select the interval with the maximum standard deviation
    [~,interval] = max(ir);	% select cross power
    % if another interval than the first has the max cross power, the interval is wrong, since work.signal always carries the
    % signal interval in the first column
%     if interval ~= 1
%         interval = 0;
%     end
elseif isfield(cfg,'detect_EPSM') && isyes(cfg.detect_EPSM)
    % compute envelope power
    for ialter = 1:Nalter
        Epow(ialter) = mean(in{ialter}(:).^2,1);
    end
    Epow_template = mean(template(:).^2,1);
    
    
        % in the case of EPSM, 'in' is not used and 'template' is the
        % internal representation to be analyzed. Furthermore, the
        % target-present response is assumed to be 2
        ratio_Epow = 10*log10(Epow_template);%10*log10(Epow(1)/Epow(2));
        if ratio_Epow >= cfg.detect_EPSM_threshold
            interval = 2;
        else
            interval = 1; %randi(2)-1;
        end
        ir = Epow_template;
        %ir = Epow;
        
%         %%% correlation-based EPSM
%         [~,interval] = min(abs(Epow-Epow_template));
elseif isfield(cfg,'detect_PEMOthreshold') && isyes(cfg.detect_PEMOthreshold)
    % in the case of PEMO_threshold, 'template' is the internal
    % representation to be analyzed by correlation with 'in{2}'.
    % Furthermore, the target-present response is assumed to be 2
    if Nalter~=1
        %error('detect_PEMOthreshold requires one single interval');
    end
    cctemp = xcorr(in{2}(:),template(:));
    r = max(cctemp);
    %cctemp = corrcoef(in{2}(:),template(:));
    %r=cctemp(1,2);
    if r >= cfg.detect_PEMOthreshold_r
        interval = 2;
    else
        interval = 1; %randi(2)-1;
    end
    %     % xcorr of internal representation and template
    %     cc = sum(in.*repmat(template,1,size(in,2)));
    %     % now select the interval with the maximum standard deviation
    ir = r;% [~,interval] = max(ir);	% select cross power
    
end

end

