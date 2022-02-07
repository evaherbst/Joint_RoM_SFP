function [p,q,f] = parseROMTrialData(specimen, subset)

%parseROMTrialData
%
%   [p,q] = parseROMTrialData(specimen) returns the marker positions and
%       orientations for all subtrials matching the specimen filename.
%
%   p is Nx3 XYZ positions over N frames, and q is Nx4 quaternions.
%
%   For example, the specimen SS01_FRZN_RH_KNEE has many subtrials, e.g.
%       SS01_FRZN_RH_KNEE_1_L.mat
%       SS01_FRZN_RH_KNEE_2_F.mat
%       SS01_FRZN_RH_KNEE_3_A.mat
%       SS01_FRZN_RH_KNEE_4_AFL.mat
%       ...
%
%   [p,q] = parseROMTrialData('SS01_FRZN_RH_KNEE') will return data from
%       all subtrails named 'SS01_FRZN_RH_KNEE*.mat'. The positions and
%       orientations are concatenated back-to-back.
%
%   [p,q] = parseROMTrialData(specimen, subset) returns the positions and
%       orientations for the subtrials specified by subset only. If subset
%       is empty (subset = []), it returns all trials.
%
%   subset is an array of desired trial indices. For example:
%
%   [p,q] = parseROMTrialData('SS01_FRZN_RH_KNEE', [1 2 3]) will return
%       data from only the first three subtrials
%
%   2018 Enrico Eberhard

root = '../Data/';

files = dir([root specimen '*.mat']);


trials = {files.name};
fprintf('%i files found matching %s*.mat.\n', numel(trials), specimen);

%#ok<*AGROW>
p = zeros(0,3);
q = zeros(0,4);
f = zeros(0,6);

if nargin < 2 || isempty(subset)
    subset = 1:numel(trials);
    fprintf('Gathering data from all files:\n');
else

    %limit subset to available range of trials
    subset = subset(subset<=numel(trials));
    subset = subset(subset>=1);
    
    fprintf('Gathering data from following files:\n');
end

for t = subset
    
    
    %ignore high frequency data format for now
    if strfind(trials{t},'_Q') == length(trials{t}) - 5
        continue
    end
    
    fprintf('%s\n',trials{t});
    
    data = load([root trials{t}]);
    
%     strains = data.trial.strains.data - data.trial.strains.offset;
%     
%     fs = 150;
%     s_r = reshape(strains,fs,round(length(strains)/fs),6);
%     
%     f_t = squeeze(mean(s_r))*data.trial.strains.mat;
    

    strains = data.trial.strains.data;
    fs = 150;
	s_r = reshape(strains,fs,round(length(strains)/fs),6);
    f_t = squeeze(mean(s_r));
    
    d6D = [data.trial.states.d6D];
    
    body = [d6D.body];
    
    p_t = [[body.X]' [body.Y]' [body.Z]'];
    
    if any(p_t==0)
        continue
    end
    
    q_t = zeros(length(body),4);
    for bod = 1:length(body)
        q_t(bod,:) = quatFromRot(body(bod).R);
    end
    
    %get rid of NaNs
    p_t = p_t(~isnan(q_t(:,1)),:);
    q_t = q_t(~isnan(q_t(:,1)),:);
    f_t = f_t(~isnan(q_t(:,1)),:);
    
    
    %filter params
    wnd = 10;
    b = (1/wnd)*ones(1,wnd);
    a = 1;
    
    
    p_f  = [ones(wnd,3).*p_t(1,:); p_t;...
                        ones(wnd,3).*p_t(end,:)];
    q_f  = [ones(wnd,4).*q_t(1,:); q_t;...
                        ones(wnd,4).*q_t(end,:)];
    f_f  = [ones(wnd,6).*f_t(1,:); f_t;...
                        ones(wnd,6).*f_t(end,:)];
                    

    for ii = 1:3
        p_f(:,ii)  = filter(b,a,p_f(:,ii));
        q_f(:,ii)  = filter(b,a,q_f(:,ii));
    end
    q_f(:,4)  = filter(b,a,q_f(:,4));
    
    for ii = 1:6
        f_f(:,ii)  = filter(b,a,f_f(:,ii));
    end

    p_ff = p_f(wnd+ceil(wnd/2):end-ceil(wnd/2),:);
    q_ff = q_f(wnd+ceil(wnd/2):end-ceil(wnd/2),:);
    f_ff = f_f(wnd+ceil(wnd/2):end-ceil(wnd/2),:);
    
    
    p = [p; p_ff]; 
    q = [q; q_ff]; 
    f = [f; f_ff]; 

end

%now, offset by the average overall strains and convert to forces
f = (f - mean(f))*data.trial.strains.mat;


end