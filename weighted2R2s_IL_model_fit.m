% nargin inputs of form:
%   weighted.data (NvoxelsX x NvoxelsY x ... x Nechoes)
%   weighted.TEs
% varargout contains weightings extrapolated to TE=0 in the same order as 
% the input.

function [R2s,varexpl,varargout]=weighted2R2s(varargin)

dims=size(varargin{1}.data);

Nweighted=nargin;
Nvoxels=prod(dims(1:end-1));

%% Build design matrix
D=[];
for w=1:Nweighted
    d=zeros(length(varargin{w}.TEs),Nweighted+1);
    d(:,1)=-varargin{w}.TEs;
    d(:,w+1)=1;
    D=[D;d]; %#ok<AGROW>
end

%% Solve for R2s and coefficients
varargout=varargin; % preserve old fields
lnS=[];
for w=1:Nweighted
    
    rData=reshape(abs(varargin{w}.data),Nvoxels,length(varargin{w}.TEs));
    
    lnS=[lnS;log(rData).']; %#ok<AGROW>
end

% Well posed as sufficient rank; avoids backslash for data matrix below
DtDSlashDt=(D'*D)\D';

% ols solution
% extra unity in reshape argument to avoid problems if size(dims)==2.
R2s.data=reshape(DtDSlashDt(1,:)*lnS,[dims(1:end-1),1]);
for w=1:(nargout-2)    
    varargout{w}.TEs=0; % extrapolation to zero echo time

    extrapolatedData=reshape(exp(DtDSlashDt(w+1,:)*lnS),[dims(1:end-1),1]);
        
    % Minimal TE value used when there are errors from poor R2* fitting
    [~,minTEidx]=min(varargin{w}.TEs);
    weightedData=reshape(varargin{w}.data,Nvoxels,length(varargin{w}.TEs));
    weightedData=reshape(weightedData(:,minTEidx),[dims(1:end-1),1]);
    
    mask=isnan(extrapolatedData)|isinf(extrapolatedData)|(R2s.data<0);
    
    extrapolatedData(mask)=weightedData(mask);
    
    varargout{w}.data=extrapolatedData;
end

%%% IL: calculate model fit
all_betas = DtDSlashDt*lnS;
predicted_vals = D * all_betas;

%%% manual implementation of correlation coefficient, so that it is quicker
%%% than loop
top = sum((predicted_vals - mean(predicted_vals)) .* (lnS - mean(lnS)));
bottom = sqrt(sum((predicted_vals - mean(predicted_vals)).^2)) .* sqrt(sum((lnS - mean(lnS)).^2));
r2_quick = (top ./ bottom).^2;

% %% probably super inefficient
% display('starting inefficient r2s')
% parfor vox = 1:size(predicted_vals,2)
%     r = corrcoef(lnS(:,vox),predicted_vals(:,vox));
%     r2(vox) = r(1,2)^2;
% end
varexpl = reshape(r2_quick,[dims(1:end-1),1]);

end
