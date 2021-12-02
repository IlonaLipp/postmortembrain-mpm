% cell array of weighted data structs:
%   weighted.data (NvoxelsX x NvoxelsY x ... x Nechoes)
%   weighted.TEs
%   weighted.fa
%
% Linear fit performed using linearisation proposed in Dathe and Helms,
% Physics in Medicine and Biology (2010)

function [A,R1]=weighted2AR1_separate_B1s(weightedStructs,TR,method,relativeB1s) %%% just used original function of luke and changed from one map to one map per contrast as input

if ~exist('relativeB1','var')
    relativeB1=1;
end

Nweighted=length(weightedStructs);

TEs=weightedStructs{1}.TEs;
if length(TEs)>1
    warning('Multiple TEs detected; echoes will be averaged and A will be non-quantitative!')
end

%% Build design matrix
y=[];
D=[];
for w=1:Nweighted
    
    assert(isequal(weightedStructs{w}.TEs,TEs),'Echoes must match between different contrasts!')
    
    % Deal with multiple echoes by averaging;note that this will incorrectly
    % scale the A parameter!
    weighted=averageEchoes(weightedStructs{w},length(weightedStructs{w}.TEs));
    
    switch method
        case 'Helms2008'
            t=weighted.fa;
        case 'smallFlipAngle'
            t=relativeB1.*weighted.fa;
        case {'exact','smallTR*R1'}
            t=2*tan(relativeB1s{w}.*weighted.fa/2);
        otherwise
            error('unknown method')
    end
    
    y=[y,weighted.data(:)./t(:)]; %#ok<AGROW>
    
    D=[D,-weighted.data(:).*t(:)/2]; %#ok<AGROW>
end

mask=(weighted.data>0)&(relativeB1s{1}>1e-1);

Nvoxels=nnz(mask);

y=y(mask,:);
D=D(mask,:);

%% Solve for R2s and coefficients
a=zeros(Nvoxels,1);
recrho=zeros(Nvoxels,1);
parfor n=1:Nvoxels
    %b=robustfit(D(n,:),y(n,:),'bisquare');
    b=y(n,:)/[ones(1,Nweighted);D(n,:)];
    a(n)=b(1);
    recrho(n)=b(2);
end

A=weighted; % preserve fields
A.data=0*A.data;
A.fa=[]; % no well-defined flip angle
A.data(mask)=a;

switch method 
    case 'Helms2008'        
        if length(relativeB1)>1, relativeB1=relativeB1(mask); end
        A.data(mask)=A.data(mask)./relativeB1;
end

if nargout>1
    % Remove values which can't give sensible results
    recrho(isnan(recrho))=inf;
    
    R1=weighted; % preserve fields
    R1.data=0*R1.data;
    R1.fa=[]; % no well-defined flip angle
    
    switch method
        case 'Helms2008'
            R1.data(mask)=relativeB1(mask).^2./(recrho*TR);
        case 'smallTR*R1'
            R1.data(mask)=1./(recrho*TR);
        case {'exact','smallFlipAngle'}
            % Derived from Equation (11) of Dathe and Helms, Physics in Medicine and Biology (2010)
            recrho(recrho<1)=inf; % Remove values which can't give sensible results
            R1.data(mask)=(2/TR)*acoth(2*recrho);
        otherwise
            error('unknown method')
    end
end

end
