% cell array of weighted data structs:
%   weighted.data (NvoxelsX x NvoxelsY x ... x Nechoes)
%   weighted.TEs
%   weighted.fa
%
% Linear fit performed using linearisation proposed in Dathe and Helms,
% Physics in Medicine and Biology (2010)

function [A,R1,b1_corr_fa_maps,actual_fa_maps,R1error]=weighted2AR1_slice_profile_correction(weightedStructs,TR,method,relativeB1,slice_profile,sim_fas)

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
signalmatrix=[];
taumatrix=[];

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
            t=2*tan(relativeB1.*weighted.fa/2);
        case 'spr'
            %% IL: make B1 separately for echos: TO DO -> get error of estimation with polyfit to feed into R1 error calc
            slice_profile = slice_profile(:,1:19); %%% exclude highest FA
            sim_fas = sim_fas(1:19);
            s=size(slice_profile);
            for z=1:s(1)
                p(z,:) = polyfit(sim_fas,slice_profile(z,:),4);
                Synth(z,:)=p(z,4)*sim_fas+p(z,3)*sim_fas.^2+p(z,2)*sim_fas.^3+p(z,1)*sim_fas.^4;
                b1corr_fa_deg = rad2deg(relativeB1.*weighted.fa);
                spc_fa_deg(:,:,z)=p(z,4)*b1corr_fa_deg(:,:,z)+p(z,3)*b1corr_fa_deg(:,:,z).^2+p(z,2)*b1corr_fa_deg(:,:,z).*3+p(z,1)*b1corr_fa_deg(:,:,z).^4;
            end;
            t = 2*tan(deg2rad(spc_fa_deg)/2);
            b1_corr_fa_maps(:,:,:,w) = b1corr_fa_deg;
            actual_fa_maps(:,:,:,w) = spc_fa_deg; %%% spit out actual fa maps
        otherwise
            error('unknown method')
    end
    
    y=[y,weighted.data(:)./t(:)]; %%% signal div tau
    
    D=[D,-weighted.data(:).*t(:)/2]; %%% signal mul tau
    
    signalmatrix=[signalmatrix,weighted.data(:)];
    
    taumatrix=[taumatrix,t(:)];
end

mask=(weighted.data>0)&(relativeB1>1e-1);

Nvoxels=nnz(mask);

y=y(mask,:);
D=D(mask,:);
signalmatrix=signalmatrix(mask,:);
taumatrix=taumatrix(mask,:);

%% Solve linearly to get rho
a=zeros(Nvoxels,1);
recrho=zeros(Nvoxels,1);
parfor n=1:Nvoxels
    %b=robustfit(D(n,:),y(n,:),'bisquare');
    b=y(n,:)/[ones(1,Nweighted);D(n,:)]; %%% just fits sig div tau by sig mul tau
    a(n)=b(1); %%% offset
    recrho(n)=b(2); %%% slope
end

%% IL: calculate error margins
if size(signalmatrix,2) == 2 %%% only then we have calculated
    S1 = signalmatrix(:,1);
    T1 = taumatrix(:,1);
    S2 = signalmatrix(:,2);
    T2 = taumatrix(:,2);
    alpha1 = weightedStructs{1}.fa;
    alpha2 = weightedStructs{2}.fa;
    %%% use actual FA
    alpha1 = deg2rad(actual_fa_maps(:,:,:,1)); alpha1 = alpha1(:); alpha1 = alpha1(mask);
    alpha2 = deg2rad(actual_fa_maps(:,:,:,2)); alpha2 = alpha2(:); alpha2 = alpha2(mask);
    %%% define errors
    error_alpha1 = .02 * alpha1; %%% is 2% error in FA1 calc
    error_alpha2 = .02 * alpha2;
    %%% calculate derivatives
    dRho_dT1 = (S1.*T2 .* (-2*S1.*T1.*T2 + S2.*(T1.^2 + T2.^2)) ) ./ (2*(S2.*T1 - S1.*T2).^2 ); %%% this is partial derivative of Dathe Helms 2010 equation17, as calculated by PS with mathematica
    dRho_dT2 = (S2.*T1 .* (-2*S2.*T1.*T2 + S1.*(T1.^2 + T2.^2)) ) ./ (2*(S2.*T1 - S1.*T2).^2 ); %%% this is partial derivative of Dathe Helms 2010 equation17, as calculated by PS with mathematica
    %%% wrong original:
    %d_T1 = (sec(alpha1))^2 * error_alpha1; %%% should be one value
    %d_T2 = (sec(alpha2))^2 * error_alpha2;
    d_T1 = (sec(alpha1/2)).^2 .* error_alpha1; %%% should be one value
    d_T2 = (sec(alpha2/2)).^2 .* error_alpha2;
    %%% calculate error in R1
    d_R1 = abs(dRho_dT1) .* abs(d_T1) + abs(dRho_dT2) .* abs(d_T2);
    R1error = weighted;
    R1error.data = 0*R1error.data;
    R1error.data(mask) = d_R1;
end
%%

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
        case {'exact','smallFlipAngle','spr'}
            % Derived from Equation (11) of Dathe and Helms, Physics in Medicine and Biology (2010)
            recrho(recrho<1)=inf; % Remove values which can't give sensible results
            R1.data(mask)=(2/TR)*acoth(2*recrho); %%% technically only necessary if TR/T1 < 1?
        otherwise
            error('unknown method')
    end
end

end
