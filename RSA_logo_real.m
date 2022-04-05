%%
%Representational Analysis
clear
clc

load('RDM_logo_Real.mat');
N = length(RDM_logo_Real);
vec_RDM_logo = zeros([1,N*(N-1)/2]);
for i=1:N
    l = (i-1)*(2*N-i)/2+1;
    vec_RDM_logo(l:(l+N-i-1)) = RDM_logo_Real((i+1):end,i);
end

RSA_logoReal = cell(51,2);
for sub = 1:51
    if sub == 7 || sub == 11 || sub == 21 || sub == 31 || sub == 23
        continue;
    end
    
    fname = sprintf('Sub%03d_neuralRDM_Real.mat',sub);
    load(fname);

    [r1,p1] = corr( neuralRDM',vec_RDM_logo','type','Spearman','tail','right');
    
    RSA_logoReal{sub,1} = r1;RSA_logoReal{sub,2} = p1;

    fprintf('Sub%03d RSA has been calculated.\n',sub);
    
end

save RSA_logoReal RSA_logoReal
%%
clear
clc

load('RSA_logoReal.mat');

Zval = zeros(51,length(RSA_logoReal_parsg_fm{1,1}));
for sub = 1:51
        temp = RSA_logoReal_parsg_fm{sub,1};
        Zval(sub,:) = 0.5*(log((1+temp)./(1-temp)));       
end

[~,p,~,stat] = ttest(Zval,0,'Tail','Right');
Pval_fMRSK_real_fm = p';
Tval_fMRSK_real_fm = stat.tstat';

load('gfMask.mat');
rmask = load_untouch_nii('rVOT_aal.nii');
mask_ref = spm_select('List',pwd,'^rgrayTPM.*\.nii$');
mask_hdrs = spm_vol(mask_ref);
[maskdata ~] = spm_read_vols(mask_hdrs);
maskINDs = find((maskdata>0.2).*(logical(rmask.img)~=0).*(gfMask~=0));
aa = rmask.img(maskINDs);

data = load_untouch_nii('rgrayTPM_mask.nii');
[x,y,z] = ind2sub(size(data.img),maskINDs);
Orth = zeros(size(data.img));
logoRepre_real = zeros(size(data.img));
real = Tval_fMRSK_real_fm.*(Pval_fMRSK_real_fm<0.05);
real(isnan(real)) = 0;
for i = 1:length(Repres_orth)   
    logoRepre_real(x(i),y(i),z(i)) = real(i,1);
end
data.img = logoRepre_real;
save_untouch_nii(data,'RSA_logo_RW_VOT.nii');
