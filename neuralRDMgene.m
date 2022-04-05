clear
clc
% Generate neural RDM

cwd = '~mainfolder';
cd(cwd);
load('ctrRelSphereSUBs.mat');
load('trilSUB.mat');

load('gfMask.mat');
rmask = load_untouch_nii('rVOT_aal.nii');
mask_ref = spm_select('List',pwd,'^rgrayTPM.*\.nii$');
mask_hdrs = spm_vol(mask_ref);
[maskdata,~] = spm_read_vols(mask_hdrs);
maskINDs = find((maskdata>0.2).*(logical(rmask.img)~=0).*(gfMask~=0));
aa = rmask.img(maskINDs);

for sub = 1:51  
    cd([cwd sprintf('/Sub%03d',sub)]);% First-level results for each trial
    gray_stimuli_neuralap = cell(length(maskINDs),40);
    
    for stimi = 1:40
        file_filter = sprintf('spmT_%04d',stimi);
        ref = spm_select('List',pwd,['^' file_filter '.*\.nii$']);
        hdrs = spm_vol(ref);
        [data,~] = spm_read_vols(hdrs);
        
        volSize_vox = size(data);
              
        for searchlightVoxel = 1:length(maskINDs)
        
                fprintf('\n Sub%03d  Stimuli %03d  voxel %01d percent is now calculated\n',sub,stimi,floor(searchlightVoxel*100/length(maskINDs)));
                
                [x,y,z] = ind2sub(volSize_vox,maskINDs(searchlightVoxel));
                
                % compute (sub)indices of (vox)els (c)urrently (ill)uminated by the spherical searchlight
                cIllVoxSUBs=repmat([x,y,z],[size(ctrRelSphereSUBs,1) 1])+ctrRelSphereSUBs;

                % exclude out-of-volume voxels
                outOfVolIs=(cIllVoxSUBs(:,1)<1 | cIllVoxSUBs(:,1)>volSize_vox(1)|...
                            cIllVoxSUBs(:,2)<1 | cIllVoxSUBs(:,2)>volSize_vox(2)|...
                            cIllVoxSUBs(:,3)<1 | cIllVoxSUBs(:,3)>volSize_vox(3));

                cIllVoxSUBs=cIllVoxSUBs(~outOfVolIs,:);

                % list of (IND)ices pointing to (vox)els (c)urrently (ill)uminated by the spherical searchlight
                cIllVox_volINDs=sub2ind(volSize_vox,cIllVoxSUBs(:,1),cIllVoxSUBs(:,2),cIllVoxSUBs(:,3));

                gray_stimuli_neuralap{searchlightVoxel,stimi} = data(cIllVox_volINDs);
                
        end
    end
    
    neuralRDM = zeros(length(maskINDs),780);
    for i = 1:length(maskINDs)
        
        for j = 1:780
            [r,p] = corr(gray_stimuli_neuralap{i,trilSUB(j,1)},...
                         gray_stimuli_neuralap{i,trilSUB(j,2)},...
                         'Type','Spearman');
            neuralRDM(i,j) = r;
        end
    
    end
    filename = sprintf('Sub%03d_neuralRDM_real',sub);
    save(filename,'neuralRDM');
    clear i sti j gray_stimuli_neuralap neuralRDM
    
end