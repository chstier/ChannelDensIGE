%% This scripts runs the spin-test for the correlation of results across channel densities
% Based on the toolbox by Alexander-Bloch et al., 2018: https://doi.org/10.1016/j.neuroimage.2018.05.070
%
% Step 1: 
% Load the statistical map for 256-channel HD-EEG as a reference (here t-stats for increased power/connectivity (C2)
% Spin data for the left and right hemisphere seperately (cs_SpinPermuFS)

% Step 2:
% Load the statistical maps for other channel densities and compute the real Spearman's correlation with the statistical map 
% for the 256-channel HD-EEG. Test Spearman's correlation coefficient against the null distribution created using SpinPermuFS and save output
% 
% written and edited by Christina Stier, 2022

%% STEP 1 
% prepare permutation-files for individual headmodels (which needs to be done for 2 metrics, and 8 densities)
metric = {'power', 'coh_img'}; % 
density = {'full', '192','128', '64', '48','32','IFCN1020_ext','IFCN1020'};

% prepare all the files and folders needed
suma_map = fullfile(pwd, 'conf/suma-all-fsaverage-10.mat'); % get SUMA incl. subcortical nuclei (ld 10)
suma_all_all = load(suma_map);
suma_cortex = fullfile(pwd, 'conf/suma-cortex-fsaverage-10.mat'); % get SUMA without subcortical nuclei
load(suma_cortex);

mkdir 'Rotation_corrmaps_hdm' 
rootfolder = ('/home/AllGGE_redux/results/');

% loop over metrics and densities and create permutations
for m = 1:length(metric)
 for d = 1:length(density)
  
  % load results for individual headmodels for the respective density
  indiv = load_mgh(fullfile(rootfolder, 'individual' , density{d} ,'/Cx-GGE/', metric{m}, '/palm_surface/palm_out_dpv_cohen_m2_c2.mgz'));

  % only use cortical nuclei
  indiv = indiv(1:2004,:);

  % apply mask (to make sure that all medial points are NaN) 
  indiv(suma_all.msk == 0) = NaN;

  indiv_left = indiv(suma_all.hemi == 1,:); % subcortical values are already zeroed out
  indiv_right = indiv(suma_all.hemi == 2,:);

  readleft = indiv_left;
  readright = indiv_right;

  % choose number of permutations:
  permno = 1000;
  perm = '1000';

  % spin and generate 0 distribution
  wsname = fullfile(pwd, 'Rotation_corrmaps_hdm', [density{d} '_cohen_m2_c2_individual_hdm_' metric{m} '_' perm '_perm_spearman.mat']); % change filename of output if needed
  cs_SpinPermuFS(readleft, readright, permno, wsname) 
  
 end
end

%% STEP 2
% load maps for the canonical headmodels and do spin-correlations for each density and metric
rootfolder = '/home/AllGGE_redux/results/';
metric = {'power', 'coh_img'}; % 'power'{'power', 'coh_img'}
density = {'full', '192','128', '64', '48','32','IFCN1020_ext','IFCN1020'};
file_abb = 'palm_out_dpv_cohen_m2_c2.mgz';
perm = '1000';
permno = 1000;

results_dir = fullfile(pwd, 'Rotation_corrmaps_hdm');

for m = 1:length(metric)
 for d = 1:length(density)
   
   % load results for the canonical headmodel
   canon_file = fullfile(rootfolder, 'canonical', density{d}, '/Cx-GGE/', metric{m}, '/palm_surface/palm_out_dpv_cohen_m2_c2.mgz');
   canon_name = ['canonical_' metric{m} '_' density{d} '_cohen_m2_c2'];
   canon = load_mgh(canon_file);
   canon = canon(1:2004,:); % delete values for subcortical nuclei
   canon(suma_all.msk == 0) = NaN; % apply mask (to make sure that all medial points are NaN)
   
   % load original reference maps (individual hdm / metric separately)
   indiv = load_mgh(fullfile(rootfolder, 'individual' , density{d} ,'/Cx-GGE/', metric{m}, '/palm_surface/palm_out_dpv_cohen_m2_c2.mgz'));
   indiv = indiv(1:2004,:); % only use cortical nuclei
   indiv(suma_all.msk == 0) = NaN; % apply mask (to make sure that all medial points are NaN) 
   
   wsname = fullfile(pwd, 'Rotation_corrmaps_hdm', [density{d} '_cohen_m2_c2_individual_hdm_' metric{m} '_' perm '_perm_spearman.mat']); % change filename of output if needed

   [pval, realrho, realpval, nullroh] = cs_pvalvsNull_spearman(indiv, canon, permno,wsname); 
   
   % save results in struct and csv
   stats(d).pval = pval;
   stats(d).real_rho = realrho;
   stats(d).realpval = realpval;
   stats(d).mod = canon_name;
   
   filename_stats = fullfile(results_dir, ['hdm_corr_' metric{m} '_theta_cohen_c2_' perm '_perm_spearman.csv']);
   
   writetable(struct2table(stats), filename_stats);
   
 end
end
