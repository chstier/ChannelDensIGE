% load reference maps: output for statistical analysis of HD-EEG 
% 256 channels (cohen's d)

analysis_dir = '/home/uni10/nmri/projects/cstier/channels_analysis';
results_analysis = '/home/uni10/nmri/projects/cstier/channels_analysis/results';

% in this case only for the theta band
rootfolder = ('/home/uni10/nmri/projects/GoeEpiLongterm/nfocke/AllGGE_redux/results/');
% refmap = load_mgh('/home/uni10/nmri/projects/GoeEpiLongterm/nfocke/AllGGE_redux/results/individual/full/Cx-GGE/power/palm_surface/palm_out_dpv_cohen_m2_c2.mgz');

% load maps for respective channel density and subtract from reference map
results_dir = '/home/uni10/nmri/projects/GoeEpiLongterm/nfocke/AllGGE_redux/results/';
hdm = {'individual', 'canonical'};
metric = {'power', 'coh_img'}; % 'power'
density = {'192','128', '64', '48','32','IFCN1020_ext','IFCN1020'};
file_abb = 'palm_out_dpv_cohen_m2_c2.mgz';

% plotting settings

% load([analysis_dir '/conf/suma-all-fsaverage-10.mat'],'suma_all')
plot_surface = load(fullfile(pwd,'conf',['suma-all-fsaverage-40.mat'])); 
ref_surface = load(fullfile(pwd,'conf',['suma-all-fsaverage-10.mat'])); 
remap_matrix=[];
[~, ~, remap_matrix.vertices, remap_matrix.weights ]=nmri_suma_surf2surf_transform(ref_surface.suma_all,plot_surface.suma_all, zeros([size(ref_surface.suma_all.pos,1),1]));


opt=[];
opt.per_hemi=1;
opt.per_cortex=0;
opt.rot=[90 0 ; -90 0];
% opt.thresh=1.3;
opt.clim=[-0.6 0.6];
opt.colormap='parula';
opt.colorbar='parula';
opt.scale=1;


for m = 1:length(metric)
 for h = 1:length(hdm)
  
  output=cell(length(density)+1,4);
  output{1,1} = 'metric';
  output{1,2} = 'headmodel';
  output{1,3} = 'channel density';
  output{1,4} = 'corr with ref. map';
  
  for d = 1:length(density)
   
   % load the reference map and map of interest
   refmap = load_mgh(fullfile(rootfolder, hdm{h}, 'full/Cx-GGE/', metric{m}, '/palm_surface/palm_out_dpv_cohen_m2_c2.mgz'));
   map = load_mgh(fullfile(results_dir, hdm{h}, density{d}, 'Cx-GGE', metric{m}, 'palm_surface', file_abb));
   diff = refmap-map;
  
   % plot difference
   opt.output = fullfile(results_analysis, ['diff_', density{d}, '_', metric{m}, '_', hdm{h}, '.png']);
   
   diff_map = nmri_suma_surf2surf_transform(ref_surface.suma_all,plot_surface.suma_all, diff, remap_matrix.vertices, remap_matrix.weights );

   hFig = cs_nmri_plot_surface_suma(plot_surface.suma_all, diff_map, opt);
   
   
   % use rank correlation (to get real rho, but do significance test via spin-test)
   [realrho, realpval] = corr(refmap,map, 'rows','complete', 'type', 'Spearman');
   
   output{d+1, 1} = metric{m};
   output{d+1, 2} = hdm{h};
   output{d+1, 3} = density{d};
   output{d+1, 4} = realrho;
   

  end
  
  % save results
  output = cell2table(output);
  filename = fullfile(results_analysis, ['correlation_d_ref_channels_theta_', metric{m}, '_', hdm{h} '.csv']);
  writetable(output, filename, 'WriteVariableNames',0)
  
  % combine plots
  name = [results_analysis, '/diff_contrast_publ2_' metric{m} '_' hdm{h} '.png'];
  img1 = imread([results_analysis '/diff_' density{1} '_' metric{m} '_' hdm{h} '.png']);
  img2 = imread([results_analysis '/diff_' density{2} '_' metric{m} '_' hdm{h} '.png']);
  img3 = imread([results_analysis '/diff_' density{3} '_' metric{m} '_' hdm{h} '.png']);
  img4 = imread([results_analysis '/diff_' density{4} '_' metric{m} '_' hdm{h} '.png']);
  img5 = imread([results_analysis '/diff_' density{5} '_' metric{m} '_' hdm{h} '.png']);
  img6 = imread([results_analysis '/diff_' density{6} '_' metric{m} '_' hdm{h} '.png']);
  img7 = imread([results_analysis '/diff_' density{7} '_' metric{m} '_' hdm{h} '.png']);
  img = [img1; img2; img3; img4; img5; img6; img7];
  imwrite (img, name);

 end
end
