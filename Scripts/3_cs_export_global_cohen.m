%% This script loads the statistical output for the reference map (HD-EEG, 256 channels, Cohen's d) 
% and subtracts the results for the other channel densities

root_dir = '/home/uni10/nmri/projects/GoeEpiLongterm/nfocke/AllGGE_redux/results/';

% load maps for respective channel density and subtract from reference map
results_dir = '/home/uni10/nmri/projects/cstier/channels_analysis/results/orig/';
hdm = {'individual', 'canonical'};
metric = {'power', 'coh_img'}; % 'power'
density = {'full','192','128', '64', '48','32','IFCN1020_ext','IFCN1020'};
freq = {'m1', 'm2', 'm3', 'm4', 'm5', 'm6'};
freqname = {'delta', 'theta', 'alpha', 'beta1', 'beta2', 'gamma'};
file_abb = 'palm_out_dat_cohen_';
file_abb2 = 'palm_out_dat_tstat_fwep_';


for m = 1:length(metric)

 output_metric = [];
 
 for h = 1:length(hdm)  
  
  out_all = [];
  
  for d = 1:length(density) 
   for f = 1:length(freq)
    
   % load cohen's d and p-value for increased metrics in patients with IGE
   cohend = csvread(fullfile(root_dir, hdm{h}, density{d}, 'Cx-GGE/', metric{m}, 'palm_global', [file_abb freq{f} '_c2.csv']));
   pval = csvread(fullfile(root_dir, hdm{h}, density{d}, 'Cx-GGE/', metric{m}, 'palm_global', [file_abb2 freq{f} '_c2.csv']));
   p_orig = 10^(-(pval));
   
   resline{1} = hdm{h};
   resline{2} = density{d};
   resline{3} = freqname{f};
   resline{4} = cohend;
   resline{5} = p_orig;
   
   out_all = [out_all; resline];
   end
  end
  
  output_metric = [output_metric; out_all];
 end
   
 
 filename = [results_dir 'global_analysis_allchannels_' metric{m} '.mat'];
 save(filename, 'output_metric') 

end
  



