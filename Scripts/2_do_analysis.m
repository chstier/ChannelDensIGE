%% example streamlined analysis

% go to analysis dir
cd('/home/uni10/nmri/projects/GoeEpiLongterm/nfocke/AllGGE_redux') % this needs to me modified for your needs


if ~exist('nmri_all_subjects','file')
 % add scripts paths, if not yet present
 addpath(genpath('scripts'))
end

% call the params include
nmri_include_read_params


params.EEG.headmodel = 'openmeeg'; % fix to OpenMEEG
%params.EEG.headmodel = 'dipoli'; % fix to dipoli

% now get all subjects, cached if possible
all_subjects=nmri_all_subjects(pwd,1);


% now make sure we have vigilance scoring done
filters=[]; % make empty filter
filters.dtype='EEG'; % EEG 
filters.stamps.vigilance=true; % we need vigilance scoring done
filters.stamps.MRI_ctf=true; % and MRI_ctf done
filtered_subjects = nmri_filter_subjects ( all_subjects, filters );

%% Now do headmodelling, if not done
% in this example we also make a HDM (OpenMEEG)


% Freesurfer storage dir
if ~isempty(getenv('SUBJECTS_DIR'))
 fsdir = getenv('SUBJECTS_DIR');
else
 % default to NMRI/PATLAN behaviour
 fsdir = '/data/freesurfer/6.0.0';
end

% check SUMA, files were copied from all different directories and we want
% a valid SUMA for all
for n=1:length(filtered_subjects)
 % check if done
 subject=filtered_subjects{n};
 % check SUMA
 if isfield(subject,'suma_dir') && ~exist(subject.suma_dir,'dir')
  % not there, so unset
  subject=rmfield(subject,'suma_dir');
  if isfield(subject,'suma_surface')
   subject=rmfield(subject,'suma_surface');
  end
 end
 if isfield(subject,'suma_surface') && isempty(subject.suma_surface)
  subject=rmfield(subject,'suma_surface');
 end
 
 
 if (~isfield(subject,'suma_dir') || ~exist(fullfile(subject.fsdir,'mri','nu.mgz'),'file'))
  % search for matching fressurfer IDs
  res=dir(fullfile(fsdir,[ subject.id '*']));
  possible={};
  for i=1:length(res)
   % check SUMA dir
   if (exist(fullfile(fsdir,res(i).name,'SUMA','std.40.lh.white.gii'),'file'))
    %seems processed
    possible{end+1}=res(i).name;
   else
    disp(fullfile(fsdir,res(i).name,'SUMA','std.40.lh.white.gii'))
   end
  end
  if (length(possible)==1)
   % then take it
   subject.fsid=possible{1};
   subject.suma_dir=fullfile(fsdir,subject.fsid,'SUMA');
   subject.fsdir=fullfile(fsdir,subject.fsid);
  end
  if (length(possible)>1)
   % more than one, ask user to pick
   qcfg=[];
   qcfg.question={['Found >1 possibility, please pick one (Subject=' subject.id  ')']};
   qcfg.title='SUMA selection';
   qcfg.options=possible;
   qcfg.mode='popup';
   manselect=nf_uidialog(qcfg);
   %manselect=listdlg('Name',,'SelectionMode','single','ListSize',[300 (50+(length(possible)*10))],'ListString',possible);
   if (~isempty(manselect))
    subject.fsid=manselect;
    subject.suma_dir=fullfile(fsdir,subject.fsid,'SUMA');
    subject.fsdir=fullfile(fsdir,subject.fsid);
   end
  end
  % if nothing, skip this subject
  if ~isfield(subject,'suma_dir') || ~exist(fullfile(subject.suma_dir,'std.40.lh.white.gii'),'file')
   fprintf('No SUMA found for %s -- skipping this one\n',subject.id)
   subject={};
  end
 end
 filtered_subjects{n}=subject;
end
% delete empty
filtered_subjects(cellfun(@isempty,filtered_subjects))=[];
 

% run local since compile does not work well

% remember failed subjects
failed_subjects={};
for n=1:length(filtered_subjects)
 % check if done
 subject=filtered_subjects{n};
 if ~exist(fullfile(subject.analysis_dir,subject.id,'processed',['hdm_lead_' subject.id '_' subject.exam_id '_individual_' params.EEG.headmodel '.mat']),'file')
  % not there, try to run it
  try
   filtered_subjects{n}=nmri_make_hdm_suma(filtered_subjects{n},[],params);
  catch
   fprintf('Subject=%s has failed!\n',subject.id)
   failed_subjects=[failed_subjects filtered_subjects{n}];
  end
 else
  fprintf('Subject=%s already done\n',subject.id)
 end
 % close the windows after each run
 close all
end

% print out what subjects have failed
failed_ids=cellfun(@(x) x.id,failed_subjects,'UniformOutput',false);
fprintf('Headmodelling failed for: %s\n',failed_ids{:})

%% check the QC plots of all heamdodels

all_subjects=nmri_all_subjects(pwd,1);
% now re-filter
filters=[]; % make empty filter
filters.stamps.(['hdmSUMA_individual_' params.EEG.headmodel])=true; % select only successfull hdm lead cases 
filtered_subjects = nmri_filter_subjects ( all_subjects, filters );

checked_subjects={};
failed_subjects={};

for i=1:length(filtered_subjects)
 if ~ismember(filtered_subjects{i}.id,checked_subjects)
  hFig=figure;
  imshow(imread(fullfile(filtered_subjects{i}.QCdir,['hdm_sensors_plot_' filtered_subjects{i}.id '_' filtered_subjects{i}.exam_id '_' params.EEG.headmodel '.png'])),'InitialMagnification',80)
  resp=questdlg('Are you happy with this QC plot?');
  if strcmp(resp,'Yes')
   checked_subjects=[checked_subjects filtered_subjects{i}.id];
  elseif  strcmp(resp,'No')
   failed_subjects=[failed_subjects filtered_subjects{i}.id];
  end
  close(hFig);
 end
end

% print groups of failed
all_ids=cellfun(@(x) x.id,filtered_subjects,'UniformOutput',false);
for i=1:length(failed_subjects)
 fprintf('Failed subject=%s, group=%s\n',failed_subjects{i},filtered_subjects{find(strcmp(failed_subjects{i},all_ids))}.group)
end

% delete MRIs of failed


%% now re-load all subjects
all_subjects=nmri_all_subjects(pwd,1);

% now re-filter
filters=[]; % make empty filter
filters.stamps.(['hdmSUMA_individual_' params.EEG.headmodel])=true; % select only successfull hdm lead cases 
filtered_subjects = nmri_filter_subjects ( all_subjects, filters );

fprintf('Have found N=%d subjects with processed headmodel\n',length(filtered_subjects))


% now run processing, if needed
todo={};
for n=1:length(filtered_subjects)
 % check if done
 subject=filtered_subjects{n};
 if ~exist(fullfile(subject.analysis_dir,subject.id,'stats',['coh_img_' subject.id '_' subject.exam_id '_individual_' params.EEG.headmodel '_' params.freqsNames{end} '.mat']),'file')
  % not there, add to todo list
  todo=[todo subject];
 end
end
fprintf('Have found N=%d subjects that need processing\n',length(todo))

% setup the config
cfg=[]; % we can stick to the default config
% use a fixed compile hash
cfg.compiled_hash={'8827d0c47d6dace31d204c8cc2b1c8a7'};
cmd={'nmri_processing'}; % give the command(s) as cell array

SGEdirs=nmri_run_multiple(todo,cfg,cmd,params);


%% You may want to check the status of your jobs, can be run repeately

nmri_check_jobs(SGEdirs,'-showerr');

% or with host infos
nmri_check_jobs(SGEdirs,'-showhosts');

%% now re-load all subjects
all_subjects=nmri_all_subjects(pwd,1);

% now re-filter
filters=[]; % make empty filter
filters.stamps.(['processing_individual_' params.EEG.headmodel])=true; % select only successfull hdm lead cases 
processed_subjects = nmri_filter_subjects ( all_subjects, filters );

fprintf('Have found N=%d subjects with processing done\n',length(filtered_subjects))

% if not all subjects are processed / differnt N, go back to previous steps
% and investigate

% rememeber the processed subjects
save(['processed_subjects' params.EEG.headmodel '.mat'],'processed_subjects')


%% now we go to the "real" analysis, start with the export of the maps
% load processed_subjects from disk if needed
%if ~exist('processed_subjects','var')
 load(['processed_subjects' params.EEG.headmodel '.mat'],'processed_subjects')
%end


% now make an export of all the needed maps

opt=[];
opt.metrics={'coh_img','power','coh_real'}; % will use imaginary part of coherency
opt.scale={'none'}; % no global scaling
opt.global={'abs'}; % take absolute values
opt.dim_reduction={'mean'}; % use the mean to reduce to one value per vertex
opt.save_mgh=false; % we take care of making the .mgh files late
opt.save_grp_mean=true; % save a mean map
opt.hdm_class=['individual_' params.EEG.headmodel];
opt.output=fullfile(pwd,'export',['all_subjects_' params.EEG.headmodel]);
nmri_export_metrics(processed_subjects,opt);



%% Now deal with regressors, and generate a unified .mat file

all_regs=dir('../Project_Filters/*.xlsx');
% Excel files, downloaded from the internal NMRI server
% the principle would equally work with the .tsv files generated by the
% export_blinded approach. You need to adopt the paths and joker (*)
% extension for you case for sure

all_regs_full=cell(1,length(all_regs));
for i=1:length(all_regs)
 all_regs_full{i}=fullfile(all_regs(i).folder,all_regs(i).name);
end

% read the groups from tables
allGrps=nmri_import_group_from_tables(all_regs_full);

save('allGrps.mat','-struct','allGrps')


% now loop over all found files and read each as an NMRI regressor
all_reg_file={};
for i=1:length(all_regs)
 [pa, fi, ext]=fileparts(all_regs(i).name);
 this_reg=nmri_import_table_as_regressor(all_regs_full{i},fullfile(pwd,[fi '.mat']));
 all_reg_file=[all_reg_file {fullfile(pwd,[fi '.mat'])}];
end

% now concatenate all the regressors into one file
ureg=nf_concat_regressors(all_reg_file,'concated_regressor.mat');

% note: the tables used here only contain the "basename" as identifying
% field. This will be fine, if each subject has only one examination and
% the values in the table (like age, medication, seizures,...) are valid
% for the examination in question. In an analyis with multiple
% examinations, this logic needs to be checked. If the timepoint logic is
% used (..A,..B,...) the principle will work. If one subject/timepoint
% combination has multiple values, you need to modify this script to make
% sure each examination has a unique "fname" otherwise, you will have
% problems later. You could combine use [subject.id '_' subject.exam_id] to
% make a merged fname.


% now check if we have everything and make our custom regressor (with just
% the needed fields)
reg=[];
reg.fname={};
reg.age=[];
reg.grp_sex=[];
reg.grp_siteTuePrisma=[];
reg.grp_siteTueTrio=[];
reg.grp_siteGoePrisma=[];

for i=1:length(processed_subjects)
 if ~any(strcmp(processed_subjects{i}.id,ureg.fname))
  error('Missing subject=%s from both regressors (or multiple hits)!!\n',processed_subjects{i}.id)
 else
  idx=find(strcmp(processed_subjects{i}.id,ureg.fname));
  reg.fname{end+1,1}=ureg.fname{idx};
  reg.age(end+1,1)=ureg.age(idx);
  reg.grp_sex(end+1,1)=ureg.grp_sex(idx);
  reg.grp_siteTuePrisma(end+1,1)=ureg.grp_mr_scanner_SIEMENS_Prisma(idx);
  reg.grp_siteTueTrio(end+1,1)=ureg.grp_mr_scanner_SIEMENS_TrioTim(idx);
  reg.grp_siteGoePrisma(end+1,1)=ureg.grp_mr_scanner_SIEMENS_Prisma_fit(idx);
 end
end

% check that scanner is set for all
if sum(~(reg.grp_siteTuePrisma|reg.grp_siteTueTrio|reg.grp_siteGoePrisma))>0
 reg.fname(~(reg.grp_scannerPrisma|reg.grp_scannerTrio))
 error('Some scan info missing, check input reg files')
end


% add quadraric age (demean first)
reg.age_square=(reg.age-mean(reg.age)).^2;

% and save
save(['unified_regressor' params.EEG.headmodel '.mat'],'-struct','reg')

% Update group selection with updated regressor, assuming that this is better / more correct
for i=1:length(processed_subjects)
 subject=processed_subjects{i};
 idx=find(strcmp(allGrps.basename,subject.id));
 if length(idx)==1
  processed_subjects{i}.group=allGrps.group{idx};
 else
  warning(['No single, unique match found for id=' subject.id '. This should not happen, investigate.'])
 end
end

% and save again
save(['processed_subjects' params.EEG.headmodel '.mat'],'processed_subjects')



%% Load regressor if needed
if ~exist('reg','var')
 load(['unified_regressor' params.EEG.headmodel '.mat'],'reg')
end




%%  make mean plots per group

% load the correct suma-surface for vizualization, using the correct one is
% cricical for subcortical nuclei!

%load(fullfile(pwd,'conf/suma-all-fsaverage-10.mat'))
load(fullfile(pwd,'conf/suma-all-nmriSUMA150_fs6-10.mat'))


% set the metrics we want to see as a cell array 
mtxts={'coh_img'};
groups_selection={'Cx','IGE|Pat'};
groups_labels={'Controls','IGE'};


for m=1:length(mtxts)
 % set the one to work with
 mtxt=mtxts{m};

 % set and make an output dir
 outdir=fullfile(pwd,'results','mean_maps',params.EEG.headmodel);
 if ~exist(outdir,'dir')
  mkdir(outdir)
 end

 % loop for all frequency bands
 for i=1:length(params.freqsNames)
  freq=params.freqsNames{i};
  
  
  % Now load the exported maps, if you used a different export path, do
  % modify this line, otherwise you will get a file not found error
  load(fullfile(pwd,'export',['all_subjects_' params.EEG.headmodel],[mtxt '_' freq '_abs_not_scaled_individual_' params.EEG.headmodel  '_N' num2str(length(processed_subjects))  '.mat']))

  % prepare the visualization config
  cfg=[];
  cfg.per_hemi=1;
  cfg.per_cortex=1;
  
  for g=1:length(groups_selection)
   % now select this group only, minimum age of 16+
   sel=zeros(length(processed_subjects),1);
   for ii=1:length(processed_subjects)
    sel(ii)=any(regexpi(processed_subjects{ii}.group,groups_selection{g}))&&reg.age(strcmp(reg.fname,processed_subjects{ii}.id))>15;
   end
  
   % now we make a mean map from our selection variable
   mean_map=nanmean(cat(3,scale_metrics{sel==1}),3); % get mean

   % now plot it
   % for the first (usually controls), we make the color scaling
   if ~isfield(cfg,'clim')
    cfg.clim=[0 prctile(mean_map,99)*1.2];
   end

   cfg.title=[mtxt ' - ' groups_labels{g} ' - '  freq ', N=' num2str(sum(sel))];
   cfg.output=fullfile(outdir,[mtxt '-' freq '-' groups_labels{g} '.png']);
   hFig=nmri_plot_surface_suma(suma_all,mean_map,cfg);
   close(hFig)
  end
 end
end

%% now do the PALM group analyis


% define the metrics, groups and matching principle to use

metrics={'coh_img','coh_real','power'}; % short standard names
metrics_txt={'Coherency (imaginary)','Coherency (real)','Power'};% longer verbose labels for plots
match_str={{'full','',0},{'1-1-age_sex','age+sex',1},{'2-1-age_sex','age+sex',2}}; % different matching strategies to use

% other strageties
%match_str={{'full','',0},{'2-1-sex','sex>age',2},{'3-1-sex','sex>age',3},{'1-1-sex','sex>age',1},{'2-1-age_sex','age+sex',2},{'3-1-age_sex','age+sex',3},{'1-1-age_sex','age+sex',1}};
%match_str={{'full','',0},{'2-1-sex','sex>age',2},{'1-1-sex','sex>age',1},{'2-1-age_sex','age+sex',2},{'1-1-age_sex','age+sex',1}};


groups={'Cx','IGE|Pat'};
groups_labels={'Controls','IGE'};



% whichs groups to compare, not 1st group will be matched (reduced) to the
% second one
groups_compare=[
 1,2; % Cx vs. IGE
 ];


% loop  over matching strategies
for si=1:length(match_str)
 % set the matching strategy
 strategy=match_str{si};
 
 % make an output dir
 outdir=fullfile(pwd,'results',['palm_' params.EEG.headmodel],strategy{1});
 if ~exist(outdir,'dir')
  mkdir(outdir)
 end
 
 % loop over metrics
 for mi=1:length(metrics)
  
  % make a cell array of selections, one per group requested
  gsel={};
  for gi=1:length(groups)
   grp=groups{gi};
   % filter and merge
   sel=ones(length(processed_subjects),1);
   for ii=1:length(processed_subjects)
    % this is main selection algorithm. You may need to adapt this to your
    % needs. Here, we only need the groups (taken from the .group field in
    % the subjects structs) and the age taken from reg.age. We limit to >16
    % years of age. Some subjects in this project were younger...
    sel(ii)=any(regexpi(processed_subjects{ii}.group,grp))&&reg.age(strcmp(reg.fname,processed_subjects{ii}.id))>15;
   end
   gsel{gi}=sel==1; %make logical and remember
   fprintf('Found N=%d subjects for %s\n',sum(sel),groups_labels{gi})
  end
  
  
  % now loop over the requested comparisons
  for gc=1:size(groups_compare,1)
   fprintf('Now comparing %s with %s...\n',groups_labels{groups_compare(gc,:)})
   all_sel=any([gsel{groups_compare(gc,:)}],2); %make the joined selections
   ref_sel=any([gsel{groups_compare(gc,1)}],2); % reference/controls only
   pat_sel=any([gsel{groups_compare(gc,2)}],2); % patients only
   % make a unified folder name
   grp=strjoin({groups_labels{groups_compare(gc,:)}},'-');
   % make the output dir
   this_out=fullfile(outdir,grp,metrics{mi});
   if ~exist(this_out,'dir')
    mkdir(this_out)
   end
   
   
   % extract matching variables, group A is all patients
   this_fnames=cellfun(@(x) x.id,processed_subjects(pat_sel),'UniformOutput',false);
   regA=[];
   regA.fname=this_fnames;
   regA.age=ureg.age(cellfun(@(x) find(strcmp(x,ureg.fname)),this_fnames));
   regA.grp_sex=ureg.grp_sex(cellfun(@(x) find(strcmp(x,ureg.fname)),this_fnames));

   % group B is controls/reference, will be matched to A
   this_fnames=cellfun(@(x) x.id,processed_subjects(ref_sel),'UniformOutput',false);
   regB=[];
   regB.fname=this_fnames;
   regB.age=ureg.age(cellfun(@(x) find(strcmp(x,ureg.fname)),this_fnames));
   regB.grp_sex=ureg.grp_sex(cellfun(@(x) find(strcmp(x,ureg.fname)),this_fnames));

   fprintf('Strategy=%s\n',strategy{1})
   
   % check if enough control cases for the matching strategy 
   if sum(ref_sel)<(sum(pat_sel)*strategy{3})
    fprintf('Not enough controls for %s, %s, %s. Skipping\n',strategy{1},grp,metrics{mi})
    continue
   end
   
   % check if we want a matching = reduction of cases in reference group
   if strategy{3}>0 
    % now make a matching, if not (yet) present
    if ~exist(fullfile(outdir,grp,'case_matching.mat'),'file')  
     [idx] = nf_reg_match_ratio(regA,regB, strategy{3},strategy{2});
     save(fullfile(outdir,grp,'case_matching.mat'),'idx')
    else
     load(fullfile(outdir,grp,'case_matching.mat'))
    end
   else
    % no matching requested, take as-is but make a distribution test  
    idx=1:length(regB.fname);    
   end


   % plot distribution of age and give some info   
   [h,p_age,ks2stat] = kstest2(regA.age,regB.age(idx));
   fprintf('Age:\n%0.2f+/-%0.2f (mean+/-std.dev) in patients, %0.2f+/-%0.2f (mean+/-std.dev) in reference\n',mean(regA.age),std(regA.age),mean(regB.age(idx)),std(regB.age(idx)))   
   fprintf('Kolmogorov-Smirnov test for age coming from the same distribution is:\n');
   fprintf('p: %0.4f\n\n',p_age);
   fprintf('Sex:\n%d females, %d males for reference\n%d females, %d males for matched group\n',sum(regA.grp_sex==1),sum(regA.grp_sex==0),sum(regB.grp_sex(idx)==1),sum(regB.grp_sex(idx)==0))
   fprintf('%0.1f%% females, %0.1f%% males for reference\n%0.1f%% females, %0.1f%% males for matched group\n',sum(regA.grp_sex==1)*100/length(regA.grp_sex),sum(regA.grp_sex==0)*100/length(regA.grp_sex),...
   sum(regB.grp_sex(idx)==1)*100/length(regB.grp_sex(idx)),sum(regB.grp_sex(idx)==0)*100/length(regB.grp_sex(idx)))
   [tbl,ch2,p_sex] = crosstab([regA.grp_sex ; regB.grp_sex(idx)], [ones(length(regA.grp_sex),1); zeros(length(idx),1)]);
   fprintf('Chi-Square test for coming from the same distribution is:\n');
   fprintf('p: %0.4f\n\n',p_sex);
   
   
   hFig=figure;
   hV=violinplot([regB.age(idx); regA.age],[zeros(size(idx')); ones(size(regA.age)) ] ,'ViolinAlpha',0.4,'BoxColor',[0.2 0.2 0.2]);
   title([grp ' ' strategy{1} sprintf(', p(age)=%0.2f, p(sex)=%0.2f',p_age,p_sex)])
   xticklabels({groups_labels{groups_compare(gc,:)}})
   export_fig(hFig,fullfile(outdir,grp,'case_matching.png'),'-nocrop')
   close(hFig)
   
   
   
   % now make design and files
   ref_idx=find(ref_sel); % get the index of the reference from the original list
   controls=processed_subjects(ref_idx(idx));
   patients=processed_subjects(pat_sel);   
   this_subjects=[controls; patients];
  
   % now make design
   design_mat=zeros(length(this_subjects),2);
   % set controls / reference
   design_mat(1:length(idx),1)=1;
   % set patients
   design_mat(length(idx)+1:end,2)=1;


   % make the conf for PALM
   cfg=[];
   cfg.surface=fullfile(pwd,'conf/suma-all-fsaverage-10.gii');
   cfg.plot_surface=fullfile(pwd,'conf/suma-all-fsaverage-10.mat');
   cfg.data={};
   cfg.data_labels={};
   cfg.viz_data_back={'#FFCCBB','#FFBBFF','#DDDDFF','#AAFFFF','#BBFFDD','#DDFFBB'};
   cfg.subject_ids=cellfun(@(x) x.id,this_subjects,'UniformOutput',false); % just give subject ids
   cfg.output=this_out;
   cfg.reg=reg;
   cfg.reg_use={'age','grp_sex','grp_scannerPrismaTue','grp_scannerPrismaGoe'};

   cfg.design=design_mat;
   cfg.design_col_labels={groups_labels{groups_compare(gc,:)}};
  
   cfg.contrast=[1 -1;-1 1];
   cfg.contrast_colormap={'cool';'autumn'};
   cfg.contrast_labels={[groups_labels{groups_compare(gc,1)} ' > ' groups_labels{groups_compare(gc,2)}];[groups_labels{groups_compare(gc,1)} ' < ' groups_labels{groups_compare(gc,2)}]};
   
   % set title
   cfg.title=metrics_txt{mi};

   % check if vertex PALM is done
   if ~exist(fullfile(cfg.output,['concat_fig_uncp_' cfg.design_col_labels{1} '_less_' cfg.design_col_labels{end} '.png']),'file')

    % now loop the frequencies and write mgh files
    for i=1:length(params.freqsNames)
     freq=params.freqsNames{i};
     if ~exist(fullfile(this_out,[ freq '_merged.mgh']),'file') || ~exist(fullfile(this_out,[ freq '_global.csv']),'file')
      fprintf('Writing vertex-based MGH: %s, %s\n',metrics{mi},freq)
      load(fullfile(pwd,'export',['all_subjects_' params.EEG.headmodel],[metrics{mi} '_' freq '_abs_not_scaled_individual_' params.EEG.headmodel  '_N' num2str(length(processed_subjects))  '.mat']))
      nmri_write_mgh(fullfile(this_out,[ freq '_merged.mgh']),eye(4),[scale_metrics(ref_idx(idx));scale_metrics(pat_sel)])
      fprintf('Writing global values CSV: %s, %s\n',metrics{mi},freq)
      csvwrite(fullfile(this_out,[ freq '_global.csv']),[cellfun(@nanmean,scale_metrics(ref_idx(idx)));cellfun(@nanmean,scale_metrics(pat_sel))])
     end
     cfg.data=[cfg.data {fullfile(this_out,[ freq '_merged.mgh'])}];
     cfg.data_labels=[cfg.data_labels {[freq ' (' num2str(params.freqs(i)) ' +/- ' num2str(params.tapsmofrq(i)) 'Hz)']}];
    end

    % make the mask
    if ~exist(fullfile(this_out,'all_suma_msk.mgh'),'file')
     fprintf('Writing Mask\n')
     allC=load(fullfile(pwd,'export',['all_subjects_' params.EEG.headmodel],'all_suma_msk.mat'));
     save_mgh(allC.all_msk,fullfile(this_out,'all_suma_msk.mgh'),eye(4));
    end
    cfg.mask=fullfile(this_out,'all_suma_msk.mgh');

    % set permutations
    cfg.permut=5000;
    cfg.vis_fwe=[1.3 2.5];
    cfg.vis_uncp=[1.3 3.0];
    % check if running
    if ~exist(fullfile(pwd,'SGE_calls',['palm-' grp '-' strategy{1} '-' metrics{mi}],'queue_lock'),'file')
     nmri_qsub(struct('compile',1,'overwrite',0,'title',['palm-' grp '-' strategy{1} '-' metrics{mi}]),'nmri_palm_run',cfg)
    else
     fprintf('PALM job for %s seems to be still running -- skipped\n',['palm-' grp '-' strategy{1} '-' metrics{mi}])
    end
    %nmri_palm_run(cfg); 
   else
    fprintf('%s, %s, %s is done already\n',strategy{1},grp,metrics{mi})
   end

   
  
   % now deal with global results

   % modify the config accordingly

   cfg.surface=[];
   cfg.plot_surface=[];
   cfg.data={};
   cfg.data_labels={};
   if strcmpi(metrics{mi},'power')
    cfg.global_vis_log=1;
   else
    cfg.global_vis_log=0;
   end
   
   % check if global PALM is done
   if ~exist(fullfile(cfg.output,['Global_fig_Violin.png']),'file')
    % now loop the frequencies and make .csv files
    for i=1:length(params.freqsNames)
     freq=params.freqsNames{i};
     cfg.data=[cfg.data {fullfile(this_out,[ freq '_global.csv'])}];
     cfg.data_labels=[cfg.data_labels {[freq ' (' num2str(params.freqs(i)) ' +/- ' num2str(params.tapsmofrq(i)) 'Hz)']}];
    end

    % set permutations
    cfg.permut=5000;

    % this is fast...
    % once with confounds removed
    cfg.global_viz_remove=cfg.reg_use;
    if strcmpi(metrics{mi},'power')
     cfg.global_vis_log=1;
     cfg.global_vis_lock_scale=0;
    else
     cfg.global_vis_log=0;
     cfg.global_vis_lock_scale=1;
    end
    nmri_palm_run(cfg); 

    % and once without
    if strcmpi(metrics{mi},'power')
     cfg.global_vis_log=1;
     cfg.global_vis_lock_scale=1;
    else
     cfg.global_vis_log=0;
     cfg.global_vis_lock_scale=1;
    end
    cfg.global_viz_remove={};
    nmri_palm_run(cfg); 

   else
    fprintf('%s, %s, %s is done already\n',strategy{1},grp,metrics{mi})
   end
  end
 end
end

%% now collect the main output into a folder

metrics={'coh_img','power'}; % short standard names
strategies={'full','1-1-age_sex','2-1-age_sex'}; % different matching strategies to use
groups={'Controls-IGE','Controls-SCN1A','Controls-STX1B','IGE-SCN1A','IGE-STX1B'};
outdir=fullfile(pwd,'results',['palm_' params.EEG.headmodel],'compiled_results');
if ~exist(outdir,'dir')
 mkdir(outdir)
end


for si=1:length(strategies)
 for gi=1:length(groups)
  thisoutdir=fullfile(outdir,strategies{si});
  if ~exist(thisoutdir,'dir')
   mkdir(thisoutdir)
  end
  if exist(fullfile(pwd,'results',['palm_' params.EEG.headmodel],strategies{si},groups{gi},'case_matching.png'),'file')
   copyfile(fullfile(pwd,'results',['palm_' params.EEG.headmodel],strategies{si},groups{gi},'case_matching.png'),fullfile(thisoutdir,[groups{gi} '_case_matching_.png']))
   for mi=1:length(metrics)
    palmdir=fullfile(pwd,'results',['palm_' params.EEG.headmodel],strategies{si},groups{gi},metrics{mi});
    files=dir(fullfile(palmdir,'concat_fig_fwep*.png'));
    for fi=1:length(files)
     copyfile(fullfile(palmdir,files(fi).name),fullfile(thisoutdir,[groups{gi} '_' metrics{mi} '_' strrep(files(fi).name,'concat_fig_','')]))
    end
    files=dir(fullfile(palmdir,'concat_fig_unc*.png'));
    for fi=1:length(files)
     copyfile(fullfile(palmdir,files(fi).name),fullfile(thisoutdir,[groups{gi} '_' metrics{mi} '_' strrep(files(fi).name,'concat_fig_','')]))
    end 
    copyfile(fullfile(palmdir,'Global_fig_Violin.png'),fullfile(thisoutdir,[groups{gi} '_' metrics{mi} '_Global.png']))
    copyfile(fullfile(palmdir,'Global_fig_Violin_confounds_removed.png'),fullfile(thisoutdir,[groups{gi} '_' metrics{mi} '_Global_adjusted.png']))
   end
  end
 end
end




%% now do a simple two-group analyis

% define the metrics, groups and matching principle to use

metrics={'coh_img','coh_real','power'}; % short standard names
metrics_txt={'Coherency (imaginary)','Coherency (real)','Power'};% longer verbose labels for plots

groups={'Cx','IGE'};
groups_labels={'Controls','GGE'};


% make an output dir
outdir=fullfile(pwd,'results',['simple_palm_' params.EEG.headmodel]);
if ~exist(outdir,'dir')
 mkdir(outdir)
end
 
% loop over metrics
for mi=1:length(metrics)
  
 % make a cell array of selections, one per group requested
 gsel={};
 for gi=1:length(groups)
  grp=groups{gi};
  % filter and merge
  sel=ones(length(processed_subjects),1);
  for ii=1:length(processed_subjects)
   % this is main selection algorithm. You may need to adapt this to your
   % needs. Here, we only need the groups (taken from the .group field in
   % the subjects structs) and the age taken from reg.age. We limit to >16
   % years of age. Some subjects in this project were younger...
   sel(ii)=any(regexpi(processed_subjects{ii}.group,grp))&&reg.age(strcmp(reg.fname,processed_subjects{ii}.id))>15;
  end
  gsel{gi}=sel==1; %make logical and remember
  fprintf('Found N=%d subjects for %s\n',sum(sel),groups_labels{gi})
 end
  
  
 fprintf('Now comparing %s with %s...\n',groups_labels{1},groups_labels{2})
 all_sel=any([gsel{:}],2); %make the joined selections
 ref_sel=gsel{1}; % reference/controls only
 pat_sel=gsel{2}; % patients only
 % make a unified folder name
 grp=strjoin({groups_labels{:}},'-');
 % make the output dir
 this_out=fullfile(outdir,grp,metrics{mi});
 if ~exist(this_out,'dir')
  mkdir(this_out)
 end
   
 % now make design and files
 controls=processed_subjects(ref_sel);
 patients=processed_subjects(pat_sel);   
 this_subjects=[controls; patients];
  
 % now make design
 design_mat=zeros(length(this_subjects),2);
 % set controls / reference
 design_mat(1:length(controls),1)=1;
 % set patients
 design_mat(length(controls)+1:end,2)=1;


 % make the conf for PALM
 cfg=[];
 cfg.surface=fullfile(pwd,'conf/suma-all-fsaverage-10.gii');
 cfg.plot_surface=fullfile(pwd,'conf/suma-all-fsaverage-10.mat');
 cfg.data={};
 cfg.data_labels={};
 cfg.viz_data_back={'#FFCCBB','#FFBBFF','#DDDDFF','#AAFFFF','#BBFFDD','#DDFFBB'};
 cfg.subject_ids=cellfun(@(x) x.id,this_subjects,'UniformOutput',false); % just give subject ids
 cfg.output=this_out;
 cfg.reg=reg;
 cfg.reg_use={'age','grp_sex'};

 cfg.design=design_mat;
 cfg.design_col_labels=groups_labels;

 cfg.contrast=[1 -1;-1 1];
 cfg.contrast_colormap={'cool';'autumn'};
 cfg.contrast_labels={[groups_labels{1} ' > ' groups_labels{2}];[groups_labels{1} ' < ' groups_labels{2}]};

 % set title
 cfg.title=metrics_txt{mi};
 % check if vertex PALM is done
 if ~exist(fullfile(cfg.output,['concat_fig_uncp_' cfg.design_col_labels{1} '_less_' cfg.design_col_labels{end} '.png']),'file')

  % now loop the frequencies and write mgh files
  for i=1:length(params.freqsNames)
   freq=params.freqsNames{i};
   if ~exist(fullfile(this_out,[ freq '_merged.mgh']),'file') || ~exist(fullfile(this_out,[ freq '_global.csv']),'file')
    fprintf('Writing vertex-based MGH: %s, %s\n',metrics{mi},freq)
    load(fullfile(pwd,'export',['all_subjects_' params.EEG.headmodel],[metrics{mi} '_' freq '_abs_not_scaled_individual_' params.EEG.headmodel  '_N' num2str(length(processed_subjects))  '.mat']))
    nmri_write_mgh(fullfile(this_out,[ freq '_merged.mgh']),eye(4),[scale_metrics(ref_sel);scale_metrics(pat_sel)])
    fprintf('Writing global values CSV: %s, %s\n',metrics{mi},freq)
    csvwrite(fullfile(this_out,[ freq '_global.csv']),[cellfun(@nanmean,scale_metrics(ref_sel));cellfun(@nanmean,scale_metrics(pat_sel))])
   end
   cfg.data=[cfg.data {fullfile(this_out,[ freq '_merged.mgh'])}];
   cfg.data_labels=[cfg.data_labels {[freq ' (' num2str(params.freqs(i)) ' +/- ' num2str(params.tapsmofrq(i)) 'Hz)']}];
  end

  % make the mask
  if ~exist(fullfile(this_out,'all_suma_msk.mgh'),'file')
   fprintf('Writing Mask\n')
   allC=load(fullfile(pwd,'export',['all_subjects_' params.EEG.headmodel],'all_suma_msk.mat'));
   save_mgh(allC.all_msk,fullfile(this_out,'all_suma_msk.mgh'),eye(4));
  end
  cfg.mask=fullfile(this_out,'all_suma_msk.mgh');

  % set permutations
  cfg.permut=5000;
  cfg.vis_fwe=[1.3 2.5];
  cfg.vis_uncp=[1.3 3.0];
  
  % run PALM locally
  nmri_palm_run(cfg); 
 else
  fprintf('%s, %s is done already\n',grp,metrics{mi})
 end



 % now deal with global results

 % modify the config accordingly

 cfg.surface=[];
 cfg.plot_surface=[];
 cfg.data={};
 cfg.data_labels={};
 if strcmpi(metrics{mi},'power')
  cfg.global_vis_log=1;
 else
  cfg.global_vis_log=0;
 end

 % check if global PALM is done
 if ~exist(fullfile(cfg.output,['Global_fig_Violin.png']),'file')
  % now loop the frequencies and make .csv files
  for i=1:length(params.freqsNames)
   freq=params.freqsNames{i};
   cfg.data=[cfg.data {fullfile(this_out,[ freq '_global.csv'])}];
   cfg.data_labels=[cfg.data_labels {[freq ' (' num2str(params.freqs(i)) ' +/- ' num2str(params.tapsmofrq(i)) 'Hz)']}];
  end

  % set permutations
  cfg.permut=5000;

  % this is fast...
  % once with confounds removed
  cfg.global_viz_remove=cfg.reg_use;
  if strcmpi(metrics{mi},'power')
   cfg.global_vis_log=1;
   cfg.global_vis_lock_scale=0;
  else
   cfg.global_vis_log=0;
   cfg.global_vis_lock_scale=1;
  end
  nmri_palm_run(cfg); 

  % and once without
  if strcmpi(metrics{mi},'power')
   cfg.global_vis_log=1;
   cfg.global_vis_lock_scale=0;
  else
   cfg.global_vis_log=0;
   cfg.global_vis_lock_scale=1;
  end
  cfg.global_viz_remove={};
  nmri_palm_run(cfg); 

 else
  fprintf('Global %s, %s, %s is done already\n',grp,metrics{mi})
 end
end

  