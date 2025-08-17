%% script to generate channel reductions for the EGI 257 layout

%% load layouts
egi=ft_prepare_layout(struct('layout','conf/layouts/GSN-HydroCel-257-layout.mat'));
eeg1020=ft_prepare_layout(struct('layout','EEG1020'));
 
eeg1010=ft_prepare_layout(struct('layout','EEG1010'));
 

% now plot
figure
hold on
scatter(eeg1010.pos(:,1),eeg1010.pos(:,2),'MarkerFaceColor','b')
text(eeg1010.pos(:,1)+0.01,eeg1010.pos(:,2), eeg1010.label)

% scatter(egi.pos(:,1),egi.pos(:,2),'MarkerFaceColor','r')
% text(egi.pos(:,1)+0.01,egi.pos(:,2), egi.label, 'Color','r')

% you will notice...this does not match. We shall use another approach


%% Read in the Excel CSV 
% from https://www.egi.com/knowledge-center/item/62-relating-the-hcgsn-sensor-positions-to-the-10-10-international-electrode-placement-system
fileID = fopen('conf/EGI-EEG1010-classical.csv','r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', ';,', 'TextType', 'string', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% now create new channel sets based in the full EGI setup

egi_channels=ft_channelselection({'eeg','eeg1020'},egi.label); %get all EEG + Cz

% match with 10-20
new_lay=[];
new_lay.pos=[];
new_lay.width=[];
new_lay.height=[];
new_lay.label={};

for i=1:length(dataArray{1})
 curr=strtrim(dataArray{1}{i});
 egicurr=strtrim(dataArray{2}{i});
 if ~strcmpi(curr,'Cz')
  idx=find(strcmpi(egi.label,['E' egicurr]));
 else
  % deal with Cz differently, us the slightly posterior E
  idx=find(strcmpi(egi.label,'E81'));
 end
 if idx>0
  new_lay.pos(end+1,1:2)=egi.pos(idx,1:2);
  new_lay.width(end+1,1)=egi.width(idx,1);
  new_lay.height(end+1,1)=egi.height(idx,1); 
  new_lay.label{end+1,1}=egi.label{idx,1};
 end
end
% add the outline/mask fields
new_lay.outline=egi.outline;
new_lay.mask=egi.mask;

% now plot
figure
hold on
scatter(new_lay.pos(:,1),new_lay.pos(:,2),'MarkerFaceColor','b')
text(new_lay.pos(:,1)+0.01,new_lay.pos(:,2), new_lay.label)
text(new_lay.pos(:,1)+0.01,new_lay.pos(:,2)+0.05, dataArray{1})

% now save it
layout=new_lay;
save('conf/layouts/GSN-HydroCel-257-layout-1020-classical_test.mat','layout')


%% Read in the Excel CSV 
% from https://www.egi.com/knowledge-center/item/62-relating-the-hcgsn-sensor-positions-to-the-10-10-international-electrode-placement-system
fileID = fopen('conf/EGI-EEG1010-extended.csv','r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', ';,', 'TextType', 'string', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% now create new channel sets based in the full EGI setup

egi_channels=ft_channelselection({'eeg','eeg1020'},egi.label); %get all EEG + Cz

% match with 10-20
new_lay=[];
new_lay.pos=[];
new_lay.width=[];
new_lay.height=[];
new_lay.label={};

for i=1:length(dataArray{1})
 curr=strtrim(dataArray{1}{i});
 egicurr=strtrim(dataArray{2}{i});
 if ~strcmpi(curr,'Cz')
  idx=find(strcmpi(egi.label,['E' egicurr]));
 else
  % deal with Cz differently, us the slightly posterior E
  idx=find(strcmpi(egi.label,'E81'));
 end
 if idx>0
  new_lay.pos(end+1,1:2)=egi.pos(idx,1:2);
  new_lay.width(end+1,1)=egi.width(idx,1);
  new_lay.height(end+1,1)=egi.height(idx,1); 
  new_lay.label{end+1,1}=egi.label{idx,1};
 end
end
% add the outline/mask fields
new_lay.outline=egi.outline;
new_lay.mask=egi.mask;


% now plot
figure
hold on
scatter(new_lay.pos(:,1),new_lay.pos(:,2),'MarkerFaceColor','b')
text(new_lay.pos(:,1)+0.01,new_lay.pos(:,2), new_lay.label)
text(new_lay.pos(:,1)+0.01,new_lay.pos(:,2)+0.05, dataArray{1})


% now save it
layout=new_lay;

egi_1020=new_lay;

save('conf/layouts/GSN-HydroCel-257-layout-1020-extended_test.mat','layout')

%% Now generate reductions
fs=[192,128,64,48,32];

% now create new channel sets based in the full EGI setup
egi_channels=ft_channelselection({'eeg'},egi.label); %get only E-type channels

cfg=[];
cfg.layout=egi;
ft_layoutplot(cfg);

 new_lay=[];
 [~,fil]=intersect(egi.label,egi_channels);
 new_lay.pos=egi.pos(fil,1:2);
 new_lay.width=egi.width(fil,1);
 new_lay.height=egi.height(fil,1);
 new_lay.label=egi.label(fil,1);


neigh=ft_prepare_neighbours(struct('method','triangulation','layout',new_lay))

for f=1:length(fs)
 % how many channels to remove
 remCh=length(egi_channels)-fs(f);
 
 new_lay=[];
 [~,fil]=intersect(egi.label,egi_channels);
 new_lay.pos=egi.pos(fil,1:2);
 new_lay.width=egi.width(fil,1);
 new_lay.height=egi.height(fil,1);
 new_lay.label=egi.label(fil,1);

 % now check which channel to remove keeping minimum distances, but
 % preserve 10-20 channels

 remCount=0;
 while remCount<remCh
  
  % check number of neighbours left
  nn=zeros(length(new_lay.label),1);
  for n=1:length(new_lay.label)
   pos_neigh=neigh(find(strcmp({neigh(:).label},new_lay.label(n)))).neighblabel;
   nn(n)=length(intersect(pos_neigh,new_lay.label));
  end
  
  % now protect 10-20
  for ci=1:length(new_lay.label)
   if any(strcmpi(new_lay.label{ci},egi_1020.label))
    nn(ci)=-1; % assume negative neighbours
   end
  end
    
  % now find out the one with most neighbours
  [nm,ni]=sort(nn);
  
  sel_among=ni(nm==nm(end));
  
  % all to all distances 
  if length(sel_among)>1
   % we need to choose
   pd=pdist2(new_lay.pos(sel_among,:),new_lay.pos);
   pd(pd==0)=NaN;
 
   [m,idx]=min(pd(:));
   [x,~]=ind2sub(size(pd),idx); % take only the 1st index
   rmidx=sel_among(x);
  else
   rmidx=sel_among;
  end
ö
  
  % also take the corresponding other side channel
  taken_pos(1)=-taken_pos(1); % flip x
  pd=pdist2(taken_pos,new_lay.pos);
  [m,rmidx]=min(pd(:));
  
  % and take if clos match
  if (m<0.01 && remCount<remCh)
   new_lay.pos(rmidx,:)=[];
   new_lay.width(rmidx,:)=[];
   new_lay.height(rmidx,:)=[];
   new_lay.label(rmidx,:)=[];
   remCount=remCount+1;
  end
  
 end

 % check
 if length(new_lay.label)~=fs(f)
  error('Could not get the requested N')
 end
 
 % and re-add scale and comment
 new_lay.pos(end+1:end+2,1:2)=egi.pos(end-1:end,1:2);
 new_lay.width(end+1:end+2,1)=egi.width(end-1:end,1);
 new_lay.height(end+1:end+2,1)=egi.height(end-1:end,1); 
 new_lay.label(end+1:end+2,1)=egi.label(end-1:end,1);
 
 % plot new
 %cfg=[];
 %cfg.layout=new_lay;
 %ft_layoutplot(cfg);
 
  figure
  hold on
  scatter(egi.pos(:,1),egi.pos(:,2),'MarkerFaceColor','g')
  scatter(new_lay.pos(:,1),new_lay.pos(:,2),'MarkerFaceColor','b')
  text(new_lay.pos(:,1)+0.01,new_lay.pos(:,2), new_lay.label)
  scatter(egi_1020.pos(:,1),egi_1020.pos(:,2),'MarkerFaceColor','k')
 
 
 % now save it
 layout=new_lay;
 save(['conf/layouts/GSN-HydroCel-257-layout-reduced-' num2str(length(layout.label)-2) '.mat'],'layout')

end

%% Make special reductions-- Gö Monitoring
load('/home/uni10/nmri/projects/nfocke/dev_nmri/P.IWTS99A/processed/clean_P.IWTS99A_eeg_EEG_rest_Routine.mat', 'data')


% Read in the Excel CSV 
% from https://www.egi.com/knowledge-center/item/62-relating-the-hcgsn-sensor-positions-to-the-10-10-international-electrode-placement-system
fileID = fopen('conf/EGI-EEG1010.csv','r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', ';,', 'TextType', 'string', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

% Close the text file.
fclose(fileID);

% now create new channel sets based in the full EGI setup
egi_channels=ft_channelselection({'eeg','eeg1020'},egi.label); %get all EEG + Cz

% match with the reference (from loaded data)
new_lay=[];
new_lay.pos=[];
new_lay.width=[];
new_lay.height=[];
new_lay.label={};
new_lay.orig_label={};

all_egi=strtrim(dataArray{1});

for i=1:length(data.label)
 curr=data.label{i};
 cidx=find(strcmpi(curr,all_egi));
 if ~isempty(cidx)
  egicurr=strtrim(dataArray{2}{cidx});
  if ~strcmpi(curr,'Cz')
   idx=find(strcmpi(egi.label,['E' egicurr]));
  else
   % deal with Cz differently, us the slightly posterior E
   idx=find(strcmpi(egi.label,'E81'));
  end
  if idx>0
   new_lay.pos(end+1,1:2)=egi.pos(idx,1:2);
   new_lay.width(end+1,1)=egi.width(idx,1);
   new_lay.height(end+1,1)=egi.height(idx,1); 
   new_lay.label{end+1,1}=egi.label{idx,1};
   new_lay.orig_label{end+1,1}=curr;
  end
 end
end
% add the outline/mask fields
new_lay.outline=egi.outline;
new_lay.mask=egi.mask;



% now plot
figure
hold on
scatter(new_lay.pos(:,1),new_lay.pos(:,2),'MarkerFaceColor','b')
text(new_lay.pos(:,1)+0.01,new_lay.pos(:,2), new_lay.label)


% now save it
layout=new_lay;
save('conf/layouts/GSN-HydroCel-257-layout-1020-GoeMonitoring.mat','layout')

