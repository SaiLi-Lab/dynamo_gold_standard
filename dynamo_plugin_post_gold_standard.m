% INPUT:
%
%
%  file_card:      Name of the card file for this iteration.
%
%                  IMPORTANT
%                  When you invoke the plugin, you do NOT pass this
%                  variable.
%                  This will be automatically provided by Dynamo in runtime.
%
%
%  After this I/O related file, those are the variables explicitely
%  passed when invoking the plugin:
%
%  'threshold'   : for resolution estimation.
%
%  'pushback'    : to effect the change of lowpass for next iteration
%
%        
%  OUTPUT
%
%  As a rule, do not use output for the command line.
%  All the I/O should run through database operations.
%
%  INVOKING THE PLUGIN
%
%  This plugin can be invoked inside a Dynamo project using the project parameter
%  'plugin_post_order'
%
%  For instance, to include it after the first round of  a project 
%  represented by a virtual project vpr, you would use:
%
%  vpr=dynamo_vpr_modify(vpr,'plugin_post_order_r1','dynamo_plugin_multireference_tutorial');
%  
%  % This just changes the virtual project: you must unfold the virtual
%  % project before running it:
%  
%  dynamo_vpr_unfold(vpr);
%
%
%  Alternatively, you can use the GUI dynamo_project_manager, pressing
%  the [edit] field in the plugin area (on the right).
%
%
%  Note in both cases that the plugin will not be still active: you need to
%  switch on the parameter 'plugin_post' of the corresponding round to make
%  clear that the execution pipeline must visit the plugin.
%
%  Author: Sai Li
%  Email: sai @ tsinghua.edu.cn
%  Address: School of Life Science, Tsinghua University
%  https://www.lisailab.com
%
%  Copyright (c) 2020
%
%  Citation: if you found this function useful for your research, please
%  cite the following paper to ackownlege the author:
%  "Molecular architecture of the SARS-CoV-2 virus", Cell, 2020
%  DOI: 10.1016/j.cell.2020.09.018

 
function dynamo_plugin_post_gold_standard_sai(file_card,input_string)

tic

%%%%%%% String input parser %%%%%%%%%
inputs = strsplit(input_string, ' '); %split string around the space delimeter
%name the individual inputs
fscmask = inputs{1};
pushback = str2num(inputs{2}); %for numbers, convert the string to a number.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% reads the card text file, to parse their contents into a easily
% readable MATLAB structure

scard = dynamo_read(file_card);

% Now the plugin locates:
%   1 from which iteration it is being executed:
ite = scard.database_ite;

%   2 from which project it is being invoked
name_project=scard.name_project;

% particles size
dim=scard.dim;

% what kind of symmetry is applied to the structure
sym=scard.sym;

% We read the corresponding virtual project
vpr=dynamo_vpr_load(name_project);

% Optional: defines an output channel for messages
o = pparse.output();
o.leading = 'Gold_standard'; % all messages will be preceded with the mark
                      % 'Gold_standard'

% allows echoing our displayed text into the log file of the project
logfile =[name_project,filesep,'log.txt'];                      
o.intofile   = logfile; 

% informs the user that this plugin has been enterered
disp(repmat('-',[1,60]));
disp(' ');
disp('   [Gold Standard: Dynamic bandpass] entering post processing plugin.');
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        PLUGIN OPERATIONS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The goal of the plugin is to identify for each particle which is the
% reference channel that produces the best correlation 

% A small test: do not operate this plugin if the project is not
% multireference:
if scard.nref==1
    disp('Only one reference available. Desisting [Gold_standard] plugin');
    return;
end

if scard.nref>2
    o.echo('Too many references: this plugin will work for simple refinements. Desisting [Gold_standard] plugin');
    o.echo(' ... but you can always write your own plugin!');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We first locate all the availabe symmetrized average files
average_ref_files = dynamo_database_locate('average',vpr,'ite',ite,'ref','*','check','disk');

% Note that we use the 'check' / 'disk' option. This ensures that the function will produce 
% an empty list if the database items of category 'refined_average' are not available
% We can make a small control on that: it is a good programming practice;

if isempty(average_ref_files)
   
    % informs the usr that something went wrong:
    disp('Sorry: no "average" files were available for this iteration');
    
    % and errors out to prevent further damage
    error('Aborting dynamo_plugin_multireference_tutorial');
end

% we can now read the found files:
average_ref = dynamo_read(average_ref_files);

% and make certain that they are listed inside MATLAB in cell array format:
average_ref = cell_checkin(average_ref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Alignment of reference 2 into reference 1
%  
%  Using dynamo_align
%  ------------------
%  reference  1 is treated as the "template"
%  reference  2 is treated as the "particle"

current_round = dynamo_intertwin_ite2round(vpr,ite);
m = dynamo_read(fscmask);

% first mask the unsymmetrized average
average_ref1_masked = average_ref{1}.*m;
average_ref2_masked = average_ref{2}.*m;

% then symmetrized them
template_masked  = dynamo_sym(average_ref1_masked,sym);
particle_masked = dynamo_sym(average_ref2_masked,sym);

% align the masked and symmetrized two maps
sal = dynamo_align(particle_masked,template_masked,'settings_from',vpr.name_project,'settings_round',current_round);

% First re-algin the masked & symmetrized particle for fsc
average_ref_2_aligned_to_ref_1_masked = dynamo_rigid_apply(particle_masked,sal.Tp);

% Then re-algin the unsymmetrized average_ref2 to ref1 to save it
average_ref_2_aligned_to_ref_1_unsym = dynamo_rigid_apply(average_ref{2},sal.Tp);
average_ref_2_aligned_to_ref_1_file = dynamo_database_locate('average',vpr,'ref',2,'ite',ite,'append','aligned_to_ref1');
dynamo_write(average_ref_2_aligned_to_ref_1_unsym,average_ref_2_aligned_to_ref_1_file);

%%%%%%%%%%%%%%%%%%%% FSC estimation %%%%%%%%%%%%%%%%%%
threshold = 0.143;

% we compute the fsc
gs_fsc = dynamo_fsc(template_masked,average_ref_2_aligned_to_ref_1_masked);

% Optional: We can store it in some location of the project
% For instance we can just place it with the normal eo_fsc curves but prepending and 'gs' mark, i.e.
location_eo = dynamo_database_locate('eo_fsc',name_project,'ite',ite,'ref',1);
location_gs = dynamo_filename(location_eo,'prepend','gs');
dynamo_write(gs_fsc.fsc,location_gs);
o.echo('Gold standard FSC for resolution estimation saved in file:');
o.echo(location_gs);

% then we check the attained resolution with the given criterion
[fourier_pixels,check] = dynamo_fsc_resolution(gs_fsc,'threshold',threshold);

resolution_report = dynamo_database_locate('bandpass_resolution',vpr,'prepend','gold_standard');
fid = fopen(resolution_report,'a+');
fprintf(fid,'iteration_%d\n %d\n',ite,fourier_pixels);
fclose(fid);

if ~check.ok
   o.echo('Problem found when estimating the resolution in iteration %d. ',ite); 
   return 
end

% stablishes a new lowpass for new iteration
new_low_pass = fourier_pixels-pushback;

o.echo('Resolution estimated for iteration %d: %5.2f (threshold %5.2f)',{ite,fourier_pixels,threshold});
o.echo('Low pass computed for next iteration %d: %5.2f (pushback %5.2f)',{ite+1,new_low_pass,pushback});

% runs appropriate checks
if new_low_pass<2
   o.echo('New low pass is too low! You are killing all the frequencies!!');
   error('aborting!')
end

% locates the round that corresponds to the next iteration
next_round = dynamo_intertwin_ite2round(vpr,ite+1);

if isempty(next_round)
    o.echo('Could not identify round number for next iteration (%d)',ite+1);
   return 
end

% changes the project settings for next iteration
check_put = dynamo_vpr_put(name_project,'disk','inround',next_round,'low',new_low_pass);
if ~check_put.ok
    o.echo('Could not update the project parameters (to use a new low)');
   return 
end

% a partial unfolds for the next iteration
[dummy1,dummy2,check_unfold] = dynamo_vpr_unfold(name_project,'refresh_cards',ite+1);

if ~check_unfold.ok
    o.echo('Could not unfold after updating the project parameters (to use a new lowpass)');
   return 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elapsedTime = toc;
disp(' ');
fprintf(1,'Resolution successfully calculated and updated to the next iteration, the process took %f seconds.',elapsedTime);
disp(' ');
disp(repmat('-',[1,60]));
