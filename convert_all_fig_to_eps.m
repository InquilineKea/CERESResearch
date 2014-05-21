function [  ] = convert_all_fig_to_eps(  )
%CONVERT_ALL_FIG_TO_PNG will convert all the fig files in a folder to png
%files
% Searching using wildcards
fig_files=dir(fullfile(pwd,'\*.fig'));
for i=1:length(fig_files)
  open(fig_files(i).name);
  print('-depsc2','-loose', '-r300', strcat(fig_files(i).name(1:end-4),'.eps'));
  close all;
end
end
