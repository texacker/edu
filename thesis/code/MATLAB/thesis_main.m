echo off; clear; clc; format short e;
colordef white;
%whitebg('w');
thesis_data;
thesis_init;
thesis_solve;
disp(['Total time:' num2str(cputime - start_time) ' seconds']);