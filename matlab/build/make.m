% a very simple make script for cosmic

%% add the current path and the data path
wd = pwd;
addpath([wd '/../']);
addpath([wd '/../../data/']);

%% build the executable
mcc -m -v -R -nodisplay cosmic.m
