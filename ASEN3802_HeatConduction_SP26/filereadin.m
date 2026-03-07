clear;
clc;
close all;

a=dir('*mA');

for i=1:length(a)
    load(a(i).name)
% how to get voltage and amperage from file names?
% - options include strsplit, regex, etc.
% ultimately, we need to use the format of each file name
% 'material'_'volts'V_'amps'mA
b = strsplit(a(i).name,'_'); % gives a cell array (b) that is 1x3
% {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'
v = strsplit(b{2},'V'); % volts are always in the second portion
ampval= strsplit(b{3},'mA'); % amps are always in the third portion
volts(i) = str2num(v{1}); % convert string to number (vector)
amps(i) = str2num(ampval{1});
end
