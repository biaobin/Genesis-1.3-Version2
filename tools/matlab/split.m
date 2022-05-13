function [S]=split(str)
%split function
str= strtrim(str);% remove space head and tail
S= regexp(str,'\s+','split');
