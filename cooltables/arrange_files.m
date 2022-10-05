% script written to arrange the cooling data files with proper way
%-----------------------------------------------------------------------
clear all;
format short e
%---------------------------------------------
% Reading the data files
%---------------------------------------------
data1 = load ('sd_0_tef.dat');
data2 = load ('sd_m_0.5_tef.dat');
data3 = load ('sd_m_1_tef.dat');
data4 = load ('sd_m_1.5_tef.dat');
data5 = load ('sd_m_2_tef.dat');
data6 = load ('sd_m_3_tef.dat');
%-----------------------------------------------
% saving them to variables
%-----------------------------------------------
c1 = [10.^data1(:,1) 10.^data1(:,2)];
c2 = [10.^data2(:,1) 10.^data2(:,2)];
c3 = [10.^data3(:,1) 10.^data3(:,2)];
c4 = [10.^data4(:,1) 10.^data4(:,2)];
c5 = [10.^data5(:,1) 10.^data5(:,2)];
c6 = [10.^data6(:,1) 10.^data6(:,2)];
%-----------------------------------------------
% dumping on files
%-----------------------------------------------
save 'sd_0.dat' c1 -ASCII
save 'sd_m_0.5.dat' c2 -ASCII
save 'sd_m_1.dat' c3 -ASCII
save 'sd_m_1.5.dat' c4 -ASCII
save 'sd_m_2.dat' c5 -ASCII
save 'sd_m_3.dat' c6 -ASCII