clear
clc
%% 位置标定
% 指定文件路径
filename = 'data/工具坐标系标定数据.xlsx';
% 读取 Excel 文件中的数据为表格
data = readtable(filename, 'Sheet', 'Sheet2');  % 读取 Sheet1
% 显示读取的数据
disp(data);
PoseArray = table2array(data(1:7,7:12));

Rcell = cell(size(PoseArray,1),1)
pcell = cell(size(PoseArray,1),1)
for i=1:length(Rcell)
    Rcell{i,1} = vec2R(PoseArray(i,4),PoseArray(i,5),PoseArray(i,6));
    pcell{i,1} = [PoseArray(i,1);PoseArray(i,2);PoseArray(i,3)];
end

R_BEcell = cell(3,1);
p_pcell = cell(3,1);
for i=1:3
   R_BEcell{i,1} = Rcell{i,1}-Rcell{i+1,1};
   p_pcell{i,1} = pcell{i+1,1}-pcell{i,1};
end
R_matrix = cell2mat(R_BEcell);
p_matrix = cell2mat(p_pcell);

p_tcp_E = inv(R_matrix' * R_matrix) * R_matrix' * p_matrix

sigma = norm(R_matrix * p_tcp_E - p_matrix)



%% 姿态标定
R_o_BE = vec2R(PoseArray(5,4),PoseArray(5,5),PoseArray(5,6));
R_z_BE = vec2R(PoseArray(6,4),PoseArray(6,5),PoseArray(6,6));
R_x_BE = vec2R(PoseArray(7,4),PoseArray(7,5),PoseArray(7,6));
p_oEo_B = PoseArray(5,1:3)';
p_zEo_B = PoseArray(6,1:3)';
p_xEo_B = PoseArray(7,1:3)';
p_otcp_B = R_o_BE * p_tcp_E + p_oEo_B ;
p_xtcp_B = R_x_BE * p_tcp_E + p_xEo_B ;
p_ztcp_B = R_z_BE * p_tcp_E + p_zEo_B ;

v_x_B = p_xtcp_B - p_otcp_B;
v_z_B = p_ztcp_B - p_otcp_B;

n_ET = inv(R_o_BE)*v_x_B/norm(v_x_B);
a_ET = inv(R_o_BE)*v_z_B/norm(v_z_B);
o_ET = cross(a_ET,n_ET);
o_ET = o_ET / norm(o_ET);
a_ET = cross(n_ET,o_ET);
a_ET = a_ET / norm(a_ET);
R = [n_ET o_ET a_ET]
T_EF = [R p_tcp_E;0 0 0 1]
