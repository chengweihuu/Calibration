clear all
close all
%%  此函数为求解机器人手眼矩阵的主函数

%% 读取测量数据
% 1、机器人信息
% 2、相机坐标系下圆心位置
% 3、验证的机器人信息
% 4、验证的相机坐标系下圆心位置
% 5、标定球的直径
D0 = 30.006;                    % 标定球的直径（mm）
R0 = D0/2;                      % 表示标定球半径；
all_data2 = xlsread('robort_verify1.xlsx');% 标定所需的机器人六轴数据
s = size(all_data2,1);
all_data1 = cell(1,s);          % 存放线激光数据
all_data0 = all_data2;

all_data3 = xlsread('circle_verify1.xlsx');  % 存放拟合圆心位置坐标 ，拟合半径大小，正负标志等；

%% 点云数据拟合并求取圆心位置
s0 = length(all_data3);
% 提取 r 
R_data = cell(1,s0);
for i = 1:length(all_data3)
    R_data{1,i} = all_data3(i,3);
end
% 计算球心 y 方向的值
y_data = cell(1,s);
for k = 1:length(R_data)
    
    ss = all_data3(k,4);
    Rs = R_data{1,k};
    y_c = sqrt(R0^2-Rs^2);  % 这里需要注意y_c的正负值；
    if ss == 1
        y_data{1,k} = y_c;
    else
        y_data{1,k} = -y_c;
    end
    
end
y_mat = cell2mat(y_data)';
% 合并计算球心坐标
center_S = [];
for j = 1:length(all_data3)
   
    coor = [all_data3(j,1),y_mat(j),all_data3(j,2)];
    center_S = [center_S;coor];
    
end
ns = length(center_S);
ee = ones(ns,1);
data2 = [center_S,ee];

for i =1:size(all_data2,1)
    %% 机器人在各个位姿下的机器人齐次变换矩阵
    T_ur = zhengyundongxue(deg2rad(all_data2(i,1:6)),6);
    M{i}=T_ur;
    
end

syms a1 a2 a3 a4 a5 a6 a7 a8 a9 d1 d2 d3;

X=[a1 a2 a3 d1;a4 a5 a6 d2;a7 a8 a9 d3; 0 0 0 1];
m = size(data2,1)-1;
P = zeros(4,12);
q = zeros(4,1);
s = 1;

for i= 1:m 
    % 基于最小二乘，通过列方程求解设置参数。
    Pa=data2(i,:)'; Pb=data2(i+1,:)';
    A=M{1,i};B=M{1,i+1};
    result=A*X*Pa-B*X*Pb;
    a = vpa(result,4);
    [Z,b] = equationsToMatrix(a, [a1 a2 a3 a4 a5 a6 a7 a8 a9 d1 d2 d3]);
    for t = 1:4
        for j = 1:12
            P(t,j) = Z(t,j);
        end
    end
    G{i} = P;
    for s = 1:4
        q(s,1)= b(s,1);
    end
    g{i} = q;
    
end
 
%% 求解未知数

g = cell2mat(g);
B = reshape(g,[4*m,1]);
A = zeros(4*m,12);
G = cell2mat(G);
for h = 1:m
   A(4*h-3:4*h,:) =G(:,12*h-11:12*h);
end
C=A\B;        % 这是反斜杠
C=[C(1:3,1)' C(10);C(4:6,1)' C(11); C(7:9,1)' C(12); 0 0 0 1];
% save C.mat

%% 结果显示
disp("机器人手眼矩阵为：");
RR = vpa(C,4)
xlswrite('相机坐标位置.xlsx',center_S);

% 相机坐标系 转换到 机器人坐标系下的球心
data_base = [];
for jj = 1:length(center_S)
    
    temp_s = center_S(jj,:);
    temp_s = [temp_s,1];
    R0 = M{jj};
    temp_b= R0*C*temp_s';
    temp_b = temp_b';
    temp_b(end) = [];
    data_base = [data_base;temp_b];
    
end
xlswrite('机器人坐标位置.xlsx',data_base);

T_se1 = C ;
% save T_se1.mat


