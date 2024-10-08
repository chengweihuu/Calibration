clear all
close all
%%  �˺���Ϊ�����������۾����������

%% ��ȡ��������
% 1����������Ϣ
% 2���������ϵ��Բ��λ��
% 3����֤�Ļ�������Ϣ
% 4����֤���������ϵ��Բ��λ��
% 5���궨���ֱ��
D0 = 30.006;                    % �궨���ֱ����mm��
R0 = D0/2;                      % ��ʾ�궨��뾶��
all_data2 = xlsread('robort_verify1.xlsx');% �궨����Ļ�������������
s = size(all_data2,1);
all_data1 = cell(1,s);          % ����߼�������
all_data0 = all_data2;

all_data3 = xlsread('circle_verify1.xlsx');  % ������Բ��λ������ ����ϰ뾶��С��������־�ȣ�

%% ����������ϲ���ȡԲ��λ��
s0 = length(all_data3);
% ��ȡ r 
R_data = cell(1,s0);
for i = 1:length(all_data3)
    R_data{1,i} = all_data3(i,3);
end
% �������� y �����ֵ
y_data = cell(1,s);
for k = 1:length(R_data)
    
    ss = all_data3(k,4);
    Rs = R_data{1,k};
    y_c = sqrt(R0^2-Rs^2);  % ������Ҫע��y_c������ֵ��
    if ss == 1
        y_data{1,k} = y_c;
    else
        y_data{1,k} = -y_c;
    end
    
end
y_mat = cell2mat(y_data)';
% �ϲ�������������
center_S = [];
for j = 1:length(all_data3)
   
    coor = [all_data3(j,1),y_mat(j),all_data3(j,2)];
    center_S = [center_S;coor];
    
end
ns = length(center_S);
ee = ones(ns,1);
data2 = [center_S,ee];

for i =1:size(all_data2,1)
    %% �������ڸ���λ���µĻ�������α任����
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
    % ������С���ˣ�ͨ���з���������ò�����
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
 
%% ���δ֪��

g = cell2mat(g);
B = reshape(g,[4*m,1]);
A = zeros(4*m,12);
G = cell2mat(G);
for h = 1:m
   A(4*h-3:4*h,:) =G(:,12*h-11:12*h);
end
C=A\B;        % ���Ƿ�б��
C=[C(1:3,1)' C(10);C(4:6,1)' C(11); C(7:9,1)' C(12); 0 0 0 1];
% save C.mat

%% �����ʾ
disp("���������۾���Ϊ��");
RR = vpa(C,4)
xlswrite('�������λ��.xlsx',center_S);

% �������ϵ ת���� ����������ϵ�µ�����
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
xlswrite('����������λ��.xlsx',data_base);

T_se1 = C ;
% save T_se1.mat


