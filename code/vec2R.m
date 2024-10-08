function R = vec2R(Rx,Ry,Rz)
% 输入旋转矢量 (Rx, Ry, Rz)
R_vec = [Rx; Ry; Rz];

% 计算旋转角度 theta
theta = norm(R_vec);

% 检查旋转向量是否为零向量
if theta == 0
    R = eye(3);  % 若旋转向量为零，旋转矩阵为单位矩阵
else
    % 计算旋转轴 k
    k = R_vec / theta;
    
    % 计算反对称矩阵 K
    K = [  0   -k(3)  k(2);
          k(3)  0    -k(1);
         -k(2)  k(1)   0   ];
    
    % 使用罗德里格旋转公式计算旋转矩阵 R
    R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end
