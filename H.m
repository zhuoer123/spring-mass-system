function H = H(x,v,rou)
%计算H飘
%   此处显示详细说明
H=[x/rou,   0;
   v/rou-x^2*v/rou^3,   x/rou];
end

