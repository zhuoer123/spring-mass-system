function H = H(x,v,rou)
%����HƮ
%   �˴���ʾ��ϸ˵��
H=[x/rou,   0;
   v/rou-x^2*v/rou^3,   x/rou];
end

