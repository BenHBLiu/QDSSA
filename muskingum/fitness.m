function f=fitness(x)
%% 5
% I=[22 23 35 71 103 111 109 100 86 71 59 47 39 32 28 24 22 21 20 19 19 18];
% Q=[22 21 21 26 34 44 55 66 75 82 85 84 80 73 64 54 44 36 30 25 22 19];
% K=x(1);
% X=x(2);
% alpha1=x(3);
% alpha2=x(4);
% beta=x(5);
% St=K*(X*22^alpha1+(1-X)*22^alpha2)^beta;
% O(1)=22;
% for i=1:21
%     deltaS=I(i)-((1/(1-X))*(St/K)^(1/beta)-(X/(1-X))*I(i)^alpha1)^(1/alpha2);
%     St=St+deltaS*6;
%     O(i+1)=((1/(1-X))*(St/K)^(1/beta)-(X/(1-X))*I(i)^alpha1)^(1/alpha2);
% end
%% 
Q=[102 140 169 190 209 218 210 194 172 149 136 228 303 366 456 615 830 969 665 519 444 321 208 176 148 125 114 106 97 89 81 76 71 66];
I=[154 150 219 182 182 192 165 150 128 168 260 471 717 1092 1145 600 365 277 277 187 161 143 126 115 102 93 88 82 76 73 70 67 63 59];
K=x(1);
X=x(2);
alpha1=x(3);
alpha2=x(4);
beta=x(5);
O(1)=102;
St=K*(X*154^alpha1+(1-X)*154^alpha2)^beta;
for i=1:33
    deltaS=I(i)-((1/(1-X))*(St/K)^(1/beta)-(X/(1-X))*I(i)^alpha1)^(1/alpha2);
    St=St+deltaS*6;
    O(i+1)=((1/(1-X))*(St/K)^(1/beta)-(X/(1-X))*I(i)^alpha1)^(1/alpha2);
end
%%
% Q=[118.4 197.4 241.1 402.1 518.2 523.9 603.1 829.7 1124.2 1379 1509.3 1379 1050.6 1013.7 1013.7 1013.7 1209.1 1248.8 1002.4 713.6 464.4 325.6 265.6 222.6];
% I=[166.2 263.6 365.3 580.5 594.7 662.6 920.3 1568.8 1775.5 1489.5 1223.3 713.6 645.6 1166.7 1427.2 1282.8 1098.7 764.6 458.7 351.1 288.8 228.8 170.2 143];
% K=x(1);
% X=x(2);
% alpha1=x(3);
% alpha2=x(4);
% beta=x(5);
% O(1)=118.4;
% St=K*(X*166.2^alpha1+(1-X)*166.2^alpha2)^beta;
% for i=1:23
%     deltaS=I(i)-((1/(1-X))*(St/K)^(1/beta)-(X/(1-X))*I(i)^alpha1)^(1/alpha2);
%     St=St+deltaS*12;
%     O(i+1)=((1/(1-X))*(St/K)^(1/beta)-(X/(1-X))*I(i)^alpha1)^(1/alpha2);
% end
%% SSQ
% f=sum(abs(O-Q));
f=(norm(Q-O))^2 ;%SSQ
end