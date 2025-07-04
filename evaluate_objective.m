% function f = evaluate_objective(x, M, V)%ZDT1,去main.m里把维数V改成30
% f = [];
% f(1) = x(1);
% g = 1;
% sum = 0;
% for i = 2:V
%     sum = sum + x(i);
% end
% sum = 9*(sum / (V-1));
% g = g + sum;
% f(2) = g * (1 - sqrt(x(1) / g));
% 
% end

% function f = evaluate_objective(x, M, V)%ZDT2，去main.m里把维数V改成30
% f = [];
% f(1) = x(1);
% g = 1;
% sum = 0;
% for i = 2:V
%     sum = sum + x(i);
% end
% sum = 9*(sum / (V-1));
% g = g + sum;
% f(2) = g .* (1 - (x(1)./ g).^2); 
% end

% function f = evaluate_objective(x, M, V)%ZDT3，去main.m里把维数V改成30
% f = [];
% f(1) = x(1);
% g = 1;
% sum = 0;
% for i = 2:V
%     sum = sum + x(i);
% end
% sum = 9*(sum / (V-1));
% g = g + sum;
% f(2) = g * (1 - sqrt(x(1) / g)-(x(1) / g).*sin(10*pi.*x(1)));
% end

% function f = evaluate_objective(x, M, V)  % ZDT4，去main.m里把维数V改成10
%     f = [];
%     % f1 计算
%     f(1) = x(1);
%     
%     % g 函数计算 (注意: x1∈[0,1], x2..xn∈[-5,5])
%     g = 1;
%     sum_val = 0;  % 避免使用 sum 作为变量名
%     for i = 2:V
%         sum_val = sum_val + (x(i)^2 - 10*cos(4*pi*x(i)));
%     end
%     g = g + 10*(V-1) + sum_val;  % 标准ZDT4公式
%     
%     % f2 计算
%     f(2) = g * (1 - sqrt(x(1)/g));
% end

% function f = evaluate_objective(x, M, V)  % ZDT6，去main.m里把维数V改成10
%     f = [];
%     % f1 计算 (非线性)
%     f(1) = 1 - exp(-4*x(1)) * (sin(6*pi*x(1)))^6;
%     
%     % g 函数计算 (带指数0.25)
%     g = 0;
%     for i = 2:V
%         g = g + x(i);
%     end
%     g = 1 + 9 * (g/(V-1))^0.25;  % 关键指数0.25
%     
%     % f2 计算
%     f(2) = g * (1 - (f(1)/g)^2);
% end

% function f = evaluate_objective(x, M, V)%DTLZ1
% f = [];
% sum =0;
% for i = 1:V
%     sum = sum + (x(i)-0.5).^2-cos(20*pi*(x(i))-0.5);
% end
% g = 100*(V+sum);
% f(1) = 0.5*x(1).*x(2).*(1+g);
% f(2) = 0.5*x(1).*(1-x(2)).*(1+g);
% f(3) = 0.5*(1-x(1)).*(1+g);
% end
% function f = evaluate_objective(x, M, V)%DTLZ1
%     % 计算g函数（只使用位置变量 x(M:V)）
%     sum_val = 0;
%     for i = M:V
%         % 修正括号位置：在cos函数内添加括号
%         term = (x(i) - 0.5)^2 - cos(20 * pi * (x(i) - 0.5));  % ✅ 关键修正
%         sum_val = sum_val + term;
%     end
%     g = 100 * ((V - M + 1) + sum_val); % 添加常数项 (V-M+1)
%     
%     % 初始化目标函数向量
%     f = zeros(1, M);
%     
%     % 计算第一个目标
%     f(1) = 0.5 * prod(x(1:M-1)) * (1 + g);
%     
%     % 计算中间目标 (2 到 M-1)
%     for j = 2:(M-1)
%         f(j) = 0.5 * prod(x(1:M-j)) * (1 - x(M-j+1)) * (1 + g);
%     end
%     
%     % 计算最后一个目标
%     f(M) = 0.5 * (1 - x(1)) * (1 + g);
% end

% function f = evaluate_objective(x, M, V)%DTLZ2
%     % 计算g函数（位置变量 x(M:V)）
%     g = sum((x(M:V) - 0.5).^2);
%     
%     % 初始化目标函数向量
%     f = ones(1, M) * (1 + g);  % 预乘(1+g)项
%     
%     % 计算乘积项（角度变量 x(1:M-1)）
%     for i = 1:M
%         % 前M-i个变量取余弦
%         for j = 1:(M - i)
%             f(i) = f(i) * cos(x(j) * pi / 2);
%         end
%         % 第(M-i+1)个变量取正弦（当i>1时）
%         if i > 1
%             f(i) = f(i) * sin(x(M - i + 1) * pi / 2);
%         end
%     end
% end
% function f = evaluate_objective(x, M, V)%DTLZ3
%     % 计算g函数（位置变量 x(M:V)）
%     sum_val = 0;
%     for i = M:V
%         term = (x(i) - 0.5)^2 - cos(20 * pi * (x(i) - 0.5));
%         sum_val = sum_val + term;
%     end
%     g = 100 * ((V - M + 1) + sum_val);  % 包含常数项
%     
%     % 初始化目标函数向量
%     f = ones(1, M) * (1 + g);  % 预乘(1+g)项
%     
%     % 计算乘积项（角度变量 x(1:M-1)）
%     for i = 1:M
%         % 前M-i个变量取余弦
%         for j = 1:(M - i)
%             f(i) = f(i) * cos(x(j) * pi / 2);
%         end
%         % 第(M-i+1)个变量取正弦（当i>1时）
%         if i > 1
%             f(i) = f(i) * sin(x(M - i + 1) * pi / 2);
%         end
%     end
% end
% 
% function f = evaluate_objective(x, M, V)%DTLZ4
%     alpha = 100;  % 密度控制参数
%     g = sum((x(M:V) - 0.5).^2);
%     
%     f = ones(1, M);
%     for i = 1:M
%         for j = 1:(M - i)
%             f(i) = f(i) * cos(x(j)^alpha * pi/2);
%         end
%         if i > 1
%             f(i) = f(i) * sin(x(M - i + 1)^alpha * pi/2);
%         end
%         f(i) = f(i) * (1 + g);
%     end
% end

% function f = evaluate_objective(x, M, V)
%     alpha = 100;  % 保持高密度参数
%     
%     % 增强位置变量的alpha效果
%     transformed_x = zeros(1, M-1);
%     for j = 1:(M-1)
%         % 使用sigmoid函数增强边界效应
%         sig = 1 / (1 + exp(-20*(x(j)-0.5))); % 强化边界吸引
%         transformed_x(j) = sig^alpha;
%     end
%     
%     g = sum((x(M:V) - 0.5).^2);
%     
%     f = ones(1, M);
%     for i = 1:M
%         for j = 1:(M - i)
%             f(i) = f(i) * cos(transformed_x(j) * pi/2);
%         end
%         if i > 1
%             f(i) = f(i) * sin(transformed_x(M - i + 1) * pi/2);
%         end
%         f(i) = f(i) * (1 + g);
%     end
% end

% function f = evaluate_objective(x, M, V)%DTLZ5
%     g = sum((x(M:V) - 0.5).^2);
%     
%     % 计算特殊角度映射
%     theta(1) = x(1) * pi/2;
%     for i = 2:M-1
%         theta(i) = pi/(4*(1 + g)) * (1 + 2*g*x(i));
%     end
%     
%     f = ones(1, M);
%     for i = 1:M
%         for j = 1:(M - i)
%             f(i) = f(i) * cos(theta(j));
%         end
%         if i > 1
%             f(i) = f(i) * sin(theta(M - i + 1));
%         end
%         f(i) = f(i) * (1 + g);
%     end
% end

% function f = evaluate_objective(x, M, V)%DTLZ6
%     g = sum(x(M:V).^0.1);  % 不同的g函数
%     
%     % 计算特殊角度映射 (同DTLZ5)
%     theta(1) = x(1) * pi/2;
%     for i = 2:M-1
%         theta(i) = pi/(4*(1 + g)) * (1 + 2*g*x(i));
%     end
%     
%     f = ones(1, M);
%     for i = 1:M
%         for j = 1:(M - i)
%             f(i) = f(i) * cos(theta(j));
%         end
%         if i > 1
%             f(i) = f(i) * sin(theta(M - i + 1));
%         end
%         f(i) = f(i) * (1 + g);
%     end
% end
% 
% function f = evaluate_objective(x, M, V)%DTLZ7
%     g = 1 + 9/(V - M + 1) * sum(x(M:V));
%     
%     % 前M-1个目标
%     f(1:M-1) = x(1:M-1);
%     
%     % 计算h函数
%     h_sum = 0;
%     for i = 1:M-1
%         h_sum = h_sum + f(i)/(1 + g) * (1 + sin(3 * pi * f(i)));
%     end
%     h = M - h_sum;
%     
%     % 最后一个目标
%     f(M) = (1 + g) * h;
% end

% function f = evaluate_objective(x, M, V)%UF1
% f = zeros(1, M);
% count1 = 0; count2 = 0;
% sum1 = 0; sum2 = 0;
% 
% for j = 3:V
%     yj = x(j) - sin(6*pi*x(1) + j*pi/V);
%     if mod(j, 2) == 1
%         sum1 = sum1 + yj^2;
%         count1 = count1 + 1;
%     else
%         sum2 = sum2 + yj^2;
%         count2 = count2 + 1;
%     end
% end
% 
% f(1) = x(1) + 2 * sum1 / max(count1, 1e-5);
% f(2) = 1 - sqrt(x(1)) + 2 * sum2 / max(count2, 1e-5);
% end

function f = evaluate_objective(x, M, V)
f = zeros(1, M);
prod1 = 1; prod2 = 1;
sum1 = 0; sum2 = 0;

for j = 2:V
    yj = x(j) - x(1)^(0.5*(1 + 3*(j-2)/(V-2)));
    pj = cos(20*yj*pi/sqrt(j));
    if mod(j, 2) == 1
        prod1 = prod1 * pj;
        sum1 = sum1 + yj^2;
    else
        prod2 = prod2 * pj;
        sum2 = sum2 + yj^2;
    end
end

f(1) = x(1) + 2*(4*sum1 - 2*prod1 + 2) / max(1, V/2);
f(2) = 1 - sqrt(x(1)) + 2*(4*sum2 - 2*prod2 + 2) / max(1, V/2);
end