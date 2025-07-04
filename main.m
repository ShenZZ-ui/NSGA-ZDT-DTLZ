
function nsga_2_optimization
clc
clear all


pop = 100; %��Ⱥ��С
gen = 3000; %��������
M = 2;     %Ŀ�꺯��������
V = 10;    %ά��
min_range = zeros(1, V); 
max_range = ones(1,V); 

%ZDT4ʱ�����ȡֵ��Χ
% min_range = [0, -5*ones(1, V-1)];  % x1��[0,1], ������[-5,5]
% max_range = [1, 5*ones(1, V-1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chromosome = initialize_variables(pop, M, V, min_range, max_range);
chromosome = non_domination_sort_mod(chromosome, M, V);


for i = 1 : gen
    pool = round(pop/2);
    tour = 2;
    parent_chromosome = tournament_selection(chromosome, pool, tour);
    
    mu = 20;
    mum = 20;
    offspring_chromosome = genetic_operator(parent_chromosome,M, V, mu, mum, min_range, max_range);
    [main_pop,~] = size(chromosome);
    [offspring_pop,~] = size(offspring_chromosome);
    
    clear temp
    intermediate_chromosome(1:main_pop,:) = chromosome;
    intermediate_chromosome(main_pop + 1 : main_pop + offspring_pop,1 : M+V) = offspring_chromosome;
    intermediate_chromosome = non_domination_sort_mod(intermediate_chromosome, M, V);
    chromosome = replace_chromosome(intermediate_chromosome, M, V, pop);
    if ~mod(i,100)
        clc;
        fprintf('%d generations completed\n',i);
    end
end
sv1 =chromosome(:,V+1);
sv2 =chromosome(:,V+2);
save('popvalue.mat','chromosome')

if M == 2
    plot(chromosome(:,V + 1),chromosome(:,V + 2),'*');
    xlabel('f_1'); ylabel('f_2');
    title('UF1 Pareto Optimal Front');
%      xlim([0.2,1]);ylim([0,1]);
elseif M == 3
        % ��ȡĿ��ֵ
    f1 = chromosome(:, V+1);
    f2 = chromosome(:, V+2);
    f3 = chromosome(:, V+3);
    %����ֵ
%    DTLZ1_lilunzhi;
%    DTLZ2_lilunzhi;
%    DTLZ3_lilunzhi;
   % ����
    figure('Color', 'white', 'Position', [100, 100, 1200, 500]);
    
    % ʵ�ʽ⼯
    plot3(f1, f2, f3, '.', ...       % '.' ��ʾСԲ��
      'MarkerSize', 12, ...       % ���Ƶ�Ĵ�С���Ƽ�6-12��
      'Color', [0.2, 0.6, 0.9]); % ���õ����ɫ

    hold on;
%     surf(x_plane, y_plane, z_plane, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    title('DTLZ7 Pareto Front');
    xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
%     legend('NSGA||�����', '����������ǰ��ƽ��', 'Location', 'best');
    grid on;
    view(135, 30);
%     axis equal;

%     xlim([0, 1]); ylim([0, 1]); zlim([0, 3.5]);
end