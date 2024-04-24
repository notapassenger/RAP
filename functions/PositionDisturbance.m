function [x,y]= PositionDisturbance(x0,y0,R,numPoints, method)
    if method == 0
        % rand方法, 点均匀分布
        % 输入依次是：圆心横纵坐标，半径和点的数量
        % 随机生成numPoints个半径
        rng(1);
        r=R*sqrt(rand(1,numPoints));
        % 随机生成numPoints个角度
        rng(2);
        seta=2*pi*rand(1,numPoints);
        % 得到点的坐标
        x=x0+r.*cos(seta);
        y=y0+r.*sin(seta);
    elseif method == 1
        % 坐标点等间隔分布（半径和角度等间隔）
        % 生成sqrt(numPoints)个等间隔半径
        grid = R/sqrt(numPoints);
        r = grid:grid:R;
        % 生成sqrt(numPoints)个等间隔角度
        seta=2*pi.*(r/R);
        % 得到点的坐标
        xtempo = zeros(size(seta, 2), size(r, 2));
        ytempo = xtempo;
        for i = 1:size(r, 2)
            xtempo(:, i)=x0+r.*cos(seta(i));
            ytempo(:, i)=y0+r.*sin(seta(i));  
        end
        x = xtempo(:);
        y = ytempo(:);
    elseif method == 2
        x2 = (-R:2*R/(floor(sqrt(numPoints))-1):R)';
        x = repmat(x2, size(x2, 1), 1);
%         for i = length(y2)
%             y() = repmat(y2(i), length(y2), 1);
%         end
        y2 = repmat(x2', size(x2, 1), 1);
        y = y2(:);
    else
        disp('error');
    end
    theta=0:0.001:360;
    % 利用极坐标得到圆的坐标
    Circle1 = x0 + R*cos(theta);
    Circle2 = y0 + R*sin(theta);
    % 画圆
    figure(2);
    plot(Circle1, Circle2,'r')
    hold on
    % 画出点
    plot(x, y, '*');
    hold off;
end