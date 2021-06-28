function [] = plotConstraints(A,b,x1,x2)
    [~, m] = size(A);
    axis([x1(1) x1(end) x2(1) x2(end)])

    hold on
    for i = 1:m
        if (A(2,i)) == 0
            if (A(1,i)) < 0
                fill([0 0 x1(end) x1(end)],[x2(1) x2(end) x2(end) x2(1)],[0.7 0.7 0.7],'facealpha',0.2)
            elseif (A(1,i)) > 0
                fill([x1(1) x1(1) 0 0],[x2(1) x2(end) x2(end) x2(1)],[0.7 0.7 0.7],'facealpha',0.2)
            end
        else
            x2ci = @(x1) (-A(1,i)*x1+b(i))/A(2,i);
            if (A(2,i)) < 0
                fill([x1(1) x1 x1(end)],[x2(end) x2ci(x1) x2(end)],[0.7 0.7 0.7],'facealpha',0.2)
            elseif (A(2,i)) > 0
                fill([x1(1) x1 x1(end)],[x2(1) x2ci(x1) x2(1)],[0.7 0.7 0.7],'facealpha',0.2)
            end

        end
    end
    hold off
end