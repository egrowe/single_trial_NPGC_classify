function plot_JackKnife(arg,col,solid,t_axis)
% adapted from JackKnife written by Chandrasekaran used by Erie
% adapted by Mona
% adapted by Marta Garrido
% arg: matrix to plot (subjexts x timepoints)
% col: color
% solid: linestyle


if solid
    linest = '-';
else
    linest = '--';
end
y = mean(arg,1);

% x=1:length(arg);
x=t_axis;

% se = std((arg))./sqrt(size(arg,1)-1);
se = std((arg))./sqrt(size(arg,1));

L = y - se;
U = y + se;
Xcoords = [x x(end:-1:1)];
Ycoords = [U L(end:-1:1)];

xlabel('Time (ms)');
ylabel('Amplitude (uV)');

%shadow color - add white noise (0.6 * diff to [1 1 1])
col2 = col;
col2(find(col<1)) = col2(find(col<1))+ 0.6*[1-col(find(col<1))];

hold on;
Pa = patch(Xcoords,Ycoords,col2);
set(Pa,'linestyle','none','linewidth',3);
Li = plot(x,y,'color',col,'linewidth',3,'LineStyle',linest);
alpha(0.5)


