N=64;% l=0,1,2,3
q=ones(1,64);
a_x=linspace(1,2,8);
x=[a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),a_x(1),...
    a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),a_x(2),...
    a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),a_x(3),...
    a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),a_x(4),...
    a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),a_x(5),...
    a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),a_x(6),...
    a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),a_x(7),...
    a_x(8),a_x(8),a_x(8),a_x(8),a_x(8),a_x(8),a_x(8),a_x(8)];
y=[a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8),...
    a_x(1),a_x(2),a_x(3),a_x(4),a_x(5),a_x(6),a_x(7),a_x(8)];
% Q_cell   is checked
% N_cell   is checked
% z        is checked
% a_1      is checked
% step 1   is checked