function [] = RemovePlotWhiteArea(gca)
% TightInset��λ��
inset_vectior = get(gca, 'TightInset');
inset_x = inset_vectior(1);
inset_y = inset_vectior(2);
inset_w = inset_vectior(3);
inset_h = inset_vectior(4);

% OuterPosition��λ��
outer_vector = get(gca, 'OuterPosition');
pos_new_x = outer_vector(1) + inset_x; % ��Position��ԭ���Ƶ���TightInset��ԭ��
pos_new_y = outer_vector(2) + inset_y;
pos_new_w = outer_vector(3) - inset_w - inset_x; % ����Position�Ŀ�
pos_new_h = outer_vector(4) - inset_h - inset_y; % ����Position�ĸ�

% ����Position
set(gca, 'Position', [pos_new_x, pos_new_y, pos_new_w, pos_new_h]);







