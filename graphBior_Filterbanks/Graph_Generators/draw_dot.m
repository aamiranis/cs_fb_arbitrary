function [x,y] = draw_dot(cm)
  bg = biograph(cm);
  dolayout(bg)
  xy = cell2mat(get(bg.Nodes,'Position'));
  x = xy(:,1);
  y = xy(:,2);
end % function draw_d