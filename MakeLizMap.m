x = -1.8:.01:1.8;

r = max(min(1.5-abs(x-0.5),1),0);
b = max(min(1.5-abs(x+0.5),1),0);
g = min(r,b);

lizmap = [r' g' b'];
