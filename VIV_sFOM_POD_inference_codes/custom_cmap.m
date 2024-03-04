vec = [100;  80;  70;   50;   30;  20;   0];
hex = [ '#0323AA';  '#0572CD';  '#4DBECE'; '#FFFFFF'; '#DA9777';'#C54833'; '#BF1008'];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 512;
map = interp1(vec,raw,linspace(100,0,N),'pchip');

vec = [100; 85;  50;  25;  10;   0];
hex = ['#4DBECE';'#FFFFFF'; '#DA9777'; '#C54833'; '#BF1008'; '#A70708';];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 512; 
map2 = interp1(vec,raw,linspace(100,0,N),'pchip');
