function thesis_init

thesis_data;

b_init_matrix = zeros(3,3,3);
b_init_matrix(:,:,1) = [2 1 0;1 2 0;0 0 0];
b_init_matrix(:,:,2) = [0 0 0;0 2 1;0 1 2];
b_init_matrix(:,:,3) = [2 0 1;0 0 0;1 0 2];

d_init_matrix = zeros(3,3,3);
d_init_matrix(:,:,1) = [1 1 0;1 1 0;0 0 0];
d_init_matrix(:,:,2) = [0 0 0;0 1 1;0 1 1];
d_init_matrix(:,:,3) = [1 0 1;0 0 0;1 0 1];

Prec = 0.001;
Gauss = [-0.5773502692, 0.5773502692];

SpheroidDomain = 0;
RegularMesh = 0;
CavityChoice = 3;
MethodChoice = 1;

nodes = 1;
nodes_x = 4;
nodes_y = 9;%2n+1!
Fea_Len = 1;
alpha = 0.75;
w0 = 3;
G0 = 0;
R = 1;
H = 0.8;
M0 = 1;

B_mesh = 0;
B_spectrum = 0;