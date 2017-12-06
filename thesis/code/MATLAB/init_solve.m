function init_solve

thesis_data;

i = size(t_p,2);
j = sum(free_point,2);
freq_c = zeros(4*i,1);
freq_d = zeros(4*(i+j),1);

for M0=[-1 1],
   
A = zeros(size(t_p,2));
B = zeros(size(t_p,2));
C = zeros(size(t_p,2));
D = zeros(size(t_p,2));

for i=1:size(t_t,2),
   gen_ac(i);
   gen_b(i);
   gen_d(i);

   for j=1:3,
      for k=1:3,
         A(t_t(j,i),t_t(k,i))=A(t_t(j,i),t_t(k,i))+Ae(j,k);
         B(t_t(j,i),t_t(k,i))=B(t_t(j,i),t_t(k,i))+Be(j,k);
         C(t_t(j,i),t_t(k,i))=C(t_t(j,i),t_t(k,i))+Ce(j,k);
         D(t_t(j,i),t_t(k,i))=D(t_t(j,i),t_t(k,i))+De(j,k);
      end
   end
end

i = size(t_p,2);
j = sum(free_point,2);

G = zeros(2*(i+j));

A11 = A(1:j,1:j);
A12 = A(1:j,j+1:i);
A21 = A(j+1:i,1:j);
A22 = A(j+1:i,j+1:i);

B11 = B(1:j,1:j);
B12 = B(1:j,j+1:i);
B21 = B(j+1:i,1:j);
B22 = B(j+1:i,j+1:i);

C11 = C(1:j,1:j);
C12 = C(1:j,j+1:i);
C21 = C(j+1:i,1:j);
C22 = C(j+1:i,j+1:i);

D11 = D(1:j,1:j);
D(1:j,j+1:i) = 0;
D(j+1:i,1:j) = 0;
D(j+1:i,j+1:i) = 0;

a41 =     A22\C21;
a42 =     A22\C22;
a43 = -M0*A22\B21;
a44 = -M0*A22\B22;
a45 =    -A22\A21;

a61 =   -D11\(C11-A12*(A22\C21))/4;
a62 =   -D11\(C12-A12*(A22\C22))/4;
a63 = M0*D11\(B11-A12*(A22\B21))/4;
a64 = M0*D11\(B12-A12*(A22\B22))/4;
a65 =    D11\(A11-A12*(A22\A21))/4 + eye(j);

G(1:i+j,i+1:i+i+j) = eye(i+j);
G(i+i+1:i+i+j,i+i+j+1:i+i+j+j) = eye(j);

G(i+j+1:i+i,1:j) = a41;
G(i+j+1:i+i,j+1:i) = a42;
G(i+j+1:i+i,i+1:i+j) = a43;
G(i+j+1:i+i,i+j+1:i+i) = a44;
G(i+j+1:i+i,i+i+1:i+i+j) = a45;

G(i+i+j+1:i+i+j+j,1:j) = a61;
G(i+i+j+1:i+i+j+j,j+1:i) = a62;
G(i+i+j+1:i+i+j+j,i+1:i+j) = a63;
G(i+i+j+1:i+i+j+j,i+j+1:i+i) = a64;
G(i+i+j+1:i+i+j+j,i+i+1:i+i+j) = a65;

[e_vector, e_value] = eig(G);

[AA, BB, QQ, ZZ, s_mode] = qz(D11,A11-A12*(A22\A21));
s_mode(j+1:i,:) = -(A22\A21)*s_mode(1:j,:);
s_freq = diag(BB)./diag(AA);
s_freq = sqrt(s_freq)/2;

[AA, BB, QQ, ZZ, u_mode] = qz([A22 zeros(i-j);zeros(i-j) C22],[-M0*B22 C22;C22 zeros(i-j)]);
u_freq = diag(BB)./diag(AA);
u_mode = [zeros(j,2*(i-j)); u_mode(1:i-j,:)];

freq_c(1+(M0+1)*i:(M0+3)*i) = sort(real([s_freq;-1*s_freq;u_freq]));
freq_d(1+(M0+1)*(i+j):(M0+3)*(i+j)) = sort(real(diag(e_value)));

end

freq_c = sort(freq_c);
freq_d = sort(freq_d);

%disp('are you still alive ?'); pause;
init_plot;

return