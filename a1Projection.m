function a1Projection
a0 = [63/64 0;0 3/8];
a1 = [1/2 -15/16;1/16 -7/32];
a2 = [1/128 0;-7/256 1/16];
an1 = [1/2 15/16;-1/16 -7/32];
an2 = [1/128 0;7/256 1/16];

% figure, sphere(1000),shading interp, axis equal, hold on;
x = linspace(0,1,21);
P = F(x);
v = Fp(x);
% plot3(P(1,:), P(2,:), P(3,:), 'k.', 'MarkerSize', 10);

for j = 1:8
    Pn = zeros(3,2*length(P)-1);
    vn = zeros(3,2*length(v)-1);
    head1x = a2'*[P(1,end-1);v(1,end-1)]+a0'*[P(1,1);v(1,1)]+an2'*[P(1,2);v(1,2)];
    head1y = a2'*[P(2,end-1);v(2,end-1)]+a0'*[P(2,1);v(2,1)]+an2'*[P(2,2);v(2,2)];
    head1z = a2'*[P(3,end-1);v(3,end-1)]+a0'*[P(3,1);v(3,1)]+an2'*[P(3,2);v(3,2)];
    [Pn(:,1),vn(:,1)] = ProjectionToS2([head1x(1);head1y(1);head1z(1)],[head1x(2);head1y(2);head1z(2)]);
    head2x = a1'*[P(1,1);v(1,1)]+an1'*[P(1,2);v(1,2)];
    head2y = a1'*[P(2,1);v(2,1)]+an1'*[P(2,2);v(2,2)];
    head2z = a1'*[P(3,1);v(3,1)]+an1'*[P(3,2);v(3,2)];
    [Pn(:,2),vn(:,2)] = ProjectionToS2([head2x(1);head2y(1);head2z(1)],[head2x(2);head2y(2);head2z(2)]);
    for k = 1:length(P)-2
        n1x = a2'*[P(1,k);v(1,k)]+a0'*[P(1,k+1);v(1,k+1)]+an2'*[P(1,k+2);v(1,k+2)];
        n1y = a2'*[P(2,k);v(2,k)]+a0'*[P(2,k+1);v(2,k+1)]+an2'*[P(2,k+2);v(2,k+2)];
        n1z = a2'*[P(3,k);v(3,k)]+a0'*[P(3,k+1);v(3,k+1)]+an2'*[P(3,k+2);v(3,k+2)];
        [Pn(:,2*k+1),vn(:,2*k+1)] = ProjectionToS2([n1x(1);n1y(1);n1z(1)],[n1x(2);n1y(2);n1z(2)]);
        n2x = a1'*[P(1,k+1);v(1,k+1)]+an1'*[P(1,k+2);v(1,k+2)];
        n2y = a1'*[P(2,k+1);v(2,k+1)]+an1'*[P(2,k+2);v(2,k+2)];
        n2z = a1'*[P(3,k+1);v(3,k+1)]+an1'*[P(3,k+2);v(3,k+2)];
        [Pn(:,2*k+2),vn(:,2*k+2)] = ProjectionToS2([n2x(1);n2y(1);n2z(1)],[n2x(2);n2y(2);n2z(2)]);
    end
    Pn(:,end) = Pn(:,1);
    vn(:,end) = vn(:,1);
    P = Pn;
    v = vn;
end
% plot3(P(1,:), P(2,:), P(3,:)); title('Extrinsic Hermite Projection');
% figure, plot(P(1,:)); title('First Component');
y = diff(P(1,:),2)/2^(-j*2);
figure, plot(linspace(0,1,length(y)),y); title(['2nd Order Divided Difference for the First Component, Subdivided ' num2str(j) ' times']);

function y = F(x)
y = [sin(cos(2*pi*x)).*cos(sin(2*pi*x) + cos(6*pi*x)); ...
    sin(cos(2*pi*x)).*sin(sin(2*pi*x) + cos(6*pi*x)); ...
    cos(cos(2*pi*x))];

function v = Fp(x)
v = [-2*pi*sin(2*pi*x).*cos(cos(6*pi*x) + sin(2*pi*x)).*cos(cos(2*pi*x)) - 2*pi*sin(cos(6*pi*x) + sin(2*pi*x)).*sin(cos(2*pi*x)).*(cos(2*pi*x) - 3*sin(6*pi*x));...
    2*pi*cos(cos(6*pi*x) + sin(2*pi*x)).*sin(cos(2*pi*x)).*(cos(2*pi*x) - 3*sin(6*pi*x)) - 2*pi*sin(2*pi*x).*cos(cos(2*pi*x)).*sin(cos(6*pi*x) + sin(2*pi*x));...
    2*pi*sin(2*pi*x).*sin(cos(2*pi*x))];
v = v/norm(v); % ???

function [px,pv] = ProjectionToS2(x,v)
px = x/norm(x);
normx = (x(1)^2+x(2)^2+x(3)^2)^(1/2);
normx_32 = (x(1)^2+x(2)^2+x(3)^2)^(3/2);
d = [ 1/normx-x(1)^2/normx_32, -x(1)*x(2)/normx_32, -x(1)*x(3)/normx_32; ...
    -x(2)*x(1)/normx_32, 1/normx-x(2)^2/normx_32, -x(2)*x(3)/normx_32; ...
    -x(3)*x(1)/normx_32, -x(3)*x(2)/normx_32, 1/normx-x(3)^2/normx_32];
pv = d*v;
