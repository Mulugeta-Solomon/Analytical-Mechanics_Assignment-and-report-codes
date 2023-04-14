% dynamic deformation of a beam (constant cross-sectional area�j
% time [ 0, tpush ] : bottom end fixed. force fpush pushes top end
% time [ tpush, tend ] : release bottom end. no force to top.
%                        bottom contacting to floor causes reaction force
%                        floor reaction force is described by penatly method
% g, cm, msec

% �r�[���̓��I�ό`�i�f�ʐψ��j
% ���� [ 0, tpush ] : ���[�Œ�@��[�͗� fpush �ŉ���
% ���� [ tpush, tend ] : ���[�̐�����O���D��[�̗͂��Ȃ����D
%                        ���[�����ɐڐG���Ă���Ƃ��́C�����甽�͂��󂯂�
%                        �����͂̓y�i���e�B�@�ŕ\��
% g, cm, sec

L = 10; A = 2;  % beam length and cross-sectional area (constant)
                % �r�[���̒����ƒf�ʐρi���j
n = 6;          % number of nodal points	% �ߓ_�̌�
h = L/(n-1);

% E = 50 kPa; c = 0.2 kPa s; rho = 1 g/cm^2
E = 0.5*1e+6;   % elastic modulus (Young's modulus)	% �e�����i�����O���j
c = 2.0*1e+3;   % viscous modulus			% �S����
rho = 1;        % density                   % ���x
alpha = 2000;   % positive constant for CSM	% ������艻�@�̒萔

% fpush = 2 N; kfloor = 0.1 MPa/cm
fpush = 2.0e+5;     % pushing force			% �������ݗ�
kfloor = 1.0*1e+6;  % floor elasticity		% ���̒e��
tpush = 0.2; tend = 0.4;

e0 = (4/6)*ones(n,1); e0(1) = (2/6); e0(n) = (2/6);
e1 = (1/6)*ones(n,1);
M = (rho*A*h)*spdiags([e1 e0 e1], -1:1, n, n);  % inertia matrix   % �����s��

e0 = 2*ones(n,1); e0(1) = 1; e0(n) = 1;
e1 = (-1)*ones(n,1);
K = (E*A/h)*spdiags([e1 e0 e1], -1:1, n, n);    % stiffness matrix % �����s��
B = (c*A/h)*spdiags([e1 e0 e1], -1:1, n, n);    % damping matrix   % �S���s��

% time interval [0,tpush]
% ���ԋ�� [0,tpush]
constraints = zeros(n,1);
constraints(1) = 1;  % constraint matrix (bottom fixed) % ����s��i���[�Œ�j
p = n;               % force application point (top)    % �͂̍�p�_�i��[�j
external_force = @(t,n,un,vn) pushing_force_param (t, n,un,vn, p,fpush,tpush); % external force % �O��
beam_dynamic_equation = @(t,q) beam_supported_dynamic_equation_param (t,q, n, M,B,K, constraints, alpha, external_force);
interval = [0,tpush];
qinit = [ zeros(n,1); zeros(n,1) ];
[timepush,qpush] = ode45(beam_dynamic_equation, interval, qinit);

% time interval [tpush,tend]
% ���ԋ�� [tpush,tend]
p = 1;               % force application point (bottom)  % �͂̍�p�_�i���[�j
external_force = @(t,n,un,vn) floor_force_param(t,n,un,vn, p,A,kfloor); % floor reaction force to bottom % ���[�̏�����
beam_dynamic_equation = @(t,q) beam_free_dynamic_equation_param (t,q, n,M,B,K,external_force);
interval = [tpush,tend];
qinit = qpush(end, :);
[timefree,qfree] = ode45(beam_dynamic_equation, interval, qinit);

% whole time interval
% �S���ԋ��
time = [ timepush; timefree ];
q = [ qpush; qfree ];

figure('position', [0, 0, 400, 400]);
set(0,'defaultAxesFontSize',16);
set(0,'defaultTextFontSize',16);

clf;
hold on;
for k=1:n
    plot(time,q(:,k)+(k-1)*h);
end
hold off;
xlabel("time");
ylabel("position");
xlim([0,tend]);
ylim([-2,28]);
grid on;
saveas(gcf,'beam_pushed_released_position.png');
fprintf("nodal point position / �ߓ_�̈ʒu\n");

clf;
hold on;
for t = 0:0.02:tend
    i = nearest_index(time,t);
    plot([time(i);time(i)], [q(i,1)+(1-1)*h;q(i,n)+(n-1)*h], 'b-o');
end
%plot(time, q(:,1), time, q(:,n)+(n-1)*h);
hold off;
xlabel("time");
ylabel("position");
xlim([0,tend]);
ylim([-2,28]);
grid on;
saveas(gcf,'beam_pushed_released_body.png');
fprintf("deformation and motion of body / �{�f�B�̕ό`�Ɖ^��\n");

%plot(time, q(:,n)-q(:,1)+L, time, L*ones(size(time,1),1),'--');
%xlabel("time");
%ylabel("length");
%ylim([6,12]);
%grid on;
%saveas(gcf,'beam_pushed_released_body_length.png');
%fprintf("body length / �{�f�B�̒���\n");

