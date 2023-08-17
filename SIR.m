% Compute R0 = beta/gamma of US,FR,JP,TW
% from SIR model we get:
%
% dI/dt = beta*I*(N-I-R)/N - gamma*I
% dR/dt = gamma*I
% 
% Let dI/dt = I_1 - I_0 = beta*I_0*(N-I_0-R_0)/N-gamma*I_0 
%     dR/dt = R_1 - R_0 = gamma*I
% We can use fsolve to solve
% F1(beta, gamma) = beta*I_0*(N-I_0-R_0)/N-gamma*I_0 - dI/dt = 0
% F2(beta, gamma) = gamma*I - dR/dt = 0
%

% US
N = 329227746; % population (constant)
% import data
DATA = csvread('US.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
days = length(I);
% start to compute every day's R_0
count = 0;
sum = 0;
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = Inf;
    else
        R_0(k) = X(1)/X(2);
        if R_0(k) >= 10
            R_0(k) = 10;
        elseif R_0(k) <= 0
            R_0(k) = 0;
        end
        count = count+1;
        sum = sum + R_0(k);
    end    
end
US_a = sum/count;
d = 1:days-1;
subplot(4,1,1)
plot(d, R_0,'b')
hold on
for k = 2:days-2
    if isinf(R_0(k-1))&&isinf(R_0(k+1))
        R_1(k) = R_0(k);
    else
        R_1(k) = NaN;
    end    
end
R_1(days-1) = R_0(days-1);
plot(d,R_1,'bx')
hold on
for k = 1:days-1
    if isinf(R_0(k))
        R_2(k) = 0;
    else
        R_2(k) = NaN;
    end    
end
plot(d,R_2,'*')
grid on
legend('US')
xlabel('date')          
ylabel('R_0')
axis([0, 325, 0, 10])

% FR
N = 68809894; % population (constant)
% import data
DATA = csvread('FR.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
% start to compute every day's R_0
count = 0;
sum = 0;
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = Inf;
    else
        R_0(k) = X(1)/X(2);
         if R_0(k) >= 10
            R_0(k) = 10;
        elseif R_0(k) <= 0
            R_0(k) = 0;
         end
        count = count+1;
        sum = sum + R_0(k);
    end    
end
FR_a = sum/count;
subplot(4,1,2)
plot(d, R_0,'g')
hold on
for k = 2:days-2
    if isinf(R_0(k-1))&&isinf(R_0(k+1))
        R_1(k) = R_0(k);
    else
        R_1(k) = NaN;
    end    
end
plot(d,R_1,'gx')
hold on
for k = 1:days-1
    if isinf(R_0(k))
        R_2(k) = 0;
    else
        R_2(k) = NaN;
    end    
end
plot(d,R_2,'*')
grid on
legend('FR')
xlabel('date')          
ylabel('R_0')   
axis([0, 325, 0, 10])

% JP
N = 124271318; % population (constant)
% import data
DATA = csvread('JP.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
days = length(I);
% start to compute every day's R_0
count = 0;
sum = 0;
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = Inf;
    else
        R_0(k) = X(1)/X(2);
         if R_0(k) >= 10
            R_0(k) = 10;
        elseif R_0(k) <= 0
            R_0(k) = 0;
         end
        count = count+1;
        sum = sum + R_0(k);
    end    
end
JP_a = sum/count;
subplot(4,1,3)
plot(d, R_0,'r')
hold on
for k = 2:days-2
    if isinf(R_0(k-1))&&isinf(R_0(k+1))
        R_1(k) = R_0(k);
    else
        R_1(k) = NaN;
    end    
end
plot(d,R_1,'rx')
hold on
for k = 1:days-1
    if isinf(R_0(k))
        R_2(k) = 0;
    else
        R_2(k) = NaN;
    end    
end
plot(d,R_2,'*')
grid on
legend('JP')
xlabel('date')          
ylabel('R_0')   
axis([0, 325, 0, 10])

% TW
N = 23563356; % population (constant)
% import data
DATA = csvread('TW.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
% start to compute every day's R_0
count = 0;
sum = 0;
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = Inf;
    else
        R_0(k) = X(1)/X(2);
         if R_0(k) >= 10
            R_0(k) = 10;
        elseif R_0(k) <= 0
            R_0(k) = 0;
         end
        count = count+1;
        sum = sum + R_0(k);
    end    
end
TW_a = sum/count;
subplot(4,1,4)
plot(d, R_0,'m')
hold on
for k = 2:days-2
    if isinf(R_0(k-1))&&isinf(R_0(k+1))
        R_1(k) = R_0(k);
    else
        R_1(k) = NaN;
    end    
end
plot(d,R_1,'mx')
hold on
for k = 1:days-1
    if isinf(R_0(k))
        R_2(k) = 0;
    else
        R_2(k) = NaN;
    end    
end
plot(d,R_2,'*')
grid on
legend('TW')
xlabel('date')          
ylabel('R_0')   
axis([0, 325, 0, 10])
disp(['US average R0 = ', num2str(US_a)]);
disp(['FR average R0 = ', num2str(FR_a)]);
disp(['JP average R0 = ', num2str(JP_a)]);
disp(['TW average R0 = ', num2str(TW_a)]);

% solve:
% F1(beta, gamma) = beta*I_0*(N-I_0-R_0)/N-gamma*I_0 - dI/dt = 0
% F2(beta, gamma) = gamma*I - dR/dt = 0
function F = IR_equations(X, N, I, R, dIdt,dRdt)
    % beta = X(1)
    % gamma = X(2)
    F(1,1) = X(1)*I*(N-I-R)/N-X(2)*I-dIdt;
    F(2,1) = X(2)*I-dRdt;
end