%% a)
lambda = 1500; %pps
C = 10; %link capacity Mbps
f = 1000000; %queue size (bytes)
P = 10000; %stopping criterium
n = [40 30 20 10]; %numero de pacotes VoIP
alfa = 0.1; %confidence intervals
N = 50; %numero de vezes que o simulator3 e chamado
size = length(n); %numero de valores de n

PLmatriz = zeros(size,N); %data
PLVoIPmatriz = zeros(size,N); %voip
APDmatriz = zeros(size,N);
APDVoIPmatriz = zeros(size,N);
MPDmatriz = zeros(size,N);
MPDVoIPmatriz = zeros(size,N);
TTmatriz = zeros(size,N);

APDmean = zeros(1,size);
APDVoIPmean = zeros(1,size);
APDterm = zeros(1,size);
APDVoIPterm = zeros(1,size);

for j = 1:size %para cada valor de n
   for i = 1:N %chamamos o simulator3 50vezes para cada valor de n
       [PLmatriz(j,i), PLVoIPmatriz(j,i), APDmatriz(j,i), APDVoIPmatriz(j,i), MPDmatriz(j,i), MPDVoIPmatriz(j,i), TTmatriz(j,i)] = Simulator3(lambda,C,f,P,n(j)); %n(0)=100000, n(1)=20000, ...
   end
   
   APDmean(j) = mean(APDmatriz(j,:)); %media dos valores de cada valor de n (n(0) fica no APDmean(0), n(1) fica no APDmean(1)...
   APDVoIPmean(j) = mean(APDVoIPmatriz(j,:));
   APDterm(j) = norminv(1-alfa/2)*sqrt(var(APDmatriz(j,:))/N);
   APDVoIPterm(j) = norminv(1-alfa/2)*sqrt(var(APDVoIPmatriz(j,:))/N);
   
end

%grafico do packet delay data
figure(1)
bar(n, APDmean)
ylabel('Average packet delay data (ms)')
xlabel('number of VoIP packets')
hold on
er = errorbar(n, APDmean, APDterm);
er.Color = [0 1 0];
er.LineStyle = 'none';
title('Avg Packet Delay data (ms)');
hold off

%grafico do packet delay voip
figure(2)
bar(n, APDVoIPmean)
ylabel('Average packet delay VoIP (ms)')
xlabel('number of VoIP packets')
hold on
er = errorbar(n, APDVoIPmean, APDVoIPterm);
er.Color = [0 1 0];
er.LineStyle = 'none';
title('Avg Packet Delay VoIP (ms)');
hold off

%% b

lambda = 1500; %pps
C = 10; %link capacity Mbps
f = 1000000; %queue size (bytes)
P = 10000; %stopping criterium
n = [40 30 20 10]; %numero de pacotes VoIP
alfa = 0.1; %confidence intervals
N = 50; %numero de vezes que o simulator3 e chamado
size = length(n); %numero de valores de n

PLmatriz = zeros(size,N); %data
PLVoIPmatriz = zeros(size,N); %voip
APDmatriz = zeros(size,N);
APDVoIPmatriz = zeros(size,N);
MPDmatriz = zeros(size,N);
MPDVoIPmatriz = zeros(size,N);
TTmatriz = zeros(size,N);

APDmean = zeros(1,size);
APDVoIPmean = zeros(1,size);
APDterm = zeros(1,size);
APDVoIPterm = zeros(1,size);

for j = 1:size %para cada valor de n
   for i = 1:N %chamamos o simulator3 50vezes para cada valor de n
       [PLmatriz(j,i), PLVoIPmatriz(j,i), APDmatriz(j,i), APDVoIPmatriz(j,i), MPDmatriz(j,i), MPDVoIPmatriz(j,i), TTmatriz(j,i)] = Simulator4(lambda,C,f,P,n(j)); %n(0)=100000, n(1)=20000, ...
   end
   
   APDmean(j) = mean(APDmatriz(j,:)); %media dos valores de cada valor de n (n(0) fica no APDmean(0), n(1) fica no APDmean(1)...
   APDVoIPmean(j) = mean(APDVoIPmatriz(j,:));
   APDterm(j) = norminv(1-alfa/2)*sqrt(var(APDmatriz(j,:))/N);
   APDVoIPterm(j) = norminv(1-alfa/2)*sqrt(var(APDVoIPmatriz(j,:))/N);
   
end

%grafico do packet delay data
figure(1)
bar(n, APDmean)
ylabel('Average packet delay data (ms)')
xlabel('number of VoIP flows')
ylim([0 8])
hold on
er = errorbar(n, APDmean, APDterm);
er.Color = [0 1 0];
er.LineStyle = 'none';
title('Avg Packet Delay data (ms)');
hold off

%grafico do packet delay voip
figure(2)
bar(n, APDVoIPmean)
ylabel('Average packet delay VoIP (ms)')
xlabel('number of VoIP flows')
ylim([0 8])
hold on
er = errorbar(n, APDVoIPmean, APDVoIPterm);
er.Color = [0 1 0];
er.LineStyle = 'none';
title('Avg Packet Delay VoIP (ms)');
hold off

%% c

n = [40 30 20 10];
size = length(n);
WVoIPms = zeros(1,size); %atraso VoIP em ms
WDatams = zeros(1,size); %atraso Data em ms

for i = 1:size
    lambdaData = 1500;
    lambdaVoIP = 50*n(i); % 1/0.02
    avgVoIPPacketSize = 120; %media de 110 a 130, todos os tamanhos tem a mesma prob
    avgDataPacketSize = 620; %valor calculado no ex4 do guiao1
    miuVoIP = 10e6/(120*8); %pps
    miuData = 10e6/(620*8); %pps
    ESVoIP = 1/miuVoIP; %seg
    ESData = 1/miuData; %seg
    ES2VoIP = 2/(miuVoIP^2); %seg2
    ES2Data = 2/(miuData^2); %seg2
    roVoIP = lambdaVoIP*ESVoIP; 
    roData = lambdaData*ESData;
    WVoIP = (((lambdaVoIP*ES2VoIP)+(lambdaData*ES2Data))/(2*(1-roVoIP))) + ESVoIP;
    WData = (((lambdaVoIP*ES2VoIP)+(lambdaData*ES2Data))/(2*(1-roVoIP)*(1-roVoIP-roData))) + ESData;
    WVoIPms(i) = WVoIP*1000;
    WDatams(i) = WData*1000;
end

%grafico do packet delay VoIP
figure(1)
bar(n, WVoIPms)
ylabel('Average packet delay VoIP (ms)')
xlabel('number of VoIP flows')
ylim([0 9])
hold on
title('Avg VoIP Packet Delay (ms) - theoretical values');
hold off

%grafico do packet delay Data
figure(2)
bar(n, WDatams)
ylabel('Average packet delay Data (ms)')
xlabel('number of Data flows')
ylim([0 9])
hold on
title('Avg Data Packet Delay (ms) - theoretical values');
hold off

%% d
N=50; %running simulator N times 
per1 = zeros(4,N);
per2 = zeros(4,N);
per3 = zeros(4,N);
per4 = zeros(4,N);
per5 = zeros(4,N);
per6 = zeros(4,N);
per7 = zeros(4,N);

lambda= 1500;
C=10;
f=10000;
P=10000;
n=[10 20 30 40];
for i =1:4
    for it = 1 : N
        [per1(i,it), per2(i,it), per3(i,it), per4(i,it), per5(i,it), per6(i,it), per7(i,it)] = Simulator3(lambda,C,f,P,n(i));
    end
end
alfa = 0.1;
term1=zeros(1,4);
term2=zeros(1,4);
term3=zeros(1,4);
term4=zeros(1,4);

avgpacketDelayData= zeros(1,4);
avgpacketDelayVoip= zeros(1,4);
avgPacketLossData= zeros(1,4);
avgPacketLossVOIP= zeros(1:4);


y = [];     %storing average packet delay values
y1 = [];    %storing packet loss values
errorAvgDelay = []; %storing avg packet delay errors
errorPacketLoss = []; %storing packet loss errors
for i =1:4
    avgpacketDelayData(i) = mean(per3(i, 1:N));
    term1(i) = norminv(1-alfa/2)*sqrt(var(per3(i, 1:N))/N);
    fprintf('Av. Packet Delay (ms) of DATA with     n= %d   => %.2e +- %.2e  &',n(i), avgpacketDelayData(i),term1(i))
    
    avgpacketDelayVoip(i) = mean(per4(i, 1:N));
    term2(i) = norminv(1-alfa/2)*sqrt(var(per4(i, 1:N))/N);
    fprintf('  Av. Packet Delay (ms) of VOIP with   n= %d   => %.2e +- %.2e\n',n(i), avgpacketDelayVoip(i),term2(i))

    avgPacketLossData(i) = mean(per1(i, 1:N));
    term3(i) = norminv(1-alfa/2)*sqrt(var(per1(i, 1:N))/N);
    fprintf('PacketLoss of DATA with                n= %d   => %.2e +- %.2e  &  ',n(i), avgPacketLossData(i),term3(i))
    
    avgPacketLossVOIP(i) = mean(per2(i, 1:N));
    term4(i) = norminv(1-alfa/2)*sqrt(var(per2(i, 1:N))/N);
    fprintf('PacketLoss of voip with              n= %d   => %.2e +- %.2e\n',n(i), avgPacketLossVOIP(i),term4(i))
    
    errorAvgDelay =[errorAvgDelay; term1(i) term2(i)];
    errorPacketLoss =[errorPacketLoss; term3(i) term4(i)];

    y = [y; avgpacketDelayData(i) avgpacketDelayVoip(i)];
    y1 = [y1; avgPacketLossData(i) avgPacketLossVOIP(i)];
    
end

%figure("Average packet delay for data")
figure(1)
b = bar(n,y, 'grouped');
hold on
%determining error
[ngroups,nbars] = size(y);
x = nan(nbars, ngroups);

for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

errorbar(x',y,errorAvgDelay,'k','linestyle','none');
    
title('Average packet delays for Data and VOIP packets')
ylabel('Average packet delay (ms)')
xlabel('Number of VOIP packets flows')
legend('Data Packets', 'VOIP packets', Location='northwest')
grid on
hold off

%2nd figure for VOIP avg packet delay 
figure(2)
b1 = bar(n,y1, 'grouped');
hold on
%determining error
[ngroups,nbars] = size(y1);
x1 = nan(nbars, ngroups);

for i = 1:nbars
    x1(i,:) = b1(i).XEndPoints;
end

errorbar(x1',y1,errorPacketLoss,'k','linestyle','none');
    
title('Packet Loss for Data and VOIP packets')
ylabel('Packet Loss (%)')
xlabel('Number of VOIP packets flows')
legend('Data Packets', 'VOIP packets', Location='northwest')
grid on
hold off

%% e)

N=50; %running simulator N times 
per1 = zeros(4,N);
per2 = zeros(4,N);
per3 = zeros(4,N);
per4 = zeros(4,N);
per5 = zeros(4,N);
per6 = zeros(4,N);
per7 = zeros(4,N);

lambda= 1500;
C=10;
f=10000;
P=10000;
n=[10 20 30 40];
for i =1:4
    for it = 1 : N
        [per1(i,it), per2(i,it), per3(i,it), per4(i,it), per5(i,it), per6(i,it), per7(i,it)] = Simulator4(lambda,C,f,P,n(i));
    end
end
alfa = 0.1;
term1=zeros(1,4);
term2=zeros(1,4);
term3=zeros(1,4);
term4=zeros(1,4);

avgpacketDelayData= zeros(1,4);
avgpacketDelayVoip= zeros(1,4);
avgPacketLossData= zeros(1,4);
avgPacketLossVOIP= zeros(1:4);


y = [];     %storing average packet delay values
y1 = [];    %storing packet loss values
errorAvgDelay = []; %storing avg packet delay errors
errorPacketLoss = []; %storing packet loss errors
for i =1:4
    avgpacketDelayData(i) = mean(per3(i, 1:N));
    term1(i) = norminv(1-alfa/2)*sqrt(var(per3(i, 1:N))/N);
    fprintf('Av. Packet Delay (ms) of DATA with     n= %d   => %.2e +- %.2e  &',n(i), avgpacketDelayData(i),term1(i))
    
    avgpacketDelayVoip(i) = mean(per4(i, 1:N));
    term2(i) = norminv(1-alfa/2)*sqrt(var(per4(i, 1:N))/N);
    fprintf('  Av. Packet Delay (ms) of VOIP with   n= %d   => %.2e +- %.2e\n',n(i), avgpacketDelayVoip(i),term2(i))

    avgPacketLossData(i) = mean(per1(i, 1:N));
    term3(i) = norminv(1-alfa/2)*sqrt(var(per1(i, 1:N))/N);
    fprintf('PacketLoss of DATA with                n= %d   => %.2e +- %.2e  &  ',n(i), avgPacketLossData(i),term3(i))
    
    avgPacketLossVOIP(i) = mean(per2(i, 1:N));
    term4(i) = norminv(1-alfa/2)*sqrt(var(per2(i, 1:N))/N);
    fprintf('PacketLoss of voip with              n= %d   => %.2e +- %.2e\n',n(i), avgPacketLossVOIP(i),term4(i))
    
    errorAvgDelay =[errorAvgDelay; term1(i) term2(i)];
    errorPacketLoss =[errorPacketLoss; term3(i) term4(i)];

    y = [y; avgpacketDelayData(i) avgpacketDelayVoip(i)];
    y1 = [y1; avgPacketLossData(i) avgPacketLossVOIP(i)];
    
end

%figure("Average packet delay for data")
figure(1)
b = bar(n,y, 'grouped');
hold on
%determining error
[ngroups,nbars] = size(y);
x = nan(nbars, ngroups);

for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

errorbar(x',y,errorAvgDelay,'k','linestyle','none');
    
title('Average packet delays for Data and VOIP packets')
ylabel('Average packet delay (ms)')
xlabel('Number of VOIP packets flows')
legend('Data Packets', 'VOIP packets', Location='northwest')
grid on
hold off

%2nd figure for VOIP avg packet delay 
figure(2)
b1 = bar(n,y1, 'grouped');
hold on
%determining error
[ngroups,nbars] = size(y1);
x1 = nan(nbars, ngroups);

for i = 1:nbars
    x1(i,:) = b1(i).XEndPoints;
end

errorbar(x1',y1,errorPacketLoss,'k','linestyle','none');
    
title('Packet Loss for Data and VOIP packets')
ylabel('Packet Loss (%)')
xlabel('Number of VOIP packets flows')
legend('Data Packets', 'VOIP packets', Location='northwest')
grid on
hold off

%% f
