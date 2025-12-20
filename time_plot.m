clc,clear

set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 12);

data1 = readmatrix("data\time_ntt_cpp.csv");
data2 = readmatrix("data\time_ntt_c.csv");
data1=data1./1000;

len = linspace(1000,10000000,251);
figure(Color='w')

subplot(3,1,1)
plot(len, data1(:,1),'Color','r','LineWidth',1.1)
hold on
plot(len, data2(:,1),'Color','b','LineWidth',1.1)
ylabel("time (ms)", "FontName", "Times New Roman")
title("convolution1 time",  "FontName", "Times New Roman", "FontSize", 12)
subplot(3,1,2)
plot(len, data1(:,2),'Color','r','LineWidth',1.1)
hold on
plot(len, data2(:,2),'Color','b','LineWidth',1.1)
ylabel("time (ms)", "FontName", "Times New Roman")
title("convolution2 time",  "FontName", "Times New Roman", "FontSize", 12)
subplot(3,1,3)
plot(len, data1(:,3),'Color','r','LineWidth',1.1)
hold on
plot(len, data2(:,3),'Color','b','LineWidth',1.1)
ylabel("time (ms)", "FontName", "Times New Roman")
xlabel("u64 array len", "FontName", "Times New Roman")
title("convolution3 time",  "FontName", "Times New Roman", "FontSize", 12)

figure(Color='w')
plot(len, data1(:,4),'Color','r','LineWidth',1.1)
hold on
plot(len, data2(:,4),'Color','b','LineWidth',1.1)
ylabel("time (ms)", "FontName", "Times New Roman")
xlabel("u64 array len", "FontName", "Times New Roman")
title("crt time",  "FontName", "Times New Roman", "FontSize", 12)

figure(Color='w')
plot(len, data1(:,5),'Color','r','LineWidth',1.1)
hold on
plot(len, data2(:,5),'Color','b','LineWidth',1.1)
ylabel("time (ms)", "FontName", "Times New Roman")
xlabel("u64 array len", "FontName", "Times New Roman")
title("total time",  "FontName", "Times New Roman", "FontSize", 12)
