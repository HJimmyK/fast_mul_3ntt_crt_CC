clc,clear

set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 12);

data1 = readmatrix("cc_ntt-crt_times.csv");
data2 = readmatrix("cpp_ntt-crt_times.csv");

len = linspace(1000,10000000,301);
figure(Color='w')

plot(len, data1,'Color','r','LineWidth',1.1)
hold on
plot(len, data2,'Color','b','LineWidth',1.1)
ylabel("time (us)", "FontName", "Times New Roman")
xlabel("u64 array len", "FontName", "Times New Roman")
title("total time",  "FontName", "Times New Roman", "FontSize", 12)

figure(Color='w')
y = data2./data1;
plot(len,y,"Color",'r','LineWidth',1)
ylabel("speed up", "FontName", "Times New Roman")
xlabel("u64 array len", "FontName", "Times New Roman")
title("speed up",  "FontName", "Times New Roman", "FontSize", 12)

