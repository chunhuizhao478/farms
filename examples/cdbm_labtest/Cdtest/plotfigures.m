xiCd40e5 = readmatrix("xiCd40e5.txt");
xiCd40e4 = readmatrix("xiCd40e4.txt");
xiCd40e3 = readmatrix("xiCd40e3.txt");
xiCd40e2 = readmatrix("xiCd40e2.txt");
xiCd40e1 = readmatrix("xiCd40e1.txt");

alphaCd40e5 = readmatrix("alphaCd40e5.txt");
alphaCd40e4 = readmatrix("alphaCd40e4.txt");
alphaCd40e3 = readmatrix("alphaCd40e3.txt");
alphaCd40e2 = readmatrix("alphaCd40e2.txt");
alphaCd40e1 = readmatrix("alphaCd40e1.txt");

figure();
plot(xiCd40e5,alphaCd40e5,'r-'); hold on;
plot(xiCd40e4,alphaCd40e4,'b-'); hold on;
plot(xiCd40e3,alphaCd40e3,'m-'); hold on;
plot(xiCd40e2,alphaCd40e2,'g-'); hold on;
plot(xiCd40e1,alphaCd40e1,'k-'); hold on;

legend("Cd=40e5,Cb=100Cd","Cd=40e4,Cb=100Cd","Cd=40e3,Cb=100Cd","Cd=40e2,Cb=100Cd","Cd=40e1,Cb=100Cd")

title("damage variable vs strain invariant ratio")
xlabel("strain invariant ratio")
ylabel("damage variable")

% xline(-0.3)

xiCd40e1_eta0dot8 = readmatrix("xiCd40e1eta0dot8.txt");
alphaCd40e1_eta0dot8 = readmatrix("alphaCd40e1eta0dot8.txt");
xiCd40e1_eta0dot1 = readmatrix("xiCd40e1eta0dot1.txt");
alphaCd40e1_eta0dot1 = readmatrix("alphaCd40e1eta0dot1.txt");

figure();
plot(xiCd40e1_eta0dot8,alphaCd40e1_eta0dot8,'r-'); hold on;
plot(xiCd40e1,alphaCd40e1,'k-'); hold on;
plot(xiCd40e1_eta0dot1,alphaCd40e1_eta0dot1,'b-'); hold on;

legend("Cd=40e1,Cb=100Cd,eta=0.8","Cd=40e4,Cb=100Cd,eta=0.5","Cd=40e3,Cb=100Cd,eta=0.1")

title("damage variable vs strain invariant ratio")
xlabel("strain invariant ratio")
ylabel("damage variable")

xiCd40e1_xi01dot1 = readmatrix("xiCd40e1xi01dot0.txt");
alphaCd40e1_xi01dot1 = readmatrix("alphaCd40e1xi01dot0.txt");

figure();
plot(xiCd40e1_xi01dot1,alphaCd40e1_xi01dot1,'r-'); hold on;
plot(xiCd40e1,alphaCd40e1,'k-'); hold on;
