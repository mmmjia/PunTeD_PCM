 month={'2020-08','2020-09','2020-10','2020-11','2020-12','2021-01','2021-02','2021-03','2021-04','2021-05','2021-06','2021-07','2021-08','2021-09','2021-10','2021-11','2021-12','2022-01','2022-02','2022-03','2022-04','2022-05','2022-06','2022-07','2022-08','2022-09','2022-10','2022-11','2022-12','2023-01','2023-02','2023-03','2023-04','2023-05','2023-06','2023-07','2023-08','2023-09','2023-10','2023-11','2023-12','2024-01','2024-02'};
ins_fre=[];
site_sele=[];
ins_en=[];
ms=0;
ncut=200;
%nsel=1000;

%nsel=500;
%[cov_var,pair]=xlsread(['ilbert\cov_hilbert\',month{tt-7},'cov_hilbert.xlsx'],'Sheet1');
%[cov_var,pair]=xlsread(['hilbert\cov_hilbert_aa\',month{tt-7},'cov_hilbert.xlsx'],'Sheet1');
[cov_var,pair]=xlsread(['hilbert\cov_hilbert_aa_hydrocalss_pair\',month{tt-7},'cov_hilbert.xlsx'],'Sheet1');

[~,nsel]=size(cov_var);

%nsel=707;

cov_var0=cov_var(:,nsel-ncut:nsel);
column_names = pair(nsel-ncut:nsel);
pos=1:size(column_names,2);

for site=1:size(column_names,2)
%for site=siteseta
%for site=66+36:200%66:66+36%1:66
     cov_evo=cov_var0(:,site);
     fo=emd(cov_evo(1:tt));
     if isempty(fo)==0
       [hs,f,t,imfinsf,imfinse]=hht(fo);
       
       marginal_spectrum = sum(hs, 2); 
       imfinse(:,imfinsf(tt,:)<0)=[];
       imfinsf(:,imfinsf(tt,:)<0)=[];
       %imfinsf(:,imfinse(tt,:)<0.0001)=[];
       %imfinse(:,imfinse(tt,:)<0.0001)=[];
       
       ins_en=[ins_en,imfinse(tt,:)];
       ins_fre=[ins_fre,imfinsf(tt,:)];
       [s,ss]=size(imfinsf(tt,:));
       s=pos(site)*ones([1,ss]);
       site_sele=[site_sele,s];
       ms=ms+marginal_spectrum;

     end
     
end
edges = linspace(0, 3, 100);
histogram(ins_fre,edges);
%xlim([0,3])
%xlabel('f (frequency)','FontSize', 15)
%ylabel('number of modes','FontSize', 15)
[f0,x]=ksdensity(ins_fre,'Bandwidth',0.08);
data=[ins_fre;ins_en]';

%data_to_save = [column_names(sele_p); num2cell(cov_var(tt, sele_p))];
%xlswrite(['C:\Users\user\OneDrive - HKUST Connect\hilbert\plotcov\delta_sele',num2str(tt),'.xlsx'],data_to_save,'Sheet1');
%data_to_save = [column_names(sele_p); num2cell(cov_var(tt, sele_p))];
%xlswrite(['C:\Users\user\OneDrive - HKUST Connect\hilbert\plotcov\delta_sele',num2str(tt),'.xlsx'],data_to_save,'Sheet1');
