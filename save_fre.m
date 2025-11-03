save_data={};
for tt=7:52
    cov_dos
    save_data(tt).x=x;
    save_data(tt).f0=f0;
    save_data(tt).ins_fre=ins_fre;
    save_data(tt).ins_en=ins_en;
    save_data(tt).site_sele=site_sele;
    save_data(tt).ms=ms;
end