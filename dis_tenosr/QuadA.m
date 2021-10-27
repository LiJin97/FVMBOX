function volume = QuadA(Pdfl,Pdfr,Pdbl,Pdbr,Pufl,Pufr,Publ,Pubr,Pcc)
volume = CalTriH(Pcc,Pdfl,Pdbl,Publ,Pufl);
volume = volume+CalTriH(Pcc,Pdfr,Pdbr,Pubr,Pufr);
volume = volume+CalTriH(Pcc,Pdfl,Pdfr,Pufr,Pufl);
volume = volume+CalTriH(Pcc,Pdbl,Pdbr,Pubr,Publ);
volume = volume+CalTriH(Pcc,Pufl,Pufr,Pubr,Publ);
volume = volume+CalTriH(Pcc,Pdfl,Pdfr,Pdbr,Pdbl);

