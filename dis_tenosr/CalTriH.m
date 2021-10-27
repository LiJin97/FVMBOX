function volume = CalTriH(Pcc,Pdfl,Pdbl,Publ,Pufl)
volume = (norm(dot(cross(Pcc-Pdfl,Pcc-Pdbl),Pcc-Publ))+norm(dot(cross(Pcc-Pdfl,Pcc-Publ),Pcc-Pufl)))/6;
