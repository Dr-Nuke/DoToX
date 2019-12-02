function [P] = f42_ChannelPhanTransform_w_film(lft)

P=f4_ChannelPhan2_w_film(lft);

P4=f42_RotPhan(P,pi/2);

P5=f42_AddPhan(P,P4);

P6=f42_RotPhan(P5,pi);

P=f42_AddPhan(P5,P6);

end