function [P] = f42_ChannelPhanTransform()

P=f4_ChannelPhan2();

P4=f42_RotPhan(P,pi/2);

P5=f42_AddPhan(P,P4);

P6=f42_RotPhan(P5,pi);

P=f42_AddPhan(P5,P6);

end