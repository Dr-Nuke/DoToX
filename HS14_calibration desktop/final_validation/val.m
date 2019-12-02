classdef val % data validation toolbox
methods(Static)

function [energy,att]=get_xs()   % prepares the xs data
    % energy in keV
    % att in 1/cm

    xsData=[ % energy MeV, chloro cm2/g, water, alu
        1.100E-02, 3.905E+01, 4.026E+00, 1.981E+01;
        1.200E-02, 3.043E+01, 3.126E+00, 1.534E+01;
        1.300E-02, 2.417E+01, 2.487E+00, 1.212E+01;
        1.400E-02, 1.951E+01, 2.021E+00, 9.744E+00;
        1.500E-02, 1.598E+01, 1.672E+00, 7.955E+00;
        1.600E-02, 1.326E+01, 1.408E+00, 6.583E+00;
        1.700E-02, 1.112E+01, 1.203E+00, 5.513E+00;
        1.800E-02, 9.421E+00, 1.042E+00, 4.667E+00;
        1.900E-02, 8.055E+00, 9.134E-01, 3.990E+00;
        2.000E-02, 6.944E+00, 8.098E-01, 3.442E+00;
        2.100E-02, 6.031E+00, 7.253E-01, 2.993E+00;
        2.200E-02, 5.274E+00, 6.557E-01, 2.623E+00;
        2.300E-02, 4.642E+00, 5.978E-01, 2.314E+00;
        2.400E-02, 4.109E+00, 5.493E-01, 2.055E+00;
        2.500E-02, 3.657E+00, 5.082E-01, 1.836E+00;
        2.600E-02, 3.272E+00, 4.733E-01, 1.649E+00;
        2.700E-02, 2.941E+00, 4.433E-01, 1.489E+00;
        2.800E-02, 2.655E+00, 4.175E-01, 1.351E+00;
        2.900E-02, 2.407E+00, 3.951E-01, 1.232E+00;
        3.000E-02, 2.190E+00, 3.756E-01, 1.128E+00;
        3.100E-02, 2.001E+00, 3.584E-01, 1.037E+00;
        3.200E-02, 1.834E+00, 3.433E-01, 9.576E-01;
        3.300E-02, 1.687E+00, 3.300E-01, 8.873E-01;
        3.400E-02, 1.557E+00, 3.181E-01, 8.250E-01;
        3.500E-02, 1.441E+00, 3.075E-01, 7.696E-01;
        3.600E-02, 1.337E+00, 2.980E-01, 7.202E-01;
        3.700E-02, 1.244E+00, 2.894E-01, 6.760E-01;
        3.800E-02, 1.161E+00, 2.817E-01, 6.364E-01;
        3.900E-02, 1.086E+00, 2.747E-01, 6.007E-01;
        4.000E-02, 1.019E+00, 2.683E-01, 5.684E-01;
        4.100E-02, 9.574E-01, 2.625E-01, 5.392E-01;
        4.200E-02, 9.019E-01, 2.571E-01, 5.128E-01;
        4.300E-02, 8.514E-01, 2.523E-01, 4.886E-01;
        4.400E-02, 8.053E-01, 2.478E-01, 4.667E-01;
        4.500E-02, 7.632E-01, 2.436E-01, 4.466E-01;
        4.600E-02, 7.247E-01, 2.398E-01, 4.281E-01;
        4.700E-02, 6.893E-01, 2.362E-01, 4.112E-01;
        4.800E-02, 6.568E-01, 2.329E-01, 3.957E-01;
        4.900E-02, 6.269E-01, 2.298E-01, 3.814E-01;
        5.000E-02, 5.993E-01, 2.269E-01, 3.681E-01;
        5.100E-02, 5.738E-01, 2.242E-01, 3.559E-01;
        5.200E-02, 5.501E-01, 2.217E-01, 3.446E-01;
        5.300E-02, 5.282E-01, 2.193E-01, 3.340E-01;
        5.400E-02, 5.079E-01, 2.171E-01, 3.242E-01;
        5.500E-02, 4.890E-01, 2.149E-01, 3.151E-01;
        5.600E-02, 4.714E-01, 2.129E-01, 3.066E-01;
        5.700E-02, 4.550E-01, 2.110E-01, 2.987E-01;
        5.800E-02, 4.396E-01, 2.092E-01, 2.913E-01;
        5.900E-02, 4.253E-01, 2.075E-01, 2.843E-01;
        6.000E-02, 4.119E-01, 2.059E-01, 2.778E-01;
        6.100E-02, 3.993E-01, 2.043E-01, 2.717E-01;
        6.200E-02, 3.876E-01, 2.028E-01, 2.659E-01;
        6.300E-02, 3.765E-01, 2.014E-01, 2.605E-01;
        6.400E-02, 3.661E-01, 2.000E-01, 2.554E-01;
        6.500E-02, 3.563E-01, 1.987E-01, 2.506E-01;
        6.600E-02, 3.470E-01, 1.974E-01, 2.460E-01;
        6.700E-02, 3.383E-01, 1.962E-01, 2.417E-01;
        6.800E-02, 3.301E-01, 1.951E-01, 2.376E-01;
        6.900E-02, 3.223E-01, 1.939E-01, 2.338E-01;
        7.000E-02, 3.149E-01, 1.929E-01, 2.301E-01;
        7.100E-02, 3.080E-01, 1.918E-01, 2.266E-01;
        7.200E-02, 3.014E-01, 1.908E-01, 2.233E-01;
        7.300E-02, 2.951E-01, 1.898E-01, 2.202E-01;
        7.400E-02, 2.891E-01, 1.888E-01, 2.172E-01;
        7.500E-02, 2.835E-01, 1.879E-01, 2.143E-01;
        7.600E-02, 2.781E-01, 1.870E-01, 2.116E-01;
        7.700E-02, 2.730E-01, 1.861E-01, 2.089E-01;
        7.800E-02, 2.681E-01, 1.853E-01, 2.064E-01;
        7.900E-02, 2.635E-01, 1.845E-01, 2.041E-01;
        8.000E-02, 2.590E-01, 1.837E-01, 2.018E-01;
        8.100E-02, 2.548E-01, 1.829E-01, 1.996E-01;
        8.200E-02, 2.507E-01, 1.821E-01, 1.975E-01;
        8.300E-02, 2.469E-01, 1.814E-01, 1.955E-01;
        8.400E-02, 2.431E-01, 1.806E-01, 1.935E-01;
        8.500E-02, 2.396E-01, 1.799E-01, 1.917E-01;
        8.600E-02, 2.362E-01, 1.792E-01, 1.899E-01;
        8.700E-02, 2.329E-01, 1.785E-01, 1.881E-01;
        8.800E-02, 2.298E-01, 1.779E-01, 1.865E-01;
        8.900E-02, 2.268E-01, 1.772E-01, 1.849E-01;
        9.000E-02, 2.239E-01, 1.766E-01, 1.833E-01;
        9.100E-02, 2.211E-01, 1.759E-01, 1.818E-01;
        9.200E-02, 2.185E-01, 1.753E-01, 1.804E-01;
        9.300E-02, 2.159E-01, 1.747E-01, 1.790E-01;
        9.400E-02, 2.134E-01, 1.741E-01, 1.777E-01;
        9.500E-02, 2.111E-01, 1.735E-01, 1.764E-01;
        9.600E-02, 2.088E-01, 1.729E-01, 1.751E-01;
        9.700E-02, 2.066E-01, 1.724E-01, 1.739E-01;
        9.800E-02, 2.044E-01, 1.718E-01, 1.727E-01;
        9.900E-02, 2.024E-01, 1.713E-01, 1.715E-01;
        1.000E-01, 2.004E-01, 1.707E-01, 1.704E-01;
        1.010E-01, 1.984E-01, 1.702E-01, 1.693E-01;
        1.020E-01, 1.966E-01, 1.697E-01, 1.683E-01;
        1.030E-01, 1.948E-01, 1.692E-01, 1.673E-01;
        1.040E-01, 1.930E-01, 1.686E-01, 1.663E-01;
        1.050E-01, 1.914E-01, 1.681E-01, 1.653E-01;
        1.060E-01, 1.897E-01, 1.677E-01, 1.644E-01;
        1.070E-01, 1.881E-01, 1.672E-01, 1.634E-01;
        1.080E-01, 1.866E-01, 1.667E-01, 1.625E-01;
        1.090E-01, 1.851E-01, 1.662E-01, 1.617E-01;
        1.100E-01, 1.837E-01, 1.657E-01, 1.608E-01;];

    densities=[1.411, 0.9555, 2.70]; %g/cm3, water: 2bar/ 104 °C
    energy=xsData(:,1)*1000; % energies in keV, steps of  1 keV from 11 to 110


    [rows,~] = size(xsData);
    att = xsData(:,2:end).*repmat(densities,rows,1);

    % seprate into different variables
    % macroXSenergies= xsData(:,1)*1000; % energies in keV, steps of  1 keV from 11 to 110
    % macroXSchloroform = macroXSall(:,1);
    % macroXSwater = macroXSall(:,2);
    % macroXSaluminum = macroXSall(:,3);

end

function [spec,cufilts]=get_spectrum()


    % peak energy 110 keV
    % minimum energy 11
    % 1 keV bins
    % thickness air 1000 mm
    % other thickness 0 except:
    % thicknes alu 3.2 mm
    % otherwise default SpekCalc data
    % units #/keV/cm2/mAs @ 1 meter distance
    % Energy[keV]  N[keV cm^2 mAs]^-1 @ 1 meter
    % Cu thickness	0.0 0.2 0.5 1.0 1.5 2.0 2.5 3.0
    cufilts = [0.0 0.2 0.3 0.4 0.5 1.0 1.5 2.0 2.5 3.0];
    sourceDataAll = [
        11	0.4535152	2.67E-14	6.47E-21	1.57E-27	3.81E-34	3.20E-67	2.69E-100	2.25E-133	1.89E-166	1.59E-199
        12	20.21458	5.26E-10	2.69E-15	1.37E-20	7.00E-26	2.42E-52	8.39E-79	2.90E-105	1.01E-131	3.48E-158
        13	324.3366	1.07E-06	6.12E-11	3.51E-15	2.01E-19	1.25E-40	7.74E-62	4.80E-83	2.98E-104	1.85E-125
        14	2743.748	0.0003025	1.00E-07	3.34E-11	1.11E-14	4.47E-32	1.80E-49	7.28E-67	2.94E-84	1.19E-101
        15	13652.11	0.023352	0.0000305	3.99E-08	5.22E-11	2.00E-25	7.65E-40	2.93E-54	1.12E-68	4.29E-83
        16	47394.28	0.6962384	0.0026685	0.0000102	3.92E-08	3.24E-20	2.68E-32	2.22E-44	1.83E-56	1.52E-68
        17	123253.3	9.758678	0.0868335	0.0007727	6.88E-06	3.83E-16	2.14E-26	1.19E-36	6.66E-47	3.71E-57
        18	259047.7	80.06569	1.407602	0.0247465	0.0004351	7.31E-13	1.23E-21	2.06E-30	3.46E-39	5.81E-48
        19	468564.9	439.9011	13.47869	0.4129907	0.0126541	3.42E-10	9.23E-18	2.49E-25	6.73E-33	1.82E-40
        20	751746.2	1760.22		85.17551	4.121569	0.1994392	5.29E-08	1.40E-14	3.72E-21	9.88E-28	2.62E-34
        21	1.10E+06	5454.793	384.5562	27.11074	1.911274	3.33E-06	5.80E-12	1.01E-17	1.76E-23	3.06E-29
        22	1.49E+06	14099.61	1372.381	133.5802	13.00199	0.0001136	9.92E-10	8.67E-15	7.57E-20	6.62E-25
        23	1.90E+06	30881.43	3932.814	500.8519	63.78451	0.0021367	7.16E-08	2.40E-12	8.03E-17	2.69E-21
        24	2.33E+06	60087.52	9659.697	1552.897	249.6444	0.0268051	2.88E-06	3.09E-10	3.32E-14	3.56E-18
        25	2.71E+06	103916.4	20345.63	3983.442	779.9124	0.2243791	0.0000646	1.86E-08	5.34E-12	1.54E-15
        26	3.10E+06	167005.5	38766.11	8998.577	2088.793	1.407678	0.0009487	6.39E-07	4.31E-10	2.90E-13
        27	3.45E+06	247410.1	66282.6		17757.49	4757.334	6.565592	0.0090612	0.0000125	1.73E-08	2.38E-11
        28	3.76E+06	353108		108207.9	33159.66	10161.58	27.46107	0.0742119	0.0002006	5.42E-07	1.46E-09
        29	4.03E+06	469184.2	160099.5	54630.67	18641.6		86.24143	0.3989779	0.0018458	8.54E-06	3.95E-08
        30	4.24E+06	600913.6	226288.6	85214.46	32089.57	243.0061	1.840222	0.0139355	0.0001055	7.99E-07
        31	4.42E+06	742921.5	304621.6	124904.6	51214.89	593.5881	6.879773	0.0797376	0.0009242	0.0000107
        32	4.57E+06	894316.6	395717.2	175097		77476.9		1314.137	22.28993	0.3780742	0.0064128	0.0001088
        33	4.68E+06	1.05E+06	495927.2	234690.9	111064.3	2636.123	62.56865	1.485074	0.0352484	0.0008366
        34	4.76E+06	1.20E+06	605347.3	304469.3	153137.8	4929.208	158.6616	5.107008	0.1643846	0.0052912
        35	4.80E+06	1.35E+06	718245.7	381205.9	202323.5	8520.778	358.8495	15.11281	0.6364706	0.0268047
        36	4.83E+06	1.50E+06	834478.1	464851.7	258948.8	13890.3		745.0906	39.96747	2.143898	0.115001
        37	4.83E+06	1.63E+06	949742.8	552310.9	321189.5	21362.37	1420.816	94.49874	6.28513		0.4180252
        38	4.81E+06	1.76E+06	1.07E+06	6.45E+05	389984.8	31589.6		2558.824	207.2702	16.78932	1.359971
        39	4.78E+06	1.87E+06	1.17E+06	7.36E+05	460895.2	44460.93	4288.988	413.7434	39.91236	3.850204
        40	4.73E+06	1.98E+06	1.28E+06	8.30E+05	536667.2	60830.92	6895.149	781.5612	88.58951	10.04157
        41	4.68E+06	2.07E+06	1.38E+06	9.19E+05	611696.8	80023.21	10468.77	1369.542	179.1659	23.43878
        42	4.61E+06	2.16E+06	1.47E+06	1.01E+06	689540.4	103183.1	15440.37	2310.503	345.7447	51.73738
        43	4.53E+06	2.22E+06	1.56E+06	1.09E+06	765273.9	129240.3	21826.26	3686.044	622.5034	105.1291
        44	4.45E+06	2.29E+06	1.64E+06	1.17E+06	840992.5	158859.9	30007.96	5668.377	1070.732	202.2568
        45	4.36E+06	2.33E+06	1.70E+06	1.24E+06	909763.1	189650.7	39534.91	8241.515	1718.04		358.1455
        46	4.28E+06	2.37E+06	1.77E+06	1.32E+06	979319.5	224289.6	51368.13	11764.63	2694.407	617.089
        47	4.18E+06	2.40E+06	1.82E+06	1.38E+06	1.04E+06	260190.2	64884.77	16180.6		4035.027	1006.232
        48	4.09E+06	2.42E+06	1.86E+06	1.44E+06	1.10E+06	298612.7	80720.09	21820.01	5898.319	1594.416
        49	3.99E+06	2.44E+06	1.90E+06	1.49E+06	1.16E+06	336771.4	97799.98	28401.57	8247.948	2395.242
        50	3.89E+06	2.44E+06	1.93E+06	1.53E+06	1.21E+06	375729.3	116697.9	36245.27	11257.44	3496.453
        51	3.80E+06	2.44E+06	1.96E+06	1.57E+06	1.26E+06	415303.7	137338.7	45417.16	15019.21	4966.771
        52	3.70E+06	2.43E+06	1.97E+06	1.59E+06	1.29E+06	450357.1	157156.3	54841.15	19137.33	6678.149
        53	3.60E+06	2.42E+06	1.98E+06	1.62E+06	1.33E+06	488647		179935.9	66258.28	24398.47	8984.316
        54	3.51E+06	2.40E+06	1.98E+06	1.64E+06	1.36E+06	524863.6	203035		78540.81	30382.25	11752.88
        55	3.41E+06	2.38E+06	1.99E+06	1.66E+06	1.39E+06	563465.7	228978.1	93050.88	37813.51	15366.45
        56	3.32E+06	2.35E+06	1.98E+06	1.67E+06	1.40E+06	593655.4	251172.3	106269.6	44962.06	19023.2
        57	3.22E+06	2.32E+06	1.97E+06	1.67E+06	1.42E+06	625507		275537.2	121374.7	53465.82	23551.81
        58	8.60E+06	6.30E+06	5.39E+06	4.61E+06	3.94E+06	1.81E+06	829623		380486.8	174501.2	80030.78
        59	1.26E+07	9.33E+06	8.04E+06	6.92E+06	5.96E+06	2.82E+06	1.33E+06	631449.7	298825.1	141414.9
        60	2.95E+06	2.22E+06	1.93E+06	1.67E+06	1.45E+06	710830.7	348666		171022.4	83887.31	41147.14
        61	2.87E+06	2.18E+06	1.91E+06	1.66E+06	1.45E+06	734383.7	371694.3	188126		95216.35	48191.94
        62	2.78E+06	2.14E+06	1.88E+06	1.65E+06	1.45E+06	752241.3	391105.4	203343.6	105722.4	54967.23
        63	2.70E+06	2.10E+06	1.85E+06	1.63E+06	1.44E+06	769619.2	411042.2	219531.5	117248.5	62620.69
        64	2.62E+06	2.06E+06	1.83E+06	1.62E+06	1.44E+06	788359.6	432522.5	237297.4	130189.9	71426.88
        65	2.54E+06	2.01E+06	1.79E+06	1.60E+06	1.42E+06	798485.1	448001.4	251357.5	141027.7	79125.58
        66	2.46E+06	1.97E+06	1.76E+06	1.58E+06	1.41E+06	809799.2	464641.6	266599.2	152967.6	87768.82
        67	5.56E+06	4.49E+06	4.04E+06	3.63E+06	3.26E+06	1.92E+06	1.12E+06	659369.4	386899.1	227021.3
        68	2.31E+06	1.88E+06	1.69E+06	1.53E+06	1.38E+06	823474		491927.8	293868.4	175551.4	104871.1
        69	3.06E+06	2.52E+06	2.28E+06	2.07E+06	1.87E+06	1.14E+06	698794.9	426903.2	260800.9	159326.8
        70	2.05E+06	1.69E+06	1.54E+06	1.40E+06	1.27E+06	792461.5	492879.1	306550.9	190662.3	118584.3
        71	1.99E+06	1.65E+06	1.51E+06	1.37E+06	1.25E+06	789347.5	497585.1	313665.4	197726.9	124642.1
        72	1.93E+06	1.61E+06	1.48E+06	1.35E+06	1.24E+06	792480.1	508364.8	326108.9	209194.3	134195.2
        73	1.87E+06	1.57E+06	1.44E+06	1.32E+06	1.22E+06	791997.3	515853.1	335991.5	218842		142538.8
        74	1.81E+06	1.53E+06	1.41E+06	1.30E+06	1.19E+06	789229.7	521473		344556.3	227660.9	150423.9
        75	1.75E+06	1.49E+06	1.38E+06	1.27E+06	1.17E+06	785900		526534.9	352766.2	236345.2	158345.9
        76	1.69E+06	1.45E+06	1.34E+06	1.24E+06	1.15E+06	779692.9	529207.2	359193.1	243798		165475
        77	1.64E+06	1.41E+06	1.31E+06	1.21E+06	1.12E+06	772703.2	530845.5	364689.8	250541.1	172121.2
        78	1.58E+06	1.37E+06	1.27E+06	1.18E+06	1.10E+06	764619.4	531446.3	369380.1	256736.5	178443.9
        79	1.53E+06	1.33E+06	1.24E+06	1.15E+06	1.07E+06	754967.3	530647.8	372979.1	262157.8	184264.2
        80	1.47E+06	1.29E+06	1.20E+06	1.12E+06	1.05E+06	744324.8	528822.8	375714.5	266935.1	189650.2
        81	1.42E+06	1.25E+06	1.17E+06	1.09E+06	1.02E+06	732641		525913		377517.1	270993.8	194527.9
        82	1.37E+06	1.21E+06	1.13E+06	1.06E+06	993797.1	720122.3	521812.9	378114.5	273988.2	198536.5
        83	1.32E+06	1.17E+06	1.09E+06	1.03E+06	965897.5	706521.6	516796.8	378019.5	276508.5	202256.7
        84	1.27E+06	1.12E+06	1.06E+06	9.96E+05	936808.8	691101.8	509839.1	376118.1	277469.5	204694.6
        85	1.22E+06	1.08E+06	1.02E+06	9.63E+05	907958.1	675544		502621.9	373963.5	278238.4	207016.5
        86	1.17E+06	1.05E+06	986750.6	931591.2	879515.2	659680		494792.7	371119		278357.7	208782.1
        87	1.12E+06	1.01E+06	950739		898961		850002.8	642418.3	485529.3	366955.1	277338.7	209608.1
        88	1.08E+06	966678.6	915343.7	866735		820707.6	624739.6	475564.7	362009.6	275569.2	209769
        89	1.03E+06	926944		878978.4	833494.9	790364.9	605970.2	464595.3	356203.7	273100.2	209385
        90	985566.8	888593		843745.1	801160.8	760725.7	587178.5	453223.3	349827.7	270020.2	208419.5
        91	939706.2	849373		807516.9	767723.5	729891		566922.9	440341.8	342023.4	265657.3	206342.1
        92	895344.8	811309		772296.8	735160.6	699810.1	546978.2	427523.4	334156.3	261179.8	204140.6
        93	850660.8	772616.7	736322.3	701732.8	668768.2	525768.7	413346.2	324962.4	255477.3	200849.8
        94	806319.7	734051.7	700384.2	668260.9	637610.9	504201.6	398706		315283.5	249315.8	197150.7
        95	763428.6	696501		665270.7	635440.7	606948.3	482541.8	383635		305001.2	242484.9	192782.7
        96	719539.9	657755.1	628881.7	601275.7	574881.5	459305.7	366965.5	293189.7	234246		187152.5
        97	676870.7	619970.7	593340.4	567854		543462.3	436348.1	350345.7	281294		225852.1	181337.6
        98	633161.4	581080		556668.6	533282.7	510879.3	412213.4	332602.9	268367.5	216537.8	174718
        99	590840.6	543310.3	520998.9	499603.7	479087		388470.9	314994.2	255415.2	207105.1	167932.5
        100	547777.5	504614.9	484326.3	464853.4	446163.4	363398.9	295987.5	241081.1	196359.9	159934.6
        101	505225.5	466167.1	447785.2	430128.2	413167.4	337883.4	276317		225968.8	184794.6	151122.8
        102	461758.7	426748.4	410251.6	394392.5	379146.5	311314.2	255617.6	209885.6	172335.4	141503.2
        103	418615.1	387500.3	372821.2	358698.2	345110.2	284512.1	234554.5	193368.9	159415.2	131423.4
        104	373229.3	346045.6	333205.5	320841.8	308936.9	255719.4	211669.2	175207.1	145026		120043.8
        105	326181.5	302858.3	291829.8	281202.8	270962.9	225092.1	186986.8	155332.2	129036.3	107192
        106	280408.8	260732.1	251417.7	242436.1	233775.4	194897.3	162484.8	135462.7	112934.6	94152.96
        107	234662.9	218509.3	210854.4	203467.6	196339.7	164275		137447		115000.2	96219.32	80505.55
        108	181053.8	168802.1	162990.7	157379.4	151961.3	127543.6	107049.4	89848.24	75411.07	63293.72
        109	88982.31	83065.11	80255.75	77541.4		74918.85	63078.1		53108.74	44715.02	37647.91	31697.75
        110	0	0	0	0	0	0	0	0	0	0];


    spec=sourceDataAll(:,2:end);
end

function [phan,Coord]=make_phantom(lft,pixwidth,res)
% script to produce bitmaps of the channel cross section
% lft = film thikcness in mm
% width = phantom image width in pixel 

% make an image!
% CAD parameters
r1=2;
r2=3;
r4=4.14;
r5=5.14;
r6=14.604;

if lft==0 % the lft=0 case is special
    zeroflag=1;
    lft=0.1;
else
    zeroflag=0;
end
    

nx=pixwidth;  % number of x-pixels
ny=pixwidth;  % number of y-pixels
% xcoord=linspace(-pixwidth*res/2,pixwidth*res/2,nx);
% ycoord=linspace(-pixwidth*res/2,pixwidth*res/2,ny);
xcoord=linspace(-(pixwidth-1)*res/2,(pixwidth-1)*res/2,nx);
ycoord=linspace(-(pixwidth-1)*res/2,(pixwidth-1)*res/2,ny);
[xx,yy]=ndgrid(xcoord,ycoord); % coordinates matrix
Coord=cat(3,xx,yy);
% generate phantom, focus on 1 quarter
Ph = f4_ChannelPhan2_w_film(lft);
Ph = f42_RotPhan(Ph,pi/4);

% air is where...:
air1=xx.^2+yy.^2>(r6)^2; % outer circle
air2=sqrt((xx-Ph.c(5,1)).^2+(yy-Ph.c(5,2)).^2)<r1; % small cirlce outside

s=4; %segment number
air31=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
air32=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s));
air33=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r2;
air3=air31 & air32 & air33;

s=6;
air4=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
air4(atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s)))=0;
air4(sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r2)=0;
% air4=air41 & air42 & air43;

air5=(yy)>Ph.b(18);
air5(xx<Ph.Lx2(18,1))=0;
%air5(yy<Ph.Lx1(18,2))=0;
air5(xx>Ph.Lx1(18,1))=0;

air=air1|air2|air3|air4|air5;

s=16;
d2o=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r4;
d2o((yy)>Ph.b(17))=0;
s=2;
d2o((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)-2*pi) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s)-2*pi)) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r1)=0;

s=1;
liq=zeros(nx,ny);

s=13;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lft))=1;

s=1;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lft))=1;

s=15;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)-2*pi) |...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<(Ph.r(s)+lft))=1;

%assemble mask
mask=air+2*d2o+4*liq;

% mirror that 1/8th to make a quarter
for i=ceil(nx/2):nx
    for j=ceil(nx/2):i
        mask(i,j)=mask(j,i);
        
    end
end
if mod(nx,2)==1
% new treatment for odd phantoms
% mirror quadrants from NE to NW
mask(1:floor(nx/2),ceil(nx/2)+0:end)=...
    flipud(mask(ceil(nx/2)+1:end,ceil(nx/2)+0:end));
% make full mask, i.e. mirror from N to S
mask(:,1:floor(nx/2))=...
    fliplr(mask(:,ceil(nx/2)+1:end));
elseif mod(nx,2)==0
% old treatment for even phantoms
% mirror quadrants
mask(1:ceil(nx/2),ceil(nx/2)+1:end)=...
    flipud(mask(ceil(nx/2)+1:end,ceil(nx/2)+1:end));
% make full mask
mask(:,1:ceil(nx/2))=...
    fliplr(mask(:,ceil(nx/2)+1:end));
else
    error('wiered phantom pixel size, neither odd nor even')
end

% get vapor area and set it to "3" by filling from the central pixel
[tempmask,n]=bwlabel(~im2bw(mask));
mask(tempmask==tempmask(ceil(nx/2),ceil(nx/2)))=3;

if zeroflag==1;
    mask(mask==4)=3;
end
phan=mask;


end

function [phan,Coord]=make_phantom_resize(lft,pixwidth,res)
% script to produce bitmaps of the channel cross section
% ame as make_phantom, but considers the downscaling of 11
% lft = film thikcness in mm
% width = phantom image width in pixel 

% make an image!
% CAD parameters
r1=2;
r2=3;
r4=4.14;
r5=5.14;
r6=14.604;

if lft==0 % the lft=0 case is special
    zeroflag=1;
    lft=0.1;
else
    zeroflag=0;
end
    

nx=pixwidth;  % number of x-pixels
ny=pixwidth;  % number of y-pixels
% for later downscaling of 11
xmax=pixwidth*res/2+5*res;

xcoord=linspace(-xmax,xmax,nx);
ycoord=linspace(-xmax,xmax,ny);
[xx,yy]=ndgrid(xcoord,ycoord); % coordinates matrix
Coord=cat(3,xx,yy);
% generate phantom, focus on 1 quarter
Ph = f4_ChannelPhan2_w_film(lft);
Ph = f42_RotPhan(Ph,pi/4);

% air is where...:
air1=xx.^2+yy.^2>(r6)^2; % outer circle
air2=sqrt((xx-Ph.c(5,1)).^2+(yy-Ph.c(5,2)).^2)<r1; % small cirlce outside

s=4; %segment number
air31=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
air32=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s));
air33=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r2;
air3=air31 & air32 & air33;

s=6;
air4=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
air4(atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s)))=0;
air4(sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r2)=0;
% air4=air41 & air42 & air43;

air5=(yy)>Ph.b(18);
air5(xx<Ph.Lx2(18,1))=0;
%air5(yy<Ph.Lx1(18,2))=0;
air5(xx>Ph.Lx1(18,1))=0;

air=air1|air2|air3|air4|air5;

s=16;
d2o=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r4;
d2o((yy)>Ph.b(17))=0;
s=2;
d2o((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)-2*pi) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s)-2*pi)) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r1)=0;

s=1;
liq=zeros(nx,ny);

s=13;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lft))=1;

s=1;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lft))=1;

s=15;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)-2*pi) |...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<(Ph.r(s)+lft))=1;

%assemble mask
mask=air+2*d2o+4*liq;

% mirror that 1/8th to make a quarter
for i=ceil(nx/2):nx
    for j=ceil(nx/2):i
        mask(i,j)=mask(j,i);
        
    end
end
if mod(nx,2)==1
% new treatment for odd phantoms
% mirror quadrants from NE to NW
mask(1:floor(nx/2),ceil(nx/2)+0:end)=...
    flipud(mask(ceil(nx/2)+1:end,ceil(nx/2)+0:end));
% make full mask, i.e. mirror from N to S
mask(:,1:floor(nx/2))=...
    fliplr(mask(:,ceil(nx/2)+1:end));
elseif mod(nx,2)==0
% old treatment for even phantoms
% mirror quadrants
mask(1:ceil(nx/2),ceil(nx/2)+1:end)=...
    flipud(mask(ceil(nx/2)+1:end,ceil(nx/2)+1:end));
% make full mask
mask(:,1:ceil(nx/2))=...
    fliplr(mask(:,ceil(nx/2)+1:end));
else
    error('wiered phantom pixel size, neither odd nor even')
end

% get vapor area and set it to "3" by filling from the central pixel
[tempmask,n]=bwlabel(~im2bw(mask));
mask(tempmask==tempmask(ceil(nx/2),ceil(nx/2)))=3;

if zeroflag==1;
    mask(mask==4)=3;
end
phan=mask;


end

function [phan,Coord]=make_phantom_bak(lft,domain,res)
% script to produce bitmaps of the channel cross section
% lft = film thikcness in mm
% width = phantom image width in pixel 

% make an image!
% CAD parameters
r1=2;
r2=3;
r4=4.14;
r5=5.14;
r6=14.604;

if lft==0 % the lft=0 case is special
    zeroflag=1;
    lft=0.1;
else
    zeroflag=0;
end
    

nx=round(domain/res)+1;  % number of x-pixels
ny=round(domain/res)+1;  % number of y-pixels
xcoord=linspace(-domain/2,domain/2,nx);
ycoord=linspace(-domain/2,domain/2,ny);
[xx,yy]=ndgrid(xcoord,ycoord); % coordinates matrix
Coord=cat(3,xx,yy);
% generate phantom, focus on 1 quarter
Ph = f4_ChannelPhan2_w_film(lft);
Ph = f42_RotPhan(Ph,pi/4);

% air is where...:
air1=xx.^2+yy.^2>(r6)^2; % outer circle
air2=sqrt((xx-Ph.c(5,1)).^2+(yy-Ph.c(5,2)).^2)<r1; % small cirlce outside

s=4; %segment number
air31=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
air32=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s));
air33=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r2;
air3=air31 & air32 & air33;

s=6;
air4=atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>Ph.p0(s);
air4(atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s)))=0;
air4(sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r2)=0;
% air4=air41 & air42 & air43;

air5=(yy)>Ph.b(18);
air5(xx<Ph.Lx2(18,1))=0;
%air5(yy<Ph.Lx1(18,2))=0;
air5(xx>Ph.Lx1(18,1))=0;

air=air1|air2|air3|air4|air5;

s=16;
d2o=sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<r4;
d2o((yy)>Ph.b(17))=0;
s=2;
d2o((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)-2*pi) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s)-2*pi)) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>r1)=0;

s=1;
liq=zeros(nx,ny);

s=13;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lft))=1;

s=1;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)) &...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>(Ph.r(s)-lft))=1;

s=15;
liq((atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))<(Ph.p0(s)-2*pi) |...
    atan2(yy-Ph.c(s,2),xx-Ph.c(s,1))>(Ph.p0(s)+Ph.dp(s))) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)>Ph.r(s) &...
    sqrt((xx-Ph.c(s,1)).^2+(yy-Ph.c(s,2)).^2)<(Ph.r(s)+lft))=1;

%assemble mask
mask=air+2*d2o+4*liq;

% mirror that 1/8th to make a quarter
for i=ceil(nx/2):nx
    for j=ceil(nx/2):i
        mask(i,j)=mask(j,i);
        
    end
end
% mirror quadrants (leave out middle pixel only in y)
mask(1:ceil(nx/2),ceil(nx/2):end)=...
    flipud(mask(ceil(nx/2):end,ceil(nx/2):end));
% make full mask
mask(:,1:ceil(nx/2))=...
    fliplr(mask(:,ceil(nx/2):end));

% get vapor area and set it to "3"
[tempmask,n]=bwlabel(~im2bw(mask));
mask(tempmask==tempmask(ceil(nx/2),ceil(nx/2)))=3;

if zeroflag==1;
    mask(mask==4)=3;
end
phan=mask;


end

function [phan,phanCoord]=makePhantomSeries(lfts,pixwidth,res)
    % lfts array of liquid fim thicknesses in mm
    % widh = width of phantom area in pixel
    % r = resolution or pixel size
    
    
    phan=zeros(pixwidth,pixwidth,length(lfts));
    phanCoord=zeros(pixwidth,pixwidth,2);
    
    for i = 1:length(lfts)
        [phan(:,:,i),phanCoord]=val.make_phantom(lfts(i),pixwidth,res);
    end
end

function [phan,phanCoord]=makePhantomSeries_bak(lfts,width,res)
    
    phan=zeros(width/res+1,width/res+1,length(lfts));
    phanCoord=zeros(width/res+1,width/res+1,2);
    
    for i = 1:length(lfts)
        [phan(:,:,i),phanCoord]=val.make_phantom(lfts(i),width,res);
    end
end

function plotxs(energy,att,spec,T)
globcol=get(groot,'DefaultAxesColorOrder'); %get colors
fh=figure(3249);clf
ax=gca;

%yyaxis left
ax.YScale='log';
hold on
grid on
name={'CHCl3','H2O','Alu'};
colind=[2,1,3];

for i =1:size(att,2) % attenuations
    p=semilogy(energy,att(:,i),'Displayname',name{i},'color',globcol(colind(i),:));
    pub.BFLine(p,T.F)
end
ylh=ylabel('attenuation [1/cm]');
ylh.Color=[0 0 0];
pub.BFylab(ylh,T.F)

yyaxis right
ax=gca;
ax.YColor=[0 0 0];
plot(energy,spec(:,2)/sum(spec(:,2)),'Displayname','x-ray spectrum','color',[0 0 0])
ylh=ylabel('normalized spectrum');
yticks([]);

xlh=xlabel('Energy [keV]');
ylh.Color=[0 0 0];
lh=legend();
fname=sprintf('%sattenuation',T.Fig.saveto);

pub.BFfigure(fh,T.F)
pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)

savefig(fh,fname)
print(fh, '-dpdf',fname);
print(fh, '-dpng',fname);
end

function plotspec(energy,spec,cufilts,T)
globcol=get(groot,'DefaultAxesColorOrder'); %get colors
fh=figure(3250);clf
ax=gca;
hold on
grid on
col=copper(length(cufilts));

for i =1:length(cufilts) % attenuations
    p(i)=plot(energy,spec(:,i),'Displayname',...
        sprintf('%2.1fmm Cu',cufilts(i)),'color',col(i,:));
    pub.BFLine(p(i),T.F)
end
ylh=ylabel('normalized flux');

pub.BFylab(ylh,T.F)
yticklabels({})
xlh=xlabel('Energy [keV]');
lh=legend(p(1:5));
fname=sprintf('%spectrum1',T.Fig.saveto);

pub.BFfigure2(fh,T.F,1)
pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)

savefig(fh,fname)
print(fh, '-dpdf',fname);
print(fh, '-dpng',fname);


fh=figure(3251);clf
ax=gca;
hold on
grid on
col=copper(length(cufilts));
for i =1:length(cufilts) % attenuations
    p(i)=plot(energy,spec(:,i)/sum(spec(:,i)),'Displayname',...
        sprintf('%2.1fmm Cu',cufilts(i)),'color',col(i,:));
    pub.BFLine(p(i),T.F)
end
ylh=ylabel('normalized flux');

pub.BFylab(ylh,T.F)
yticklabels({})
xlh=xlabel('Energy [keV]');
lh=legend(p([1,6:10]),'location','northwest');
fname=sprintf('%spectrum2',T.Fig.saveto);

pub.BFfigure2(fh,T.F,1)
pub.BFlegend(lh,T.F)
pub.BFaxis(ax,T.F)
pub.BFylab(ylh,T.F)
pub.BFxlab(xlh,T.F)

savefig(fh,fname)
print(fh, '-dpdf',fname);
print(fh, '-dpng',fname);
end

function plotPhans(phan,phanCoord,T)
    % viszualizes the phantom
    fh=figure(2346);clf
    
    s=size(phan,3);
    im=zeros(size(phan,1),size(phan,2));
    ang=linspace(-pi,pi,s+1);
    dang=2*pi/s;
    
    for i=1:s
        tempphan=phan(:,:,i);
        mask=atan2(phanCoord(:,:,2),phanCoord(:,:,1))>=ang(i);
        mask(atan2(phanCoord(:,:,2),phanCoord(:,:,1))>ang(i+1))=0;
        im(mask)=tempphan(mask)+0.1*i;
    end
    
    recsize=size(phan,1);
    im=imresize(im,[recsize,recsize]);
    imshow(im',[]);set(gca,'YDir','normal')
    ax=gca;
    
    recsize=500;
    ax.Units='pixel';
    ax.Position=[0 0 recsize recsize];
    fh.Units='pixel';
    fh.Position=[50 50 recsize recsize];
    fh.Units='centimeters';
    fh.PaperPosition=fh.Position;
    fh.PaperSize=fh.Position(3:4);
    set(fh, 'PaperUnits', 'centimeters')
    set(fh, 'PaperPosition', [0 0 fh.Position(3:4)])
    
    fname=sprintf('PhantomS');
    print(sprintf('%s%s.pdf',T.F.saveto,fname),'-dpdf')
    print(sprintf('%s%s.png',T.F.saveto,fname),'-dpng')
    
    savefig(fh,fname)

end

function [pixflux,sourceBinSize]=getFluxes(geo,spectrum,energy,filterCase)
    % calculated values
exposurePerProjection = geo.beamCurrent*geo.totalMeasuringTime/geo.nProjections; % in [mAs]

pixelSizeRatio = geo.ph_n_coarse/geo.ph_n_fine; % ratio of coarse to fine
if floor(pixelSizeRatio)~=(pixelSizeRatio),disp('Error! PixelSizeRatio is non-integer.');return;end
% source energy steps of 1 keV from 10 to 110

% source output photons/keV/cm^2/mAs @ 1 meter
% Y_sourceSpectrum = spectrum(:,filterCase);
sourceBinSize = energy(2)-energy(1); % energy bin of source spectrumd data, keV
% the scripts assume 1 keV bins, if not then everything needs to be checked:
if sourceBinSize~=(1),disp('Error! Source bin size not 1 keV as expected. Check everything.');return;end
%if geo.ph_n_fine~=(0.01),disp('Error! pixelSizeFine not = 0.01 mm as expected. Check everything.');return;end

% ex-Y_sourceSpectrum in N[keV cm^2 mAs]^-1 @ 1 meter * cm *cm *mAs = [N/keV/pixel] @ 1m 
pixflux = spectrum(:,filterCase) * geo.det_vert/10 * geo.det_p/10 * exposurePerProjection ; % convert to #/keV/pixel/projection @ 1 meter
end

function [pixFluxSet,sourceBinSize]=getFluxSet(geo,spectrum,energy,cufilts)
        pixFluxSet=zeros(length(energy),length(cufilts));
    for i=1:length(cufilts)
        [pixFluxSet(:,i),sourceBinSize]=val.getFluxes(geo,spectrum,energy,i);
    end
end

function plotBinningSpec(sp,energy,spec,specbin,attbin,energyBinStart,energyBinSize,cufilts,T)
    % plots the binning stuff
    fh=figure(8462);
    %set(fh,'units','normalized','outerposition',[0 0 1 1]);
    
    
    ax=subplot(sp(1),sp(2),sp(3));cla;
    hold on
    %set(ax, 'YScale', 'log')
    plot(energy,spec(:,1),'color',[0 0 0]);
    col=copper(size(specbin,2));
    for i=1:size(specbin,2)
        p(i)=bar(energyBinStart+floor(energyBinSize/2),specbin(:,i)/energyBinSize,...
            'FaceColor',col(i,:),'Displayname',...
        sprintf('%2.1fmm Cu',cufilts(i)));
    end
    ylim([1,max(spec(:))])
    grid on
    xlabel('energy [keV]')
    ylabel('spectrum')
    legend(p)
    
end

function plotBinningAtt(sp,energy,att,attbin,energyBinStart,energyBinSize,T)
    % plots the binning stuff
    fh=figure(8463);
    %set(fh,'units','normalized','outerposition',[0 0 1 1]);
    globcol=get(groot,'DefaultAxesColorOrder'); %get colors
    colind=[2,1,3];
    
    name={'CHCl3','H2O','Alu'};
    ax=subplot(sp(1),sp(2),sp(3));cla;
    title(sprintf('Binned XS at Binsize %dkeV and 0.2mmCu',energyBinSize))
    hold on
    set(ax, 'YScale', 'log')
    
    
    for i=[3,1,2]
        plot(energy,att(:,i),'color',globcol(colind(i),:));
        p(i)=bar(energyBinStart+floor(energyBinSize/2),attbin(:,i),...
            'Displayname',name{i},'FaceColor',globcol(colind(i),:),'FaceAlpha',1);
    end
    ylim([0.1,max(att(:))])
    xlabel('energy [keV]')
    ylabel('attenuation [1/cm]')
    grid on
    legend(p)
    
end

function attphan=ImbuePhantom(phan,attbin,geo)
    % puts the attenuation coefficient into the phantom,
    % returns a (x,y,lft,enbin)-dim array
    
    apsize =size(phan,1)/geo.phanResizeRatio;% size of attphan
    if mod(apsize,1)~=0
        error('reduced phatom pixel size is not integer')
    end
    nen = size(attbin,1); % number of energy bins
    nlft=size(phan,3); % number of lfts

    % the new attenuation coefficient phantom:
    attphan=zeros(apsize,apsize,nlft,nen);
    
    for lft=1:nlft
        for en=1:nen
            tempphan=val.attToPhan(phan(:,:,lft),attbin(en,:));
            attphan(:,:,lft,en)=imresize(tempphan,...
                1/geo.phanResizeRatio,'box','Method',...
                'bilinear','Antialiasing',false);% downscale image
               
        end
    end
end

function im=attToPhan(im,atts)
    % takes a single image and puts in attenuations
    im(im==0)=atts(3);% alu
    im(im==1)=0;% air
    im(im==2)=atts(2);% h2o
    im(im==3)=atts(1)*4.475/1411;% chloro vap
    im(im==4)=atts(1);% chloro
end

function rec2=fan_new(sinogram,ang,src,det,n_r,p_r,p_d)
% sinogram= sinogram, size(sino) =[n_ang,n_pix]
% ang = angle list
% src = source - origin distance
% det = detector-origin distance
% n_r = width of recon volume in pixel
% p_r = recon volume pixel pitch in mm
%

d_r=p_r*n_r;
n_d=size(sinogram,2);
n_ang=length(ang);

vol_geom2 = astra_create_vol_geom(n_r, n_r,-d_r/2, d_r/2, -d_r/2, d_r/2);
proj_geom2 = astra_create_proj_geom('fanflat',p_d, n_d, ang, src, det);

%proj_id2 = astra_create_projector('strip_fanflat', proj_geom, vol_geom);
sinogram_id2 = astra_mex_data2d('create', '-sino', proj_geom2, sinogram);
rec_id2 = astra_mex_data2d('create', '-vol', vol_geom2);

% create configuration
cfg2 = astra_struct('FBP_CUDA');
%cfg2.ProjectorId = proj_id2;
cfg2.ReconstructionDataId = rec_id2;
cfg2.ProjectionDataId = sinogram_id2;
%cfg.FilterType = 'Ram-Lak';
alg_id2 = astra_mex_algorithm('create', cfg2);
astra_mex_algorithm('run', alg_id2);

% get the reconstruction
%rec = astra_mex_data2d('get', rec_id2);
rec2 = astra_mex_data2d('get', rec_id2)*((p_r/p_d)^2);

rec2=rec2/p_d*n_ang/(pi/2); % <= Question 1) & 2)

rec2=rec2/(src/(det+src)); % magnification correction of the astra fan beam

astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');

end

function rec2=par_new(sinogram,ang,n_r,p_r,p_d)
% sinogram= sinogram, size(sino) =[n_ang,n_pix]
% ang = angle list
% src = source - origin distance
% det = detector-origin distance
% n_r = width of recon volume in pixel
% p_r = recon volume pixel pitch in mm
%

d_r=p_r*n_r;
n_d=size(sinogram,2);

vol_geom2 = astra_create_vol_geom(n_r, n_r,-d_r/2, d_r/2, -d_r/2, d_r/2);
proj_geom2 = astra_create_proj_geom('parallel',p_d, n_d, ang);

%proj_id2 = astra_create_projector('strip_fanflat', proj_geom, vol_geom);
sinogram_id2 = astra_mex_data2d('create', '-sino', proj_geom2, sinogram);
rec_id2 = astra_mex_data2d('create', '-vol', vol_geom2);

% create configuration
cfg2 = astra_struct('FBP_CUDA');
%cfg2.ProjectorId = proj_id2;
cfg2.ReconstructionDataId = rec_id2;
cfg2.ProjectionDataId = sinogram_id2;
%cfg.FilterType = 'Ram-Lak';
alg_id2 = astra_mex_algorithm('create', cfg2);
astra_mex_algorithm('run', alg_id2);

% get the reconstruction
%rec = astra_mex_data2d('get', rec_id2);
rec2 = astra_mex_data2d('get', rec_id2)*((p_r/p_d)^2);
rec2 = rec2*p_d;  % fix 1 parallel

astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');

end

function sinogram=parforward_new(phan,n_ph,d_ph,p_d,n_d,ang)

% phan = phantom image
% n_ph = width of phantom in pixel
% d_ph = width of phantom in mm
% p_d = detector pixel pitch
% n_d = detector pixel amount
% and = angle list in radian

proj_geom = astra_create_proj_geom('parallel', p_d ,n_d, ang);
vol_geom = astra_create_vol_geom(n_ph,n_ph,-d_ph/2,d_ph/2,-d_ph/2,d_ph/2);

% store volume
volume_id = astra_mex_data2d('create', '-vol', vol_geom, phan);

% create forward projection
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.VolumeDataId = volume_id;
fp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fp_id);
sinogram = astra_mex_data2d('get', sinogram_id);

sinogram=sinogram/p_d;  % fix 1 parallel

% garbage disposal
astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');
end

function sinogram=fanforward_new(phan,d_ph,p_d,n_d,ang,src,det)

% phan = phantom image
% n_ph = width of phantom in pixel
% d_ph = width of phantom in mm
% p_d = detector pixel pitch
% n_d = detector pixel amount
% and = angle list in radian

n_ph=size(phan,1);%required square phan
proj_geom = astra_create_proj_geom('fanflat',p_d, n_d, ang, src, det);
vol_geom = astra_create_vol_geom(n_ph,n_ph,-d_ph/2,d_ph/2,-d_ph/2,d_ph/2);

% store volume
volume_id = astra_mex_data2d('create', '-vol', vol_geom, phan);

% create forward projection
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
cfg = astra_struct('FP_CUDA');
cfg.ProjectionDataId = sinogram_id;
cfg.VolumeDataId = volume_id;
fp_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', fp_id);
sinogram = astra_mex_data2d('get', sinogram_id);

% garbage disposal
astra_mex_data2d('clear');
astra_mex_data3d('clear');
astra_mex_algorithm('clear');
astra_mex_projector('clear');
astra_mex_projector3d('clear');
astra_mex_matrix('clear');
end

function HistCompare(fid,sinosim,sinoreal,filt,neBins,T)
   % compares two sinograms and their histograms
   fh=figure(fid);clf
   crange=[0 1];
   subplot(1,4,1)
   imshow(sinosim,crange);
   title('sim')
   
   subplot(1,4,2)
   imshow(sinoreal,crange)
   title('real')
   
   subplot(1,4,3:4)
   hold on
   histogram(sinosim,'displayname','sim')
   histogram(sinoreal,'displayname','real')
   grid on
   legend()
   title(sprintf('%d E-bins, %2.1fmmCu',neBins,filt))
   
   fname=sprintf('%sattenuation',T.Fig.saveto);
   fname=sprintf('histcompare %d',fid);
   savefig(fh,fname)
print(fh, '-dpdf',fname);
print(fh, '-dpng',fname);
end

function ShowRecons(recon,lfts)
    
    crange=[min(recon(:)),max(recon(:))];

    for i=2:size(recon,3)
        figure(12456+i);clf;
        imshow(recon(:,:,i)-recon(:,:,1),[crange])
    end
    
end

function ReconHistCompare(phan,recsim,recreal,filt)
    fh=figure(4267);clf;
    subplot(1,3,1);
    
    im= cat(1,recsim,recreal,phan);
    imshow(im,[]);
    subplot(1,3,2:3)
    histogram(recsim,100);
    hold on
    histogram(recreal,100);
    %histogram(phan,100,'FaceAlpha',0.5)
    legend({'sim','real','phantom'})
    ylim([0 1e4])
    grid on
    title(sprintf('%2.1f mmCu',filt))
    

    
    
end

end
end

























