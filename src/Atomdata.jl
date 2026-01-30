"""
    module Atomdata

Holds the data for isolated atoms.
"""
module Atomdata
using ..AtomicMod:makeatom
using ..ThreeBodyTB:SRCDIR
#using ..Atomic:AtomicMod



#######old













































"""
    atoms::Dict()

Periodic table information.
"""
atoms = Dict()

#      name                     name  Z  row col mass nval nsemi  orbitals total_energy    e_orb1, e_orb2
#                                                (amu)                      (Ryd)             (eV)
#
#atoms["T"] = makeatom("T",    0, 0, 0, 0.0,  1.0,   0,   [:s], 0.0, [0.0], 0.0);

#atoms["Is"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p], -9.15243175, [-10.6594, -3.9477],0.007075865418, 0.5300262399999837);

atoms["Hx" ] = makeatom("Hx" , 100, 1, 1, 1.0079,   1.0,   0,   [:s ],     -0.90531367, [-6.3233], 0.0, 0.9194820210526328);
atoms["X"] = makeatom("X",    101, 0, 0, 0.0,  1.0,   0,   [:s], 0.0, [0.0], 0.0, 0.5);
atoms["Xa"] = makeatom("Xa",    102, 0, 0, 0.0,  1.0,   0,   [:s], 3.0, [3.0* 13.605693122], 0.0, 0.5);



atoms["H" ] = makeatom("H" , 1, 1, 1, 1.0079,   1.0,   0,   [:s ],     -0.90531367, [-6.3176], 0.0, 0.8451632950940942 );
atoms["Li"] = makeatom("Li", 3, 2, 1, 6.941,    1.0,   2,   [:s, :p], -14.23107803, [-2.7537, -0.99],0.0, 0.3902587823285621 );
atoms["Be"] = makeatom("Be", 4, 2, 2, 9.012,    2.0,   2,   [:s, :p], -27.80694399, [-5.5055, -1.9313],0.0, 0.5751048209185653 );
atoms["B" ] = makeatom("B" , 5, 2, 3, 10.81,    3.0,   0,   [:s, :p], -5.84632720, [-9.3213, -3.5217], 0.0, 0.6267943641604508 );
atoms["C" ] = makeatom("C" , 6, 2, 4, 12.011,   4.0,   0,   [:s, :p], -10.78129974, [-13.5882, -5.1806],0.0, 0.7978153170555154 );
atoms["N" ] = makeatom("N" , 7, 2, 5, 14.007,   5.0,   0,   [:s, :p], -19.47573456, [-18.3876, -6.9917],0.0, 0.8874159249179556 );
atoms["O" ] = makeatom("O" , 8, 2, 6, 15.999,   6.0,   0,   [:s, :p], -31.71324769, [-23.7250, -8.9381],0.0, 1.015623980577927 );
atoms["F" ] = makeatom("F" , 9, 2, 7, 18.998,   7.0,   0,   [:s, :p], -48.22248233, [-29.6811, -11.0125],0.0, 1.1398170857521994 );
atoms["Na"] = makeatom("Na", 11, 3, 1, 22.98,  1.0,   8,   [:s, :p], -95.06749946, [-2.6455, -0.6556], 0.0, 0.38560892349289677 );
atoms["Mg"] = makeatom("Mg", 12, 3, 2, 24.31,  2.0,   8,   [:s, :p], -124.99127898, [-4.6187, -1.2416], 0.0, 0.5106562587165097 );
atoms["Al"] = makeatom("Al", 13, 3, 3, 26.98,  3.0,   0,   [:s, :p],  -6.37406439, [-7.6515, -2.6256], 0.0, 0.47560947595799796 );
atoms["Si"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p], -9.15243181, [-10.7049, -3.9932],0.0, 0.5597222300055447 );
#atoms["Sid"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p, :d], -9.15243181, [-10.7049, -3.9932, 0.61],0.0, 0.5597222300055447 );
atoms["P" ] = makeatom("P" , 15, 3, 5, 30.97,  5.0,   0,   [:s, :p], -15.07327848, [-13.8557, -5.4181], 0.0, 0.6422064900805308 );
atoms["S" ] = makeatom("S" , 16, 3, 6, 32.06,  6.0,   0,   [:s, :p], -23.70562681, [-17.1360, -6.9201], 0.0, 0.7079117127047072 );
atoms["Cl"] = makeatom("Cl", 17, 3, 7, 35.45,  7.0,   0,   [:s, :p], -33.06905525, [-20.5627, -8.5062], 0.0, 0.7476423278619385 );
atoms["K" ] = makeatom("K" , 19, 4, 1, 39.098,  1.0,   8,   [:s, :d, :p], -56.84482202 , [-2.1989, -0.1312, -0.6747], 0.0, 0.37734366682102694 );
atoms["Ca"] = makeatom("Ca", 20, 4, 2, 40.078,  2.0,   8,   [:s, :d, :p], -74.68919704, [-3.6547, -1.8454, -1.2718], 0.0, 0.4749594731162585 );
atoms["Sc"] = makeatom("Sc", 21, 4, 2.05, 44.956,  3.0,   8,   [:s, :d, :p], -93.89389249,  [-4.0850, -3.1428, -1.3456],  0.0, 0.7781746286900919 );
atoms["Ti"] = makeatom("Ti", 22, 4, 2.15, 47.867,  4.0,   8,   [:s, :d, :p], -118.63461912, [-4.3470, -4.0656, -1.3310],  0.0, 0.5344202735274272 );
atoms["V" ] = makeatom("V" , 23, 4, 2.25, 50.941,  5.0,    8,   [:s, :d, :p], -144.34826807, [-4.3727, -4.2223, -1.2153], 0.0, 0.5051179562238686 );
atoms["Cr"] = makeatom("Cr", 24, 4, 2.35, 51.996,  6.0,    8,   [:s, :d, :p], -174.69181247, [-4.3746, -4.2928, -1.0764], 0.0, 0.5023896126040788 );
atoms["Mn"] = makeatom("Mn", 25, 4, 2.45, 54.938,  7.0,    8,   [:s, :d, :p], -210.40629067, [-4.3786, -4.3505, -0.9746], 0.0, 0.5085605056527054 );
atoms["Fe"] = makeatom("Fe", 26, 4, 2.55, 55.845,  8.0,    8,   [:s, :d, :p], -249.56808439, [-4.3911, -4.4108, -0.8715], 0.0, 0.5138399863052884 );
atoms["Co"] = makeatom("Co", 27, 4, 2.65, 58.993,  9.0,    8,   [:s, :d, :p], -297.72249071, [-4.2335, -4.2995, -0.5167], 0.0, 0.5205632537537023 );
atoms["Ni"] = makeatom("Ni", 28, 4, 2.75, 58.963,  10.0,   8,   [:s, :d, :p], -342.65159683, [-4.2480, -4.3731, -0.3946], 0.0, 0.530332541897714 );
atoms["Cu"] = makeatom("Cu", 29, 4, 2.85, 63.546,  11.0,   8,   [:s, :d, :p], -403.19763559, [-4.5506, -4.9911, -0.5855], 0.0, 0.5619138137847363 );
atoms["Zn"] = makeatom("Zn", 30, 4, 2.95, 65.38,   12.0,   8,   [:s, :d, :p], -461.12493484, [-5.9305, -10.1008, -1.0049], 0.0, 0.6233810643339528 );
atoms["Ga"] = makeatom("Ga", 31, 4, 3.0,  69.72,   13.0,   6,   [:s,:d, :p], -414.75774465, [-8.8812,-19.0647, -2.4866], 0.0, 0.5103937462614931 );
atoms["Ge"] = makeatom("Ge", 32, 4, 4.0,  72.63,   4.0,    10,   [:s, :p], -212.84384318, [-11.6643, -3.8180], 0.0, 0.5643618220580027 );
#atoms["Ged"] = makeatom("Ge", 32, 4, 4.0,  72.63,   4.0,    10,   [:s, :p, :d], -212.84384318, [-11.6643, -3.8180, 0.56], 0.0, 0.5643618220580027 );
atoms["As"] = makeatom("As", 33, 4, 5.0,  74.92,   5.0,    0,   [:s, :p], -39.62501000, [-14.4349, -5.1279], 0.0, 0.6154159430456763 );
atoms["Se"] = makeatom("Se", 34, 4, 6.0,  78.97,   6.0,    0,   [:s, :p],  -43.03651432, [-17.2230, -6.4409], 0.0, 0.6730297191236639 );
atoms["Br"] = makeatom("Br", 35, 4, 7.0,  79.90,   7.0,    0,   [:s, :p], -40.38769338, [-20.0640, -7.7862], 0.0, 0.7295557399434305 );
atoms["Rb"] = makeatom("Rb", 37, 5, 1, 85.468,  1.0,   8,   [:s, :d, :p], -53.05281238, [-2.0994, -0.3, -0.5921], 0.0, 0.38651513892406236 );
atoms["Sr"] = makeatom("Sr", 38, 5, 2, 87.62,  2.0,   8,   [:s, :d, :p], -69.92604259, [-3.4017, -1.1543, -1.1318], 0.0, 0.47280011784661014 );
atoms["Y" ] = makeatom("Y" , 39, 5, 2.05, 88.906,  3.0,   8,    [:s, :d, :p], -92.30423365,  [-3.9864, -2.3907, -1.2867], 0.0, 0.6384727135899312 );
atoms["Zr"] = makeatom("Zr", 40, 5, 2.15, 91.224,  4.0,   8,    [:s, :d, :p], -98.57421435,  [-4.3347, -3.4817, -1.3201], 0.0, 0.6903649183604765 );
atoms["Nb"] = makeatom("Nb", 41, 5, 2.25, 92.906,  5.0,    8,   [:s, :d, :p], -117.19066274 ,  [-4.5053, -4.3124, -1.3159],0.0, 0.5348560719045159 );
atoms["Mo"] = makeatom("Mo", 42, 5, 2.35, 95.95,   6.0,    8,   [:s, :d, :p], -138.51328928,  [-4.2469, -4.2022, -1.0704], 0.0, 0.5098857597525097 );
atoms["Tc"] = makeatom("Tc", 43, 5, 2.45, 98.0,    7.0,    8,   [:s, :d, :p], -173.39188188,  [-3.9325, -3.9874, -0.7985], 0.0, 0.4983977286525886 );
atoms["Ru"] = makeatom("Ru", 44, 5, 2.55, 101.07,  8.0,    8,   [:s, :d, :p], -193.54935366,  [-3.5929, -3.7493, -0.4736], 0.0, 0.511729036259258 );
atoms["Rh"] = makeatom("Rh", 45, 5, 2.65, 102.91,  9.0,    6,   [:s, :d, :p], -170.64161619,  [-3.3079, -3.6154, -0.2568], 0.0, 0.6907390586058985 );
atoms["Pd"] = makeatom("Pd", 46, 5, 2.75, 106.42,  10.0,   6,   [:s, :d, :p], -198.31288985,  [-3.2109, -3.9641, -0.1564], 0.0, 0.7801811214471324 );
atoms["Ag"] = makeatom("Ag", 47, 5, 2.85, 107.87,  11.0,   8,   [:s, :d, :p], -295.31434496,  [-4.3268, -7.3569, -0.5745], 0.0, 0.5706935927525575 );
atoms["Cd"] = makeatom("Cd", 48, 5, 2.95, 112.41,  12.0,   0,   [:s, :d, :p], -114.39800378,  [-5.6202, -11.5954, -1.0602], 0.0, 0.6005115028052722 );
atoms["In"] = makeatom("In", 49, 5, 3.0,  114.82,   13.0,   0,   [:s,:d, :p], -132.97861648 , [-8.1551, -18.4828, -2.4193], 0.0, 0.49245976583127 );
atoms["Sn"] = makeatom("Sn", 50, 5, 4.0,  118.71,   4.0,   10,   [:s, :p], -158.61463088, [-10.4877, -3.5995], 0.0, 0.5456114086239429 );
atoms["Sb"] = makeatom("Sb", 51, 5, 5.0,  121.76,   5.0,   10,   [:s, :p], -184.28336729, [-12.7656, -4.7332], 0.0, 0.5973377632845377 );
atoms["Te"] = makeatom("Te", 52, 5, 6.0,  127.60,   6.0,    0,   [:s, :p], -26.19141773, [-15.0608, -5.8767], 0.0, 0.6243066065034633 );
atoms["I" ] = makeatom("I" , 53, 5, 7.0,  126.90,   7.0,    0,   [:s, :p], -65.01087258 , [-17.3546, -7.0038], 0.0, 0.6675305924093194 );
atoms["Cs"] = makeatom("Cs", 55, 6, 1, 132.9,  1.0,   8,   [:s, :d, :p ], -62.83174426, [-1.9276, -0.3, -0.5431 ], 0.0, 0.4013353474699838 );
atoms["Ba"] = makeatom("Ba", 56, 6, 2, 139.3,  2.0,   8,   [:s, :d, :p], -70.01638694, [-3.0703, -1.8523, -1.0457], 0.0, 0.48617114830668345 );
atoms["La"] = makeatom("La", 57, 6, 2.05, 138.9,  3.0,   8,    [:s, :d], -101.82206370,  [-3.2213,-2.43], 0.0, 0.5816142601437702 );
atoms["Hf"] = makeatom("Hf", 72, 6, 2.15, 178.5,  4.0,   8,    [:s, :d, :p], -157.97552412,  [-5.0187, -2.6052, -1.2994], 0.0, 0.6517565809524214 );
atoms["Ta"] = makeatom("Ta", 73, 6, 2.25, 181.0,  5.0,    8,   [:s, :d, :p], -141.00656328,  [-5.3481, -3.4993, -1.2418], 0.0, 0.6912754847554329 );
atoms["W" ] = makeatom("W" , 74, 6, 2.35, 183.8,  6.0,    8,   [:s, :d, :p], -158.16868802,  [-5.6245, -4.3914, -1.2092], 0.0, 0.7326044378416148 );
atoms["Re"] = makeatom("Re", 75, 6, 2.45, 186.2,  7.0,    8,   [:s, :d, :p], -189.48871022,  [-5.8618, -5.2763, -1.1642], 0.0, 0.7652760081244261 );
atoms["Os"] = makeatom("Os", 76, 6, 2.55, 190.2,  8.0,    8,   [:s, :d, :p], -197.71975092,  [-5.9152, -5.8238, -1.0547], 0.0, 0.6076951539696345 );
atoms["Ir"] = makeatom("Ir", 77, 6, 2.65, 192.2,  9.0,    6,   [:s, :d, :p], -180.96527345,  [-5.5876, -5.6416, -0.8041], 0.0, 0.5861813291765662 );
atoms["Pt"] = makeatom("Pt", 78, 6, 2.75, 195.1,  10.0,   6,   [:s, :d, :p], -210.06169163,  [-5.3377, -5.5187, -0.5975], 0.0, 0.5839433031650294 );
atoms["Au"] = makeatom("Au", 79, 6, 2.85, 197.0,  11.0,   0,   [:s, :d, :p], -114.46631173,  [-5.7155, -6.8176, -0.6483], 0.0, 0.611432813778228 );
atoms["Hg"] = makeatom("Hg", 80, 6, 2.95, 200.6,  12.0,   0,   [:s, :d, :p], -107.48183748,  [-6.7109, -9.7253, -0.8339], 0.0, 0.6495200693187368 );
atoms["Tl"] = makeatom("Tl", 81, 6, 3.0,  204.4,   13.0,   0,   [:s, :d, :p], -144.53037033, [ -9.4990,-15.3839,  -2.2651], 0.0, 0.5007833337462329 );
atoms["Pb"] = makeatom("Pb", 82, 6, 4.0,  207.2,   4.0,   10,   [:s, :p], -163.61141039, [-11.9396, -3.3624], 0.0, 0.5488505326347655 );
atoms["Bi"] = makeatom("Bi", 83, 6, 5.0,  209.0,   5.0,   10,   [:s, :p], -184.68674216, [-14.4080, -4.4700], 0.0, 0.597524848429184 );




atoms["B_d" ] = makeatom("B" , 5, 2, 3, 10.81,    3.0,   0,   [:s, :p, :d], -5.84632720, [-9.3213, -3.5217, 0.6], 0.0, 0.6267943641604508 );
atoms["C_d" ] = makeatom("C" , 6, 2, 4, 12.011,   4.0,   0,   [:s, :p, :d], -10.78129974, [-13.5882, -5.1806, 0.6],0.0, 0.7978153170555154 );
atoms["N_d" ] = makeatom("N" , 7, 2, 5, 14.007,   5.0,   0,   [:s, :p, :d], -19.47573456, [-18.3876, -6.9917, 0.6],0.0, 0.8874159249179556 );
atoms["O_d" ] = makeatom("O" , 8, 2, 6, 15.999,   6.0,   0,   [:s, :p, :d], -31.71324769, [-23.7250, -8.9381, 0.6],0.0, 1.015623980577927 );
atoms["F_d" ] = makeatom("F" , 9, 2, 7, 18.998,   7.0,   0,   [:s, :p, :d], -48.22248233, [-29.6811, -11.0125, 0.6],0.0, 1.1398170857521994 );
atoms["Al_d"] = makeatom("Al", 13, 3, 3, 26.98,  3.0,   0,   [:s, :p, :d],  -6.37406439, [-7.6515, -2.6256, 0.6], 0.0, 0.47560947595799796 );
atoms["Si_d"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p, :d], -9.15243181, [-10.7049, -3.9932, 0.6],0.0, 0.5597222300055447 );
atoms["P_d" ] = makeatom("P" , 15, 3, 5, 30.97,  5.0,   0,   [:s, :p, :d], -15.07327848, [-13.8557, -5.4181, 0.6], 0.0, 0.6422064900805308 );
atoms["S_d" ] = makeatom("S" , 16, 3, 6, 32.06,  6.0,   0,   [:s, :p, :d], -23.70562681, [-17.1360, -6.9201, 0.6], 0.0, 0.7079117127047072 );
atoms["Cl_d"] = makeatom("Cl", 17, 3, 7, 35.45,  7.0,   0,   [:s, :p, :d], -33.06905525, [-20.5627, -8.5062, 0.6], 0.0, 0.7476423278619385 );
atoms["Ge_d"] = makeatom("Ge", 32, 4, 4.0,  72.63,   4.0,    10,   [:s, :p, :d], -212.84384318, [-11.6643, -3.8180, 0.55], 0.0, 0.5643618220580027 );
atoms["As_d"] = makeatom("As", 33, 4, 5.0,  74.92,   5.0,    0,   [:s, :p, :d], -39.62501000, [-14.4349, -5.1279,  0.55], 0.0, 0.6154159430456763 );
atoms["Se_d"] = makeatom("Se", 34, 4, 6.0,  78.97,   6.0,    0,   [:s, :p, :d],  -43.03651432, [-17.2230, -6.4409, 0.55], 0.0, 0.6730297191236639 );
atoms["Br_d"] = makeatom("Br", 35, 4, 7.0,  79.90,   7.0,    0,   [:s, :p, :d], -40.38769338, [-20.0640, -7.7862, 0.55], 0.0, 0.7295557399434305 );
atoms["Sn_d"] = makeatom("Sn", 50, 5, 4.0,  118.71,   4.0,   10,   [:s, :p, :d], -158.61463088, [-10.4877, -3.5995, 0.55], 0.0, 0.5456114086239429 );
atoms["Sb_d"] = makeatom("Sb", 51, 5, 5.0,  121.76,   5.0,   10,   [:s, :p, :d], -184.28336729, [-12.7656, -4.7332, 0.55], 0.0, 0.5973377632845377 );
atoms["Te_d"] = makeatom("Te", 52, 5, 6.0,  127.60,   6.0,    0,   [:s, :p, :d], -26.19141773, [-15.0608, -5.8767, 0.55], 0.0, 0.6243066065034633 );
atoms["I_d" ] = makeatom("I" , 53, 5, 7.0,  126.90,   7.0,    0,   [:s, :p, :d], -65.01087258 , [-17.3546, -7.0038, 0.55], 0.0, 0.6675305924093194 );
atoms["Pb_d"] = makeatom("Pb", 82, 6, 4.0,  207.2,   4.0,   10,   [:s, :p, :d], -163.61141039, [-11.9396, -3.3624, 0.55], 0.0, 0.5488505326347655 );
atoms["Bi_d"] = makeatom("Bi", 83, 6, 5.0,  209.0,   5.0,   10,   [:s, :p, :d], -184.68674216, [-14.4080, -4.4700, 0.55], 0.0, 0.597524848429184 );

zatoms = Dict()

for key in keys(atoms)
    println("key $key typeof ", typeof(key))
    if occursin("_d", key)
        continue
    end
    zatoms[atoms[key].Z] = key
end
for key in keys(atoms)
    atoms[Symbol(key)] = atoms[key]
end

#U = IE - EA  (ionization energy - electron affinity

#stop printing
nothing

#approximate metallic radius / single covalent radius 
#for use in prerelaxing structures in pm
"""
    atom_radius::Dict()

Atomic radius (metallic), in pm
"""
atom_radius = Dict()

atom_radius[ "Hx" ] =  100

atom_radius[ "X" ] =  200
atom_radius[ "Xa" ] =  200
atom_radius[ "T" ] =  200

atom_radius[ "H" ] =  60
atom_radius[ "He" ] =  60
atom_radius[ "Li" ] =  152
atom_radius[ "Be" ] =  112
atom_radius[ "B" ] =   82
atom_radius[ "C" ] =   77
atom_radius[ "N" ] =   75
atom_radius[ "O" ] =   73
atom_radius[ "F" ] =   71
atom_radius[ "Ne" ] =  80
atom_radius[ "Na" ] =  186
atom_radius[ "Mg" ] =  160
atom_radius[ "Al" ] =  143
atom_radius[ "Si" ] =  111
atom_radius[ "P" ] =   106
atom_radius[ "S" ] =   102
atom_radius[ "Cl" ] =  99
atom_radius[ "Ar" ] =  100
atom_radius[ "K" ] =   227
atom_radius[ "Ca" ] =  197
atom_radius[ "Sc" ] =  162
atom_radius[ "Ti" ] =  147
atom_radius[ "V" ] =   134
atom_radius[ "Cr" ] =  128
atom_radius[ "Mn" ] =  127
atom_radius[ "Fe" ] =  126
atom_radius[ "Co" ] =  125
atom_radius[ "Ni" ] =  124
atom_radius[ "Cu" ] =  128
atom_radius[ "Zn" ] =  134
atom_radius[ "Ga" ] =  135
atom_radius[ "Ge" ] =  122
atom_radius[ "Ged" ] =  122
atom_radius[ "As" ] =  119
atom_radius[ "Se" ] =  116
atom_radius[ "Br" ] =  114
atom_radius[ "Kr" ] =  100
atom_radius[ "Rb" ] =  248
atom_radius[ "Sr" ] =  215
atom_radius[ "Y" ] =  180
atom_radius[ "Zr" ] =  160
atom_radius[ "Nb" ] =  146
atom_radius[ "Mo" ] =  139
atom_radius[ "Tc" ] =  136
atom_radius[ "Ru" ] =  134
atom_radius[ "Rh" ] =  134
atom_radius[ "Pd" ] =  137
atom_radius[ "Ag" ] =  144
atom_radius[ "Cd" ] =  151
atom_radius[ "In" ] =  167
atom_radius[ "Sn" ] =  141
atom_radius[ "Sb" ] =  138
atom_radius[ "Te" ] =  135
atom_radius[ "I" ] =   133
atom_radius[ "Xe" ] =  130
atom_radius[ "Cs" ] =  265
atom_radius[ "Ba" ] =  222
atom_radius[ "La" ] =  187
atom_radius[ "Ce" ] =  181
atom_radius[ "Pr" ] =  181
atom_radius[ "Nd" ] =  181
atom_radius[ "Pm" ] =  181
atom_radius[ "Sm" ] =  181
atom_radius[ "Eu" ] =  181
atom_radius[ "Gd" ] =  181
atom_radius[ "Tb" ] =  181
atom_radius[ "Dy" ] =  181
atom_radius[ "Ho" ] =  181
atom_radius[ "Er" ] =  181
atom_radius[ "Tm" ] =  181
atom_radius[ "Yb" ] =  181
atom_radius[ "Lu" ] =  173
atom_radius[ "Hf" ] =  159
atom_radius[ "Ta" ] =  146
atom_radius[ "W" ] =  139
atom_radius[ "Re" ] =  137
atom_radius[ "Os" ] =  135
atom_radius[ "Ir" ] =  135
atom_radius[ "Pt" ] =  135
atom_radius[ "Au" ] =  144
atom_radius[ "Hg" ] =  151
atom_radius[ "Tl" ] =  170
atom_radius[ "Pb" ] =  147
atom_radius[ "Bi" ] =  146
atom_radius[ "Po" ] =  140
atom_radius[ "At" ] =  140
atom_radius[ "Rn" ] =  130

for key in keys(atom_radius)
    atom_radius[key*"_d"] = atom_radius[key]
end

for key in keys(atom_radius)
    atom_radius[Symbol(key)] = atom_radius[key]
end



"""
    atom_prefered_oxidation::Dict()

Prefered oxidation states of atoms, in descending order of preference (approximate).
"""
atom_prefered_oxidation = Dict()


atom_prefered_oxidation[ "H" ] =  [1, -1, 0]
atom_prefered_oxidation[ "He" ] = [0]
atom_prefered_oxidation[ "Li" ] = [1, 0]
atom_prefered_oxidation[ "Be" ] = [2, 0]
atom_prefered_oxidation[ "B" ] =  [3, 0]
atom_prefered_oxidation[ "C" ] =  [4, -4, 0]
atom_prefered_oxidation[ "N" ] =  [-3, 0]
atom_prefered_oxidation[ "O" ] =  [-2, 0]
atom_prefered_oxidation[ "F" ] =  [-1, 0]
atom_prefered_oxidation[ "Ne" ] = [0]
atom_prefered_oxidation[ "Na" ] = [1, 0]
atom_prefered_oxidation[ "Mg" ] = [2, 0]
atom_prefered_oxidation[ "Al" ] = [3, 0]
atom_prefered_oxidation[ "Si" ] = [4, -4, 2, 0]
#atom_prefered_oxidation[ "Sid" ] = [4, -4, 2, 0]
atom_prefered_oxidation[ "P" ] =  [-3, 3, 5, 0]
atom_prefered_oxidation[ "S" ] =  [-2, 6, 0]
atom_prefered_oxidation[ "Cl" ] = [-1, 0]
atom_prefered_oxidation[ "Ar" ] = [0]
atom_prefered_oxidation[ "K" ] =  [1, 0]
atom_prefered_oxidation[ "Ca" ] = [2, 0]
atom_prefered_oxidation[ "Sc" ] = [3, 2, 0]
atom_prefered_oxidation[ "Ti" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "V" ] =  [5, 3, 0]
atom_prefered_oxidation[ "Cr" ] = [4, 3, 6, 2, 0]
atom_prefered_oxidation[ "Mn" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Fe" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Co" ] = [3, 4, 2, 0]
atom_prefered_oxidation[ "Ni" ] = [2, 4, 1, 0]
atom_prefered_oxidation[ "Cu" ] = [2, 4, 1, 0]
atom_prefered_oxidation[ "Zn" ] = [2, 0]
atom_prefered_oxidation[ "Ga" ] = [3, 0]
atom_prefered_oxidation[ "Ge" ] = [4, -4, 2, 0]
atom_prefered_oxidation[ "Ged" ] = [4, -4, 2, 0]
atom_prefered_oxidation[ "As" ] = [-3, 3, 5, 0]
atom_prefered_oxidation[ "Se" ] = [-2, 2, 4, 6, 0]
atom_prefered_oxidation[ "Br" ] = [-1, 0]
atom_prefered_oxidation[ "Kr" ] = [0]
atom_prefered_oxidation[ "Rb" ] = [1, 0]
atom_prefered_oxidation[ "Sr" ] = [2, 0]
atom_prefered_oxidation[ "Y" ] =  [3, 0]
atom_prefered_oxidation[ "Zr" ] = [4, 2, 0]
atom_prefered_oxidation[ "Nb" ] = [5, 0]
atom_prefered_oxidation[ "Mo" ] = [6, 4, 0]
atom_prefered_oxidation[ "Tc" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Ru" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Rh" ] = [4, 3, 2, 0]
atom_prefered_oxidation[ "Pd" ] = [4, 2, 0]
atom_prefered_oxidation[ "Ag" ] = [1, 2, 0]
atom_prefered_oxidation[ "Cd" ] = [2,0]
atom_prefered_oxidation[ "In" ] = [3,0]
atom_prefered_oxidation[ "Sn" ] = [4, 2, -4, 0]
atom_prefered_oxidation[ "Sb" ] = [-3, 3, 5, 0]
atom_prefered_oxidation[ "Te" ] = [-2, 6, 0]
atom_prefered_oxidation[ "I" ] =  [-1, 0]
atom_prefered_oxidation[ "Xe" ] = [0]
atom_prefered_oxidation[ "Cs" ] = [1, 0]
atom_prefered_oxidation[ "Ba" ] = [2, 0]
atom_prefered_oxidation[ "La" ] = [3, 0]
atom_prefered_oxidation[ "Ce" ] = [3, 0]
atom_prefered_oxidation[ "Pr" ] = [3, 0]
atom_prefered_oxidation[ "Nd" ] = [3, 0]
atom_prefered_oxidation[ "Pm" ] = [3, 0]
atom_prefered_oxidation[ "Sm" ] = [3, 0]
atom_prefered_oxidation[ "Eu" ] = [2, 3, 0]
atom_prefered_oxidation[ "Gd" ] = [3, 0]
atom_prefered_oxidation[ "Tb" ] = [3, 0]
atom_prefered_oxidation[ "Dy" ] = [3, 0]
atom_prefered_oxidation[ "Ho" ] = [3, 0]
atom_prefered_oxidation[ "Er" ] = [3, 0]
atom_prefered_oxidation[ "Tm" ] = [3, 0]
atom_prefered_oxidation[ "Yb" ] = [3, 0]
atom_prefered_oxidation[ "Lu" ] = [3, 0]
atom_prefered_oxidation[ "Hf" ] = [4, 2, 0]
atom_prefered_oxidation[ "Ta" ] = [5, 3, 0]
atom_prefered_oxidation[ "W" ] =  [6, 4, 0]
atom_prefered_oxidation[ "Re" ] = [4, 3, 6, 0]
atom_prefered_oxidation[ "Os" ] = [4, 3, 0]
atom_prefered_oxidation[ "Ir" ] = [4, 3, 0]
atom_prefered_oxidation[ "Pt" ] = [2, 4, 0]
atom_prefered_oxidation[ "Au" ] = [3, 1, 0]
atom_prefered_oxidation[ "Hg" ] = [2, 1, 0]
atom_prefered_oxidation[ "Tl" ] = [3, 1, 0]
atom_prefered_oxidation[ "Pb" ] = [4, 2, -4, 0]
atom_prefered_oxidation[ "Bi" ] = [3, -3, 0]
atom_prefered_oxidation[ "Po" ] = [-2, 0]
atom_prefered_oxidation[ "At" ] = [-1, 0]

for key in keys(atom_prefered_oxidation)
    atom_prefered_oxidation[key*"_d"] = atom_prefered_oxidation[key]
end

for key in keys(atom_prefered_oxidation)
    atom_prefered_oxidation[Symbol(key)] = atom_prefered_oxidation[key]
end


"""
    cutoff_dist::Dict()

Cutoff distance for Hamiltonian calculation
"""
cutoff_dist = Dict()

function get_cutoff(at1, at2)
#    println("asdf")
#    return [25.0, 20.0]
    if (at1,at2) in keys(cutoff_dist)
        return cutoff_dist[(at1,at2)]
    else
        rad1 = atom_radius[at1] / 0.529 / 100.0
        rad2 = atom_radius[at2] / 0.529 / 100.0
        
#        cutoff2X    = (rad1 + rad2) / 2.0 * 7.0
 #       cutoff_onX  = (rad1 + rad2) / 2.0 * 6.0

 #       cutoff2X   = max(min(cutoff2X,   19.01), 15.01) #2body
 #       cutoff_onX = max(min(cutoff_onX, 18.01), 14.51) #onsite

#        cutoff2X    = (rad1 + rad2) / 2.0 * 7.5
#        cutoff_onX  = (rad1 + rad2) / 2.0 * 7.0

#        cutoff2X   = max(min(cutoff2X,   21.01), 18.01) #2body
#        cutoff_onX = max(min(cutoff_onX, 20.51), 17.51) #onsite

        cutoff2X    = (rad1 + rad2) / 2.0 * 6.5
        cutoff_onX  = (rad1 + rad2) / 2.0 * 6.0

        cutoff2X   = max(min(cutoff2X,   19.51), 14.01) #2body
        cutoff_onX = max(min(cutoff_onX, 18.51), 13.01) #onsite

        cutoff2X = cutoff2X-0.25 + 1.5/2 
        cutoff_onX = cutoff_onX-0.25 + 1.5/2 
        
        cutoff_dist[(at1,at2)] = [cutoff2X, cutoff_onX]
        cutoff_dist[(at2,at1)] = [cutoff2X, cutoff_onX]
        
#        return [cutoff2X, cutoff_onX]

        return [cutoff2X, cutoff_onX]        

        #        return [25.0, 20.0]
    end
end

function get_cutoff(at1, at2, at3)
    #return 20.0
    if (at1,at2,at3) in keys(cutoff_dist)
        return cutoff_dist[(at1,at2,at3)]
    else
        rad1 = atom_radius[at1] / 0.529 / 100.0
        rad2 = atom_radius[at2] / 0.529 / 100.0
        rad3 = atom_radius[at3] / 0.529 / 100.0

#        cutoff3bX = minimum([rad1, rad2, rad3])*5.0
#        cutoff3bX  = max(min(cutoff3bX,  13.51),  10.01) #3body

#        cutoff3bX = minimum([rad1, rad2, rad3])*7.0
#        cutoff3bX  = max(min(cutoff3bX,  20.01),  17.01) #3body

        cutoff3bX = minimum([rad1, rad2, rad3])*4.5
        cutoff3bX  = max(min(cutoff3bX,  12.01),  9.01) #3body

        cutoff3bX = cutoff3bX-0.25 + 2.0/2

        
        cutoff_dist[(at1,at2,at3)] = cutoff3bX
        cutoff_dist[(at1,at3,at2)] = cutoff3bX
        cutoff_dist[(at2,at1,at3)] = cutoff3bX
        cutoff_dist[(at2,at3,at1)] = cutoff3bX
        cutoff_dist[(at3,at1,at2)] = cutoff3bX
        cutoff_dist[(at3,at2,at1)] = cutoff3bX
        
        return cutoff3bX
    end
end        


function get_cutoff(at1, at2, at3, at4)
    #return 20.0
    return 9.0
end        

"""
    min_dimer_dist_dict::Dict()

Approximate minimum distance in current fitting data.
"""
min_dimer_dist_dict = Dict()
fil = open("$SRCDIR/mindist.csv")

for line in readlines(fil)
    sp = split(line)
#    if sp[1] == "Li" && sp[2] == "Li"
#        println(sp, " xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
#    end
    min_dimer_dist_dict[(Symbol(sp[1]), Symbol(sp[2]))] = parse(Float64, sp[3])
    min_dimer_dist_dict[(sp[1], sp[2])] = parse(Float64, sp[3])
end


########
"""
    sub_list::Dict()

Most likely substitution list, approximate.
"""
sub_list = Dict()

sub_list[ "H" ] =  ["F","Li"]
sub_list[ "He" ] = []
sub_list[ "Li" ] = ["Na", "H"]
sub_list[ "Be" ] = ["Mg", "B"]
sub_list[ "B" ] =  ["Al", "C"]
sub_list[ "C" ] =  ["Si", "B"]
sub_list[ "N" ] =  ["P", "O"]
sub_list[ "O" ] =  ["S"]
sub_list[ "F" ] =  ["Cl"]
sub_list[ "Ne" ] = []
sub_list[ "Na" ] = ["K", "Li"]
sub_list[ "Mg" ] = ["Be", "Ca", "Zn"]
sub_list[ "Al" ] = ["Ga", "Sc", "Y"]
sub_list[ "Si" ] = ["Ge", "Ti", "C"]
sub_list[ "P" ] = ["As", "N"]
sub_list[ "S" ] = ["Se", "O"]
sub_list[ "Cl" ] = ["F", "Br"]
sub_list[ "Ar" ] = []
sub_list[ "K" ] = ["Rb", "Na"]
sub_list[ "Ca" ] = ["Mg", "Sr"]
sub_list[ "Sc" ] = ["Al", "Ga", "Y"]
sub_list[ "Ti" ] = ["Zr", "Hf"]
sub_list[ "V" ] = ["Nb", "Ti"]
sub_list[ "Cr" ] = ["V", "Mn", "Mo"]
sub_list[ "Mn" ] = ["Cr", "Fe", "Ru"]
sub_list[ "Fe" ] = ["Mn", "Co", "Ru"]
sub_list[ "Co" ] = ["Fe", "Ni", "Rh"]
sub_list[ "Ni" ] = ["Co", "Pd", "Cu"]
sub_list[ "Cu" ] = ["Ni", "Zn", "Ag"]
sub_list[ "Zn" ] = ["Mg", "Cd", "Ca"]
sub_list[ "Ga" ] = ["Al", "In", "Sc"]
sub_list[ "Ge" ] = ["Si", "Sn", "Ti"]
sub_list[ "As" ] = ["Sb", "P", "Ga"]
sub_list[ "Se" ] = ["S", "Te"]
sub_list[ "Br" ] = ["Cl", "I"]
sub_list[ "Kr" ] = []
sub_list[ "Rb" ] = ["Cs", "K"]
sub_list[ "Sr" ] = ["Ca", "Ba", "Pb"]
sub_list[ "Y" ] = ["Sc", "Bi", "Ga"]
sub_list[ "Zr" ] = ["Ti", "Hf", "Ge"]
sub_list[ "Nb" ] = ["Ta", "V", "Mo"]
sub_list[ "Mo" ] = ["Cr", "W", "Nb"]
sub_list[ "Tc" ] = ["Mn", "Re", "Ru"]
sub_list[ "Ru" ] = ["Fe", "Os", "Mn"]
sub_list[ "Rh" ] = ["Co", "Ir", "Re"]
sub_list[ "Pd" ] = ["Ni", "Pt", "Ag"]
sub_list[ "Ag" ] = ["Cu", "Au", "Pd"]
sub_list[ "Cd" ] = ["Zn", "Hg", "Ca"]
sub_list[ "In" ] = ["Ga", "Tl", "Sb"]
sub_list[ "Sn" ] = ["Ge", "Pb", "Sb"]
sub_list[ "Sb" ] = ["As", "Bi", "Sn"]
sub_list[ "Te" ] = ["Se", "Sb"]
sub_list[ "I" ] = ["Br", "Cl"]
sub_list[ "Xe" ] = []
sub_list[ "Cs" ] = ["Rb", "K"]
sub_list[ "Ba" ] = ["Sr", "La"]
sub_list[ "La" ] = ["Y", "Sc", "Ba"]
sub_list[ "Ce" ] = ["Ti", "Zr", "Hf"]
sub_list[ "Pr" ] = ["La",  "Y"]
sub_list[ "Nd" ] = ["La",  "Y"]
sub_list[ "Pm" ] = ["La",  "Y"]
sub_list[ "Sm" ] = ["La",  "Y"]
sub_list[ "Eu" ] = ["Sr", "Ba", "Y"]
sub_list[ "Gd" ] = ["La",  "Y", "Sc"]
sub_list[ "Tb" ] = ["La",  "Y", "Sc"]
sub_list[ "Dy" ] = ["La",  "Y", "Sc"]
sub_list[ "Ho" ] = ["La",  "Y", "Sc"]
sub_list[ "Er" ] = ["La",  "Y", "Sc"]
sub_list[ "Tm" ] = ["La",  "Y", "Sc"]
sub_list[ "Yb" ] = ["La",  "Y", "Sc"]
sub_list[ "Lu" ] = ["La", "Y", "Sc"]
sub_list[ "Hf" ] = ["Zr", "Ti", "Ru"]
sub_list[ "Ta" ] = ["Nb", "V", "W"]
sub_list[ "W" ] = ["Mo", "Re", "Ta", "Cr"]
sub_list[ "Re" ] = ["Tc", "Os", "W", "Mn"]
sub_list[ "Os" ] = ["Ru", "Ir", "Re", "Fe"]
sub_list[ "Ir" ] = ["Rh", "Os", "Pt", "Co"]
sub_list[ "Pt" ] = ["Pd", "Ir", "Au", "Ni"]
sub_list[ "Au" ] = ["Ag", "Pt", "Hg", "Cu"]
sub_list[ "Hg" ] = ["Cd", "Pb", "Zn"]
sub_list[ "Tl" ] = ["In", "Bi", "Hg"]
sub_list[ "Pb" ] = ["Sn", "Sr", "Ba"]
sub_list[ "Bi" ] = ["Sb", "In", "Tl"]
sub_list[ "Po" ] = []
sub_list[ "At" ] = []


"""
    electronegativity::Dict()

Pauling scale.
"""
electronegativity = Dict()

electronegativity["H"] =  2.20
electronegativity["He"] =  0.98
electronegativity["Be"] =  1.57
electronegativity["B"] =  2.04
electronegativity["C"] =  2.55
electronegativity["N"] =  3.04
electronegativity["O"] =  3.44
electronegativity["F"] =  3.98
electronegativity["Ne"] = 0.0
electronegativity["Na"] = 0.93
electronegativity["Mg"] = 1.31
electronegativity["Al"] = 1.61
electronegativity["Si"] = 1.90
#electronegativity["Sid"] = 1.90
electronegativity["P"] = 2.19
electronegativity["S"] = 2.58
electronegativity["Cl"] = 3.16
electronegativity["K"] = 0.82
electronegativity["Ca"] = 1.00
electronegativity["Sc"] = 1.36
electronegativity["Ti"] = 1.54
electronegativity["V"] = 1.63
electronegativity["Cr"] = 1.66
electronegativity["Mn"] = 1.55
electronegativity["Fe"] = 1.83
electronegativity["Co"] = 1.88
electronegativity["Ni"] = 1.91
electronegativity["Cu"] = 1.90
electronegativity["Zn"] = 1.65
electronegativity["Ga"] = 1.81
electronegativity["Ge"] = 2.01
electronegativity["Ged"] = 2.01
electronegativity["As"] = 2.18
electronegativity["Se"] = 2.55
electronegativity["Br"] = 2.96
electronegativity["Kr"] = 3.00
electronegativity["Rb"] = 0.82
electronegativity["Sr"] = 0.95
electronegativity["Y"] = 1.22
electronegativity["Zr"] = 1.33
electronegativity["Nb"] = 1.6
electronegativity["Mo"] = 2.16
electronegativity["Tc"] = 1.9
electronegativity["Ru"] = 2.2
electronegativity["Rh"] = 2.28
electronegativity["Pd"] = 2.20
electronegativity["Ag"] = 1.93
electronegativity["Cd"] = 1.69
electronegativity["In"] = 1.78
electronegativity["Sn"] = 1.96
electronegativity["Sb"] = 2.05
electronegativity["Te"] = 2.1
electronegativity["I"] = 2.66
electronegativity["Xe"] = 2.6
electronegativity["Cs"] = 0.79
electronegativity["Ba"] = 0.89
electronegativity["La"] = 1.10
electronegativity["Ce"] = 1.12
electronegativity["Pr"] = 1.13
electronegativity["Nd"] = 1.14
electronegativity["Pm"] = 1.13
electronegativity["Sm"] = 1.17
electronegativity["Eu"] = 1.2
electronegativity["Gd"] = 1.2
electronegativity["Tb"] = 1.22
electronegativity["Dy"] = 1.23
electronegativity["Ho"] = 1.24
electronegativity["Er"] = 1.24
electronegativity["Tm"] = 1.25
electronegativity["Yb"] = 1.1
electronegativity["Lu"] = 1.27
electronegativity["Hf"] = 1.3
electronegativity["Ta"] = 1.5
electronegativity["WT"] = 2.36
electronegativity["Re"] = 1.9
electronegativity["Os"] = 2.2
electronegativity["Ir"] = 2.2
electronegativity["Pt"] = 2.28
electronegativity["Au"] = 2.54
electronegativity["Hg"] = 2.00
electronegativity["Tl"] = 1.62
electronegativity["Pb"] = 2.33
electronegativity["Bi"] = 2.02
electronegativity["Po"] = 2.0
electronegativity["At"] = 2.2
electronegativity["Fr"] = 0.7
electronegativity["Ra"] = 0.89
electronegativity["Ac"] = 1.1
electronegativity["Th"] = 1.3
electronegativity["Pa"] = 1.5
electronegativity["U"] = 1.38
electronegativity["Np"] = 1.36
electronegativity["Pu"] = 1.28
electronegativity["Am"] = 1.3
electronegativity["Cm"] = 1.3
electronegativity["Bk"] = 1.3
electronegativity["Cf"] = 1.3
electronegativity["Es"] = 1.3


for key in keys(electronegativity)
    electronegativity[key*"_d"] = electronegativity[key]
end

formation_energy_ref = Dict()

formation_energy_ref["Hx"] =  -0.2412200680378065
formation_energy_ref["X"] =  0.0
formation_energy_ref["Xa"] =  0.0

formation_energy_ref["Cl"] =  -0.13497163207243545
formation_energy_ref["Al"] =  -0.29091361813767413
formation_energy_ref["Be"] =  -0.29201079531889107
formation_energy_ref["Re"] =  -0.8883644117346989
formation_energy_ref["Cr"] =  -0.6989093821143229
formation_energy_ref["Na"] =  -0.09691800599890144
formation_energy_ref["Sb"] =  -0.3175643486870854
formation_energy_ref["Ni"] =  -0.4364163406820012
formation_energy_ref["S"] =  -0.30197576430356793
formation_energy_ref["Ru"] =  -0.6982894041532859
formation_energy_ref["Tl"] =  -0.18202296106778704
formation_energy_ref["W"] =  -0.8747898540244421
formation_energy_ref["Zr"] =  -0.5573338529326719
formation_energy_ref["F"] =  -0.13426498419775612
formation_energy_ref["Co"] =  -0.55079152086779
formation_energy_ref["Rh"] =  -0.5191278986005443
formation_energy_ref["Ca"] =  -0.15552687954944133
formation_energy_ref["N"] =  -0.6158937937176567
formation_energy_ref["As"] =  -0.3575026746186296
formation_energy_ref["Se"] =  -0.2674979594162418
formation_energy_ref["Y"] =  -0.3527326401211184
formation_energy_ref["Pt"] =  -0.4927545253753749
formation_energy_ref["I"] =  -0.11761828669794738
formation_energy_ref["Fe"] =  -0.6295056634305922
formation_energy_ref["Ba"] =  -0.15643118139661283
formation_energy_ref["Hf"] =  -0.5592348224131172
formation_energy_ref["C"] =  -0.6939900829834258
formation_energy_ref["Zn"] =  -0.11525323820592348
formation_energy_ref["Ir"] =  -0.6963577685116604
formation_energy_ref["Rb"] =  -0.06830710806754325
formation_energy_ref["In"] =  -0.2074761713795965
formation_energy_ref["Hg"] =  -0.04145676412692012
formation_energy_ref["Te"] =  -0.243887173445439
formation_energy_ref["Bi"] =  -0.29865804886640035
formation_energy_ref["Cu"] =  -0.31327806206661535
formation_energy_ref["Tc"] =  -0.7929095477392991
formation_energy_ref["Sn"] =  -0.2937190504709122
formation_energy_ref["Mo"] =  -0.8260038917015606
formation_energy_ref["Ta"] =  -0.7460229038179591
formation_energy_ref["Li"] =  -0.13932628859861396
formation_energy_ref["Cd"] =  -0.08787569055060374
formation_energy_ref["Os"] =  -0.8517405771113147
formation_energy_ref["Ti"] =  -0.5119888847399542
formation_energy_ref["V"] =  -0.6349543354131981
formation_energy_ref["Si"] =  -0.40754479040813507
formation_energy_ref["Sid"] =  -0.40754479040813507
formation_energy_ref["Ga"] =  -0.23449803933579005
formation_energy_ref["Au"] =  -0.28341301745817304
formation_energy_ref["Mg"] =  -0.126284990107294
formation_energy_ref["Ag"] =  -0.23657551086864714
formation_energy_ref["K"] =  -0.0761247994147709
formation_energy_ref["Sc"] =  -0.349825562961243
formation_energy_ref["Ge"] =  -0.3466031556989151
formation_energy_ref["Ged"] =  -0.3466031556989151
formation_energy_ref["Nb"] =  -0.737561493627723
formation_energy_ref["Pd"] =  -0.3270244478002269
formation_energy_ref["Sr"] =  -0.13365669688262471
formation_energy_ref["Br"] =  -0.12357029567390754
formation_energy_ref["Pb"] =  -0.27802998869077555
formation_energy_ref["P"] =  -0.4052925946101844
formation_energy_ref["La"] =  -0.3659101249542971
formation_energy_ref["H"] =  -0.2412200680378065
formation_energy_ref["Cs"] = -0.06273440204358138
formation_energy_ref["Mn"] = -0.6713694611146082
formation_energy_ref["B"] = -0.5109718427571567
formation_energy_ref["O"] = -0.3767033842534424

for key in keys(formation_energy_ref)
    formation_energy_ref[key*"_d"] = formation_energy_ref[key]
end


for key in keys(formation_energy_ref)
    formation_energy_ref[Symbol(key)] = formation_energy_ref[key]
end

uniform_charge_interaction = Dict()

uniform_charge_interaction[:Hx ] = 0.0
uniform_charge_interaction[:X] = 0.0
uniform_charge_interaction[:Xa] = 0.0




for key in keys(uniform_charge_interaction)
    uniform_charge_interaction[String(key)] = uniform_charge_interaction[key]
end


data_temp = Dict()

data_temp["Ag"] =  [0.47609387973587514, 67.34725662242677]
data_temp["Al"] =  [0.4325911327022497, 82.2977471332531]
data_temp["As"] =  [0.5673213457189965, 47.83363758199629]
data_temp["Au"] =  [0.5272682730980324, 57.433069219691035]
data_temp["Ba"] =  [0.30727027800798373, 161.02090516624165]
data_temp["B"] =  [0.6381791968250066, 37.92911684769323]
data_temp["Be"] =  [0.5595157097364083, 47.92369499829685]
data_temp["Bi"] =  [0.4810473708551888, 65.74425467396907]
data_temp["Br"] =  [0.6487058021192891, 39.468124375087164]
data_temp["Ca"] =  [0.3676356073173389, 109.33714260288636]
data_temp["C"] =  [0.7702800011582039, 26.88555738793191]
data_temp["Cd"] =  [0.5192534187428902, 52.20888007549389]
data_temp["Cl"] =  [0.7293490347269327, 32.49608969334515]
data_temp["Co"] =  [0.4480666602798483, 83.39695660660308]
data_temp["Cr"] =  [0.416200320930515, 93.42941242360374]
data_temp["Cs"] =  [0.23659354875651875, 258.69326061187024]
data_temp["Cu"] =  [0.5008718896009546, 62.75229850290244]
data_temp["Fe"] =  [0.4368524489783394, 86.81408809295716]
data_temp["F"] =  [1.090821815297757, 17.626602941047544]
data_temp["Ga"] =  [0.4393075771923728, 79.2817865698265]
data_temp["Ge"] =  [0.5063737860332831, 60.1392072451489]
data_temp["Ged"] =  [0.5063737860332831, 60.1392072451489]
data_temp["Hf"] =  [0.5367826942470761, 51.64230896252188]
data_temp["Hg"] =  [0.5617844607440594, 43.50181533015791]
data_temp["H"] =  [0.8714837907568558, 18.987364618334517]
data_temp["I"] =  [0.5580268944626001, 51.30394763510125]
data_temp["In"] =  [0.40131032313960385, 94.75727949689453]
data_temp["Ir"] =  [0.45007920455852224, 80.6324472959521]
data_temp["K"] =  [0.2660316113917355, 205.0753192185277]
data_temp["La"] =  [0.4415798622981549, 78.59884529779893]
data_temp["Li"] =  [0.35031914597768776, 122.60257453578949]
data_temp["Mg"] =  [0.4641421511158119, 69.1223608908259]
data_temp["Mn"] =  [0.42071488211907343, 92.63302148648297]
data_temp["Mo"] =  [0.38263405630146924, 108.9490122327217]
data_temp["Na"] =  [0.3345987589320675, 131.7729342308781]
data_temp["Nb"] =  [0.3958741054351723, 101.17836851462616]
data_temp["Ni"] =  [0.46015558368000153, 79.78753919944495]
data_temp["N"] =  [0.9042974802119927, 20.367444331668707]
data_temp["O"] =  [1.0145543956989567, 17.765554012950126]
data_temp["Os"] =  [0.47466543877005807, 71.78916395707896]
data_temp["Pb"] =  [0.44070374135225493, 79.34400848607586]
data_temp["Pd"] =  [0.7146457120693652, 31.6515380931471]
data_temp["P"] =  [0.6023859243119037, 42.4910114468705]
data_temp["Pt"] =  [0.49388657173022066, 68.5656732488334]
data_temp["Rb"] =  [0.25908611548503346, 212.1609742850391]
data_temp["Re"] =  [0.6041979328728289, 42.67488912748284]
data_temp["Rh"] =  [0.6050159245280038, 44.595584480097735]
data_temp["Ru"] =  [0.48033579694973716, 70.3031106451712]
data_temp["Sb"] =  [0.5016010922908432, 60.637532394985215]
data_temp["Sc"] =  [0.5873089033731779, 42.726119843127684]
data_temp["Se"] =  [0.6233519032682732, 39.53330339534418]
data_temp["Si"] =  [0.5203977905150886, 56.8314240643758]
data_temp["Sid"] =  [0.5203977905150886, 56.8314240643758]
data_temp["Sn"] =  [0.45588030411946134, 73.83127678236619]
data_temp["Sr"] =  [0.3396906409023482, 126.76704291627323]
data_temp["S"] =  [0.6808063801727222, 33.850523636729974]
data_temp["Ta"] =  [0.5834697574029772, 44.1441130084451]
data_temp["Tc"] =  [0.3843221409906412, 109.13230376005151]
data_temp["Te"] =  [0.5480878800622282, 49.44330339099042]
data_temp["Ti"] =  [0.41845981847129954, 91.19261470006593]
data_temp["Tl"] =  [0.37541685512417644, 108.36528179837865]
data_temp["V"] =  [0.4092597927158933, 96.0614722016224]
data_temp["W"] =  [0.6227793180705922, 38.35008067706342]
data_temp["Y"] =  [0.525144790638078, 52.08026894508406]
data_temp["Zn"] =  [0.560635872916608, 47.117192678427145]
data_temp["Zr"] =  [0.5493393488752593, 49.14643018003017]

for key in keys(data_temp)
    data_temp[key*"_d"] = data_temp[key]
end


function prepare_atom(A, U, uni)    
    if !(String(A) in keys(atoms))
        return
    end
    
    old = atoms[String(A) ]
    atoms[String(A) ] =  makeatom(String(A) ,old.Z, old.row, old.col, old.mass,   old.nval,   old.nsemicore,
                                  old.orbitals, old.total_energy, 13.6057*[old.eigs[o] for o in old.orbitals],     0.0 ,U)
    
    atoms[Symbol(A)] = atoms[ String(A) ]
    atoms[String(A)] = atoms[ Symbol(A) ]

    uniform_charge_interaction[String(A)] = uni
    uniform_charge_interaction[Symbol(A)] = uni

end

for key in keys(data_temp)
    prepare_atom(key, data_temp[key][1], data_temp[key][2])
end


end #ends module
 
