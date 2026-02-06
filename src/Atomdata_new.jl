"""
    module Atomdata

Holds the data for isolated atoms.
"""
module Atomdata
using ..AtomicMod:makeatom
using ..ThreeBodyTB:SRCDIR
#using ..Atomic:AtomicMod


"""
    atoms::Dict()

Periodic table information.
"""
atoms = Dict()

function load_atoms()
#      name                     name  Z  row col mass nval nsemi  orbitals total_energy    e_orb1, e_orb2
#                                                (amu)                      (Ryd)             (eV)
#
#atoms["T"] = makeatom("T",    0, 0, 0, 0.0,  1.0,   0,   [:s], 0.0, [0.0], 0.0);

#atoms["Is"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p], -9.15243175, [-10.6594, -3.9477],0.007075865418, 0.5300262399999837);

    atoms["Hx" ] = makeatom("Hx" , 100, 1, 1, 1.0079,   1.0,   0,   [:s ],     -0.90531367, [-6.3233], 0.0, 0.9194820210526328);
    atoms["X"] = makeatom("X",    101, 0, 0, 0.0,  1.0,   0,   [:s], 0.0, [0.0], 0.0, 0.5);
    atoms["Xa"] = makeatom("Xa",    102, 0, 0, 0.0,  1.0,   0,   [:s], 3.0, [3.0* 13.605693122], 0.0, 0.5);


    
    atoms["H" ] = makeatom("H" , 1, 1, 1, 1.0079,   1.0,   0,   [:s ],     -0.90531367, [-0.46533567834863526], 0.0, 0.8451632950940942, adj= 15.56325018071482 );
    atoms["Li"] = makeatom("Li", 3, 2, 1, 6.941,    1.0,   2,   [:s, :p], -14.23107803, [-0.20724077061802793; -0.07697938341878834],0.0, 0.3902587823285621, adj=75.628523583303  , semicore_orbitals=[:s], semicore_eigs=[-3.77612128609209] );
    atoms["Be"] = makeatom("Be", 4, 2, 2, 9.012,    2.0,   2,   [:s, :p], -27.80694399, [-0.4073529309633574; -0.14464734287229414],0.0, 0.5751048209185653, adj= 70.89149364003049  , semicore_orbitals=[:s], semicore_eigs=[-7.7495366938864825] );
    atoms["B" ] = makeatom("B" , 5, 2, 3, 10.81,    3.0,   0,   [:s, :p], -5.84632720, [-0.6875808938239544; -0.2613131797530708], 0.0, 0.6267943641604508, adj= 67.62585046150113 , fullU=true, Uarr = [0.616882180634314, 0.6566090373908344, 0.6245206094692716] );
    atoms["C" ] = makeatom("C" , 6, 2, 4, 12.011,   4.0,   0,   [:s, :p], -10.78129974, [-1.002601890640717; -0.3846551855250298],0.0, 0.7978153170555154, adj= 60.516015221894946 , fullU=true, Uarr = [0.7423589123939367, 0.7982248012476293, 0.7641398318704328]  );
    atoms["N" ] = makeatom("N" , 7, 2, 5, 14.007,   5.0,   0,   [:s, :p], -19.47573456, [-1.3549512787495301; -0.5173659349266165],0.0, 0.8874159249179556, adj= 54.39724301582911  );
    atoms["O" ] = makeatom("O" , 8, 2, 6, 15.999,   6.0,   0,   [:s, :p], -31.71324769, [-1.7468475138908586; -0.6600313069986147],0.0, 1.015623980577927, adj= 48.942302010664385  );
    atoms["F" ] = makeatom("F" , 9, 2, 7, 18.998,   7.0,   0,   [:s, :p], -48.22248233, [-2.184368257560789; -0.8122569224622285],0.0, 1.1398170857521994, adj= 44.53471750938858  );
    atoms["Na"] = makeatom("Na", 11, 3, 1, 22.98,  1.0,   8,   [:s, :p], -95.06749946, [-0.19844190728525501; -0.05218815323693155], 0.0, 0.38560892349289677, adj= 108.92841309130345  , semicore_orbitals=[:s, :p], semicore_eigs=[-4.157649733166457, -2.105445706578367] );
    atoms["Mg"] = makeatom("Mg", 12, 3, 2, 24.31,  2.0,   8,   [:s, :p], -124.99127898, [-0.34372480877039063; -0.09551041814912471], 0.0, 0.5106562587165097, adj= 118.00457569316572  , semicore_orbitals=[:s, :p], semicore_eigs=[-5.849358490948838; -3.420898450359709] );
    atoms["Al"] = makeatom("Al", 13, 3, 3, 26.98,  3.0,   0,   [:s, :p],  -6.37406439, [-0.5669635647556415; -0.19757116543256825], 0.0, 0.47560947595799796 , adj= 120.35327370972065  , fullU=true, Uarr = [0.331290464808923, 0.4835827995950456, 0.42998890690669406]);
    atoms["Si"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p], -9.15243181, [-0.7914241887679531; -0.2981179964824184],0.0, 0.5597222300055447, adj= 124.19186186707114,  fullU=true, Uarr = [0.4246166888157803, 0.6103690583467152, 0.5612372485699727]);
    atoms["P" ] = makeatom("P" , 15, 3, 5, 30.97,  5.0,   0,   [:s, :p], -15.07327848, [-1.022828949354686; -0.40267725947367494], 0.0, 0.6422064900805308, adj= 120.25412027683043  );
    atoms["S" ] = makeatom("S" , 16, 3, 6, 32.06,  6.0,   0,   [:s, :p], -23.70562681, [-1.2637282132428946; -0.5128739318001997], 0.0, 0.7079117127047072, adj= 114.81807427734294  );
    atoms["Cl"] = makeatom("Cl", 17, 3, 7, 35.45,  7.0,   0,   [:s, :p], -33.06905525, [-1.5153800694575343; -0.6292474932451964], 0.0, 0.7476423278619385, adj= 109.3609934532405  );
    atoms["K" ] = makeatom("K" , 19, 4, 1, 39.098,  1.0,   8,   [:s, :d, :p], -56.84482202 , [-0.1688811278691164; -0.01698856198985435; -0.056855390138818754], 0.0, 0.37734366682102694, adj= 192.31468038720865  , semicore_orbitals=[:s, :p], semicore_eigs=[-2.5919695598282306; -1.3811606089650703] );
    atoms["Ca"] = makeatom("Ca", 20, 4, 2, 40.078,  2.0,   8,   [:s, :d, :p], -74.68919704, [ -0.27660048821000915; -0.14361607141470636; -0.1014606571098085 ], 0.0, 0.4749594731162585, adj= 213.27385795261185  , semicore_orbitals=[:s, :p], semicore_eigs=[-3.451119954305641; -2.0540300136124428] );
    atoms["Sc"] = makeatom("Sc", 21, 4, 2.05, 44.956,  3.0,   8,   [:s, :d, :p], -93.89389249,  [ -0.30780962808922513; -0.2385703068083227; -0.10646754352172771],  0.0, 0.7781746286900919 , adj= 202.88438490397112  , semicore_orbitals=[:s, :p], semicore_eigs=[-4.033242049845243; -2.4657682331409854]);
    atoms["Ti"] = makeatom("Ti", 22, 4, 2.15, 47.867,  4.0,   8,   [:s, :d, :p], -118.63461912, [-0.32668593549677594; -0.30600003403662074; -0.10501413617737301],  0.0, 0.5344202735274272, adj= 192.39001137597768  , semicore_orbitals=[:s, :p], semicore_eigs=[-4.579043001286401; -2.838828500088978]);
    atoms["V" ] = makeatom("V" , 23, 4, 2.25, 50.941,  5.0,    8,   [:s, :d, :p], -144.34826807, [-0.3280477165418133; -0.31699982246372754; -0.09598701231320347], 0.0, 0.5051179562238686, adj= 179.1050259538386  , semicore_orbitals=[:s, :p], semicore_eigs=[-5.063689993395585; -3.149763759434965] );
    atoms["Cr"] = makeatom("Cr", 24, 4, 2.35, 51.996,  6.0,    8,   [:s, :d, :p], -174.69181247, [-0.3277707801780594; -0.3217541659915794; -0.08535475083154154], 0.0, 0.5023896126040788 , adj= 168.13548531988167  , semicore_orbitals=[:s, :p], semicore_eigs=[-5.561070838386585; -3.4657325795837184]);
    atoms["Mn"] = makeatom("Mn", 25, 4, 2.45, 54.938,  7.0,    8,   [:s, :d, :p], -210.40629067, [-0.3276990460927028; -0.3256346732421914;-0.07750632701308155], 0.0, 0.5085605056527054, adj= 157.92334663274937  , semicore_orbitals=[:s, :p], semicore_eigs=[-6.078058704646573; -3.794718726750049]);
    atoms["Fe"] = makeatom("Fe", 26, 4, 2.55, 55.845,  8.0,    8,   [:s, :d, :p], -249.56808439, [-0.3282171355958345; -0.32966680672920284; -0.06953166509075472], 0.0, 0.5138399863052884 , adj= 148.86501881965447  , semicore_orbitals=[:s, :p], semicore_eigs=[-6.616005085303623; -4.136819859230937]);
    atoms["Co"] = makeatom("Co", 27, 4, 2.65, 58.993,  9.0,    8,   [:s, :d, :p], -297.72249071, [-0.3285562624529606; -0.33353263126199884; -0.05955024442762707], 0.0, 0.5205632537537023, adj= 140.21772403677664  , semicore_orbitals=[:s, :p], semicore_eigs=[-7.177473034997793; -4.495538901683293] );
    atoms["Ni"] = makeatom("Ni", 28, 4, 2.75, 58.963,  10.0,   8,   [:s, :d, :p], -342.65159683, [-0.32888432962862346;-0.33818091690002905; -0.049199552300187695], 0.0, 0.530332541897714, adj= 134.1024773997327  , semicore_orbitals=[:s, :p], semicore_eigs=[0.0, 0.0] );
    atoms["Cu"] = makeatom("Cu", 29, 4, 2.85, 63.546,  11.0,   8,   [:s, :d, :p], -403.19763559, [-0.3392595220980112; -0.371639309702241; -0.04783574462144305], 0.0, 0.5619138137847363, adj= 128.09549628995777  , semicore_orbitals=[:s, :p], semicore_eigs=[-8.40281251211807; -5.287392185939331]);
    atoms["Zn"] = makeatom("Zn", 30, 4, 2.95, 65.38,   12.0,   8,   [:s, :d, :p], -461.12493484, [-0.44079883980156; -0.747314944736136; -0.07877751429918178], 0.0, 0.6233810643339528, adj= 134.69292515789127  , semicore_orbitals=[:s, :p], semicore_eigs=[-9.475588017870018; -6.118805445036133] );
    atoms["Ga"] = makeatom("Ga", 31, 4, 3.0,  69.72,   13.0,   6,   [:s,:d, :p], -414.75774465, [-0.6587239854030615; -1.4072063981803296; -0.18873207882392629], 0.0, 0.5103937462614931, adj= 157.6873480049696  , semicore_orbitals=[:p], semicore_eigs=[-7.250669797283856] );
    atoms["Ge"] = makeatom("Ge", 32, 4, 4.0,  72.63,   4.0,    10,   [:s, :p], -212.84384318, [-0.8631076406670051; -0.2864129800754766], 0.0, 0.5643618220580027, adj=155.72365698047153  , semicore_orbitals=[:d], semicore_eigs=[-2.1543465449641492] );
    atoms["As"] = makeatom("As", 33, 4, 5.0,  74.92,   5.0,    0,   [:s, :p], -39.62501000, [-1.0658794604198387; -0.3818226009590801], 0.0, 0.6154159430456763, adj=132.70436281983118 );
    atoms["Se"] = makeatom("Se", 34, 4, 6.0,  78.97,   6.0,    0,   [:s, :p],  -43.03651432, [-1.270861842892195; -0.4783993632527501], 0.0, 0.6730297191236639, adj=135.46595175052084 );
    atoms["Br"] = makeatom("Br", 35, 4, 7.0,  79.90,   7.0,    0,   [:s, :p], -40.38769338, [-1.4797235464586467; -0.5773254243409252], 0.0, 0.7295557399434305, adj=136.2783003205783 );
    atoms["Rb"] = makeatom("Rb", 37, 5, 1, 85.468,  1.0,   8,   [:s, :d, :p], -53.05281238, [-0.1630615175904313;  -0.01705671661180152; -0.05227246123519389], 0.0, 0.38651513892406236, adj= 234.15893922191282  , semicore_orbitals=[:s, :p], semicore_eigs=[2.3445063882389627; -1.1740025289733453] );
    atoms["Sr"] = makeatom("Sr", 38, 5, 2, 87.62,  2.0,   8,   [:s, :d, :p], -69.92604259, [-0.25989366680607157; -0.09471289698165973; -0.09305917831962837], 0.0, 0.47280011784661014, adj= 265.3720178977989  , semicore_orbitals=[:s, :p], semicore_eigs=[-3.004109150719295; -1.6763528444073972] );
    atoms["Y" ] = makeatom("Y" , 39, 5, 2.05, 88.906,  3.0,   8,    [:s, :d, :p], -92.30423365,  [-0.30273360607984623; -0.18544528189211198; -0.1043094162182441], 0.0, 0.6384727135899312, adj= 262.03087938542785  , semicore_orbitals=[:s, :p], semicore_eigs=[-3.524199243506504; -2.0496495370639223] );
    atoms["Zr"] = makeatom("Zr", 40, 5, 2.15, 91.224,  4.0,   8,    [:s, :d, :p], -98.57421435,  [-0.32798934632854626; -0.26529614203587987; -0.10642044505922535], 0.0, 0.6903649183604765, adj= 255.09572313090078  , semicore_orbitals=[:s, :p], semicore_eigs=[-4.002140654495196; -2.383358632703009] );
    atoms["Nb"] = makeatom("Nb", 41, 5, 2.25, 92.906,  5.0,    8,   [:s, :d, :p], -117.19066274 ,  [-0.34010414026997; -0.3259260270600244; -0.10568400063639913],0.0, 0.5348560719045159, adj= 243.85653587280657  , semicore_orbitals=[:s, :p], semicore_eigs=[-4.443148381640511; -2.681053314293842]);
    atoms["Mo"] = makeatom("Mo", 42, 5, 2.35, 95.95,   6.0,    8,   [:s, :d, :p], -138.51328928,  [-0.3204435410778702; -0.3171576541432412; -0.08697670943079841], 0.0, 0.5098857597525097, adj= 225.56472882250577  , semicore_orbitals=[:s, :p], semicore_eigs=[-4.783904195672978; -2.881575995588024]);
    atoms["Tc"] = makeatom("Tc", 43, 5, 2.45, 98.0,    7.0,    8,   [:s, :d, :p], -173.39188188,  [-0.2969074909559409; -0.3009468981457615; -0.06656179029914151], 0.0, 0.4983977286525886, adj= 212.40279418247877  , semicore_orbitals=[:s, :p], semicore_eigs=[-5.125131347104067; -3.078256700043885]);
    atoms["Ru"] = makeatom("Ru", 44, 5, 2.55, 101.07,  8.0,    8,   [:s, :d, :p], -193.54935366,  [-0.27141582810651993;-0.2829153530978795; -0.04215544281902223], 0.0, 0.511729036259258, adj= 197.64473652496974  , semicore_orbitals=[:s, :p], semicore_eigs=[-5.473293479583835; -3.278420332118144] );
    atoms["Rh"] = makeatom("Rh", 45, 5, 2.65, 102.91,  9.0,    6,   [:s, :d, :p], -170.64161619,  [-0.24960384464558694; -0.27220038174558514; -0.02536252748592261], 0.0, 0.6907390586058985, adj= 175.02929050347564  , semicore_orbitals=[ :p], semicore_eigs=[-3.49302240541897]);
    atoms["Pd"] = makeatom("Pd", 46, 5, 2.75, 106.42,  10.0,   6,   [:s, :d, :p], -198.31288985,  [-0.24215163190410655; -0.2975066758725801; -0.020638263301147243], 0.0, 0.7801811214471324, adj= 168.08982023680983  , semicore_orbitals=[:p], semicore_eigs=[-3.758625872558585] );
    atoms["Ag"] = makeatom("Ag", 47, 5, 2.85, 107.87,  11.0,   8,   [:s, :d, :p], -295.31434496,  [-0.32510984035597895; -0.5478269244677325; -0.04931902836857995], 0.0, 0.5706935927525575, adj= 189.53628807305094  , semicore_orbitals=[:s, :p], semicore_eigs=[-6.979959293823642; -4.300258857682814] );
    atoms["Cd"] = makeatom("Cd", 48, 5, 2.95, 112.41,  12.0,   0,   [:s, :d, :p], -114.39800378,  [-0.4192879192482951; -0.8584546443356645; -0.0841336333019938], 0.0, 0.6005115028052722, adj=166.23163667773144);
    atoms["In"] = makeatom("In", 49, 5, 3.0,  114.82,   13.0,   0,   [:s,:d, :p], -132.97861648 , [-0.6067611541367017; -1.3658289819045248; -0.1851901155227374], 0.0, 0.49245976583127, adj=199.74865952500107);
    atoms["Sn"] = makeatom("Sn", 50, 5, 4.0,  118.71,   4.0,   10,   [:s, :p], -158.61463088, [-0.7788143666960191; -0.2725359483854868], 0.0, 0.5456114086239429, adj=215.85406673518847, semicore_orbitals=[:d], semicore_eigs=[-1.8945415352873742]);
    atoms["Sb"] = makeatom("Sb", 51, 5, 5.0,  121.76,   5.0,   10,   [:s, :p], -184.28336729, [-0.9464502611207213; -0.3560830978948372], 0.0, 0.5973377632845377, adj=221.5444021383691, semicore_orbitals=[:d], semicore_eigs=[-2.4589362341980237] );
    atoms["Te"] = makeatom("Te", 52, 5, 6.0,  127.60,   6.0,    0,   [:s, :p], -26.19141773, [-1.1136466400725715; -0.43862775966463713], 0.0, 0.6243066065034633, adj= 179.86056670213443  );
    atoms["I" ] = makeatom("I" , 53, 5, 7.0,  126.90,   7.0,    0,   [:s, :p], -65.01087258 , [-1.2823617761081663; -0.5215959483847197], 0.0, 0.6675305924093194, adj= 184.02187318011497  );
    atoms["Cs"] = makeatom("Cs", 55, 6, 1, 132.9,  1.0,   8,   [:s, :d, :p ], -62.83174426, [-0.15276252381950234; -0.03313619073851424; -0.051003702149173005], 0.0, 0.4013353474699838, adj= 295.04867135494004  , semicore_orbitals=[:s, :p], semicore_eigs=[-1.9682409685060689; -0.9950327138502739] );
    atoms["Ba"] = makeatom("Ba", 56, 6, 2, 139.3,  2.0,   8,   [:s, :d, :p], -70.01638694, [-0.23845192186212916; -0.14893061829083185; -0.08964664172218782], 0.0, 0.48617114830668345, adj= 344.03101429219186  , semicore_orbitals=[:s, :p], semicore_eigs=[-2.4767733696767706; -1.3772701526786062] );
    atoms["La"] = makeatom("La", 57, 6, 2.05, 138.9,  3.0,   8,    [:s, :d], -101.82206370,  [-0.2678101223286892; -0.2108445837481361], 0.0, 0.5816142601437702, adj= -0.00010378853535193782  , semicore_orbitals=[:s, :p], semicore_eigs=[0.0, 0.0] );
    atoms["Hf"] = makeatom("Hf", 72, 6, 2.15, 178.5,  4.0,   8,    [:s, :d, :p], -157.97552412,  [-0.3777034454259773; -0.20031069134359342; -0.10434007566183431], 0.0, 0.6517565809524214, adj= 238.99056835797487  , semicore_orbitals=[:s, :p], semicore_eigs=[-4.925134766003607; -2.615494246507261]);
    atoms["Ta"] = makeatom("Ta", 73, 6, 2.25, 181.0,  5.0,    8,   [:s, :d, :p], -141.00656328,  [-0.4019620838141361; -0.2660686877089561; -0.10016055709111202], 0.0, 0.6912754847554329, adj= 241.19196583759157  , semicore_orbitals=[:s, :p], semicore_eigs=[-5.367137429486504; -2.8971963868226824]);
    atoms["W" ] = makeatom("W" , 74, 6, 2.35, 183.8,  6.0,    8,   [:s, :d, :p], -158.16868802,  [-0.4222789939236033; -0.33165315512810734; -0.0977666949109555], 0.0, 0.7326044378416148, adj= 239.69638576729054  , semicore_orbitals=[:s, :p], semicore_eigs=[-5.813834624981691; -3.1769206362914155]);
    atoms["Re"] = makeatom("Re", 75, 6, 2.45, 186.2,  7.0,    8,   [:s, :d, :p], -189.48871022,  [-0.4395858849952568; -0.39654450031042543; -0.0943213729593109], 0.0, 0.7652760081244261, adj= 237.0845805996611  , semicore_orbitals=[:s, :p], semicore_eigs=[-6.268172789309192; -3.456147518758882]);
    atoms["Os"] = makeatom("Os", 76, 6, 2.55, 190.2,  8.0,    8,   [:s, :d, :p], -197.71975092,  [-0.4433567666290275; -0.43663725573432416; -0.08611580107480214], 0.0, 0.6076951539696345, adj= 232.57696602665388  , semicore_orbitals=[:s, :p], semicore_eigs=[-6.693284916155786; -3.7023518851575727]);
    atoms["Ir"] = makeatom("Ir", 77, 6, 2.65, 192.2,  9.0,    6,   [:s, :d, :p], -180.96527345,  [-0.41854304283231303; -0.42250421864683285; -0.06695728214496], 0.0, 0.5861813291765662, adj= 214.97667658413565  , semicore_orbitals=[ :p], semicore_eigs=[-3.8735016110543232]);
    atoms["Pt"] = makeatom("Pt", 78, 6, 2.75, 195.1,  10.0,   6,   [:s, :d, :p], -210.06169163,  [-0.3999163150830427; -0.41321823251630246; -0.05151491591403237], 0.0, 0.5839433031650294, adj= 206.58243148535067  , semicore_orbitals=[ :p], semicore_eigs=[-4.065106634197381]);
    atoms["Au"] = makeatom("Au", 79, 6, 2.85, 197.0,  11.0,   0,   [:s, :d, :p], -114.46631173,  [-0.4263323678261778; -0.5073292421824291; -0.053900174732151], 0.0, 0.611432813778228 , adj= 170.12046505665774  );
    atoms["Hg"] = makeatom("Hg", 80, 6, 2.95, 200.6,  12.0,   0,   [:s, :d, :p], -107.48183748,  [-0.5046095642903924; -0.7261634926300994; -0.07256623070784395], 0.0, 0.6495200693187368, adj= 178.3552501095153);
    atoms["Tl"] = makeatom("Tl", 81, 6, 3.0,  204.4,   13.0,   0,   [:s, :d, :p], -144.53037033, [-0.7064249187404431; -1.1389479900265096; -0.17473822050028032], 0.0, 0.50078333374623, adj= 223.04670714874857);
    atoms["Pb"] = makeatom("Pb", 82, 6, 4.0,  207.2,   4.0,   10,   [:s, :p], -163.61141039, [-0.8890428995019919; -0.2584955035922471], 0.0, 0.5488505326347655, adj= 241.63945158822114  , semicore_orbitals=[:d], semicore_eigs=[-1.5554213660160645]);
    atoms["Bi"] = makeatom("Bi", 83, 6, 5.0,  209.0,   5.0,   10,   [:s, :p], -184.68674216, [-1.0683520310241112; -0.3379213594323034], 0.0, 0.597524848429184, adj= 253.30476308498328  , semicore_orbitals=[:d], semicore_eigs=[-1.9889042834093047]);

                           


    atoms["B_d" ] = makeatom("B" , 5, 2, 3, 10.81,    3.0,   0,   [:s, :p, :d], -5.84632720, [-0.68760839796914; -0.26133917512872396; 0.05], 0.0, 0.6267943641604508, adj= 67.66802947739926  );
    atoms["C_d" ] = makeatom("C" , 6, 2, 4, 12.011,   4.0,   0,   [:s, :p, :d], -10.78129974, [-1.002529158795344; -0.3845826657687871; 0.05],0.0, 0.7978153170555154,adj= 60.71565982622933  );
    atoms["N_d" ] = makeatom("N" , 7, 2, 5, 14.007,   5.0,   0,   [:s, :p, :d], -19.47573456, [-1.3549552380156213; -0.5173699040224977; 0.05],0.0, 0.8874159249179556 , adj= 54.389500548598136  );
    atoms["O_d" ] = makeatom("O" , 8, 2, 6, 15.999,   6.0,   0,   [:s, :p, :d], -31.71324769, [-1.7468725107388499; -0.6600555484901767; 0.05],0.0, 1.015623980577927,adj= 49.061029052201036  );
    atoms["F_d" ] = makeatom("F" , 9, 2, 7, 18.998,   7.0,   0,   [:s, :p, :d], -48.22248233, [-2.184369166433599; -0.8122577461996937; 0.05],0.0, 1.1398170857521994,adj= 44.54705910604037  );
    atoms["Al_d"] = makeatom("Al", 13, 3, 3, 26.98,  3.0,   0,   [:s, :p, :d],  -6.37406439, [-0.5669635647556415; -0.19757116543256825; 0.05], 0.0, 0.47560947595799796 , adj= 120.35327370972065  );
    atoms["Si_d"] = makeatom("Si", 14, 3, 4, 28.09,  4.0,   0,   [:s, :p, :d], -9.15243181, [-0.7914081531895271; -0.2981042131884771; 0.05],0.0, 0.5597222300055447, adj= 124.00766606556414  );
    atoms["P_d" ] = makeatom("P" , 15, 3, 5, 30.97,  5.0,   0,   [:s, :p, :d], -15.07327848, [-1.0228315148577285; -0.40267982198157176; 0.05], 0.0, 0.6422064900805308 , adj= 120.26459710508915  );
    atoms["S_d" ] = makeatom("S" , 16, 3, 6, 32.06,  6.0,   0,   [:s, :p, :d], -23.70562681, [-1.2637320631089843; -0.5128774201992288; 0.05], 0.0, 0.7079117127047072, adj= 114.82659137149781  );
    atoms["Cl_d"] = makeatom("Cl", 17, 3, 7, 35.45,  7.0,   0,   [:s, :p, :d], -33.06905525, [-1.5153800694575343; -0.6292474932451964; 0.05], 0.0, 0.7476423278619385, adj= 155.47482542721207  );
    atoms["Ge_d"] = makeatom("Ge", 32, 4, 4.0,  72.63,   4.0,    10,   [:s, :p, :d], -212.84384318, [-0.8630308785518821; -0.2863400516133921; 0.05], 0.0, 0.5643618220580027, adj=155.47482542721207, semicore_orbitals=[:d], semicore_eigs=[-2.1543465449641492]  );
    atoms["As_d"] = makeatom("As", 33, 4, 5.0,  74.92,   5.0,    0,   [:s, :p, :d], -39.62501000, [-1.0658886642175263; -0.38183099697945; 0.05], 0.0, 0.6154159430456763 , adj=  132.97989766577297  );
    atoms["Se_d"] = makeatom("Se", 34, 4, 6.0,  78.97,   6.0,    0,   [:s, :p, :d],  -43.03651432, [-1.270876795104006; -0.47841368921375826; 0.05], 0.0, 0.6730297191236639, adj= 135.35647665899776  );
    atoms["Br_d"] = makeatom("Br", 35, 4, 7.0,  79.90,   7.0,    0,   [:s, :p, :d], -40.38769338, [-1.47972675653346; -0.5773272761649881; 0.05], 0.0, 0.7295557399434305 , adj= 136.30744895199706  );
    atoms["Sn_d"] = makeatom("Sn", 50, 5, 4.0,  118.71,   4.0,   10,   [:s, :p, :d], -158.61463088, [-0.7787699153070409; -0.2724932157659414; 0.05], 0.0, 0.5456114086239429, adj= 215.90551155126934  , semicore_orbitals=[:d], semicore_eigs=[-1.8945415352873742]);
    atoms["Sb_d"] = makeatom("Sb", 51, 5, 5.0,  121.76,   5.0,   10,   [:s, :p, :d], -184.28336729, [-0.9463940599724311; -0.3560272625722634; 0.05], 0.0, 0.5973377632845377 , adj= 221.4015210271326  , semicore_orbitals=[:d], semicore_eigs=[-2.4589362341980237]);
    atoms["Te_d"] = makeatom("Te", 52, 5, 6.0,  127.60,   6.0,    0,   [:s, :p, :d], -26.19141773, [-1.1135922760380588; -0.4385750533606238; 0.05], 0.0, 0.6243066065034633 , adj= 179.8013401206219  );
    atoms["I_d" ] = makeatom("I" , 53, 5, 7.0,  126.90,   7.0,    0,   [:s, :p, :d], -65.01087258 , [-1.2823860828584828; -0.5216184788129354; 0.05], 0.0, 0.6675305924093194, adj= 184.1984935967065  );
    atoms["Pb_d"] = makeatom("Pb", 82, 6, 4.0,  207.2,   4.0,   10,   [:s, :p, :d], -163.61141039, [-0.8890752912951138; -0.2585282858379871; 0.05], 0.0, 0.5488505326347655, adj= 241.50052083965022  , semicore_orbitals=[:d], semicore_eigs=[-1.5554213660160645] );
    atoms["Bi_d"] = makeatom("Bi", 83, 6, 5.0,  209.0,   5.0,   10,   [:s, :p, :d], -184.68674216, [-1.0683762029989725; -0.3379477237651605; 0.05], 0.0, 0.597524848429184 , adj= 252.9191136481541  , semicore_orbitals=[:d], semicore_eigs=[-1.9889042834093047]);

    for key in keys(atoms)
        atoms[Symbol(key)] = atoms[key]
    end
    
end

load_atoms()

zatoms = Dict()

for key in keys(atoms)
    println("key $key typeof ", typeof(key))
    if typeof(key) != String
        continue
    end
    if occursin("_d", key)
        continue
    end
    zatoms[atoms[key].Z] = key
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

        cutoff2X    = (rad1 + rad2) / 2.0 * 6.0
        #        cutoff_onX  = (rad1 + rad2) / 2.0 * 6.0
                cutoff_onX  = (rad1 + rad2) / 2.0 * 6.0

        cutoff2X   = max(min(cutoff2X,   19.01), 13.51) #2body
        cutoff_onX = max(min(cutoff_onX, 19.01), 13.51) #onsite

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
                                  old.orbitals, old.total_energy, [old.eigs[o] for o in old.orbitals],     0.0 ,U, adj=old.voladj, semicore_orbitals=old.semicore_orbitals, semicore_eigs=[old.semicore_eigs[o] for o in old.semicore_orbitals], fullU=old.fullU, Uarr=old.Uarr )
    
    atoms[Symbol(A)] = atoms[ String(A) ]
    atoms[String(A)] = atoms[ Symbol(A) ]

    uniform_charge_interaction[String(A)] = uni
    uniform_charge_interaction[Symbol(A)] = uni

end

for key in keys(data_temp)
    prepare_atom(key, data_temp[key][1], data_temp[key][2])
end


function reload()
    load_atoms()
    for key in keys(data_temp)
        prepare_atom(key, data_temp[key][1], data_temp[key][2])
    end
    
end

end #ends module
 
